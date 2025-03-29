#include "headers/xam_util.h"
#include "headers/util_structs.h"

void sam_read_to_file(std::ofstream& out, sam_read& read) {
    // sam_read is built for structured binding 
    auto [tags, qname, rname, cigar, rnext, seq, qual, flags, pos, posnext, tlen, mapq] = read;
    out << qname << "\t" 
        << flags << "\t" 
        << rname << "\t" 
        << pos << "\t" 
        << mapq << "\t"
        << cigar << "\t" 
        << rnext << "\t"
        << posnext << "\t"
        << tlen << "\t" 
        << seq << "\t" 
        << qual << "\t";

    // Trim when reading line will remove extra whitespace at end of line
    for(auto& tag : tags) {
        out << tag << "\t";
    }
    out << "\n";
}

void sam_to_file(sam_container& sam, std::string file_path) {
    // File handling
    std::ofstream file(file_path);
    if(!file.is_open()) {
        throw std::invalid_argument("FILE DOESN'T EXIST OR CANNOT OPEN FILE");
    }

    // Stream headers to file
    const auto& headers = sam.get_headers();
    for(auto& [tag, lines] : headers) {
        for(auto& line : lines) {
            file << "@" << tag << "\t" << line << "\n";
        }
    }

    // Stream reads to file
    const auto& reads = sam.get_reads();
    for(auto& read : reads) {
        sam_read_to_file(file, *read);
    }
    file.close();
}

std::string make_cigar(alignment& align) {
    std::string cigar = "";
    char prev = '\0';
    int count = 0;

    // Main processing loop
    for(size_t i = 0; i < (*align.aligned_ref).size(); ++i) {
        char rf = (*align.aligned_ref)[i];
        char rd = (*align.aligned_read)[i];
        char cur;
        if(rf == '-') {
            // Reference sequence has a gap (ref insert)
            cur = 'I';
        } else if(rd == '-') {
            // Read has a gap (ref delete)
            cur = 'D';
        } else if(rd != rf) {
            // Reference/read mismatch
            cur = 'X';
        } else {
            // Reference/read match
            cur = 'M';
        }

        // Add to chain or end chain
        if(cur == prev) {
            ++count;
        } else {
            if(count > 0) {
                cigar += std::to_string(count) + prev;
            }
            prev = cur;
            count = 1;
        }
    }

    // There may be processed but not pushed nucleotides
    if(count > 0) {
        cigar += std::to_string(count) + prev;
    }
    return cigar;
}

std::vector<sam_read*> map_reads_to_ref(std::string& ref, std::string ref_id, std::vector<fq_read*>& reads, size_t k) {
    // Using 2 bits per nucleotide, max nucleotides in a 64-bit int = 64/2 = 32
    if(k > 32 || k <= 0) {
        throw std::invalid_argument("K CANNOT BE GREATER THAN 32 OR LESS THAN 1");
    }
    std::unordered_map<uint64_t, std::vector<size_t>> ref_index;
    
    // Rolling hash for O(1) kmer key creation
    uint64_t hash = 0;

    // Mask is a string of 1 bits from 0 (lsb) to the k*2th bit. EX: k = 2, mask = 0xf
    uint64_t mask = (1 << ((k * 2) % 64)) - 1;
    for(size_t i = 0; i < ref.size(); ++i) {
        hash <<= 2;

        // base_to_bit returns 0x80 in some cases. For ACGT bases, mask out non-2-bit values
        hash |= (static_cast<char>(0b11) & base_to_bit(static_cast<char>(ref[i])));
        
        // From index 0 to k - 1, hash is not a yet a full kmer
        if(i >= k - 1) {
            ref_index[hash & mask].push_back(i - k + 1);
        }
    }

    // Read mapping loop
    std::vector<sam_read*> mapped_reads;
    for(size_t i = 0; i < reads.size(); ++i) {
        // Candidates for matching with ref based on kmer similarity
        fq_read* read = reads[i];
        const std::string& seq = read->get_seq();
        if(seq.size() < k) {
            continue;
        }
        std::unordered_set<size_t> cand;
        hash = 0;
        for(size_t j = 0; j < seq.size(); ++j) {
            hash <<= 2;

            // base_to_bit returns 0x80 in some cases. For ACGT bases, mask out non-2-bit values
            hash |= (static_cast<char>(0b11) & base_to_bit(static_cast<char>(seq[j])));
            
            // From index 0 to k - 1, hash is not a yet a full kmer
            if(j >= k - 1) {
                auto it = ref_index.find(hash & mask);
                if(it != ref_index.end()) {
                    for(auto pos : it->second) {
                        cand.insert(pos);
                    }
                }
            }
        }

        // Find best alignment from all candidates
        alignment best_align;
        int best_score = 0x80000000;
        size_t best_pos_global = 0;
        for(size_t c : cand) {
            if(c + seq.size() > ref.size()) {
                continue;
            }
            std::string window = ref.substr(c, seq.size() * 2);
            alignment align = local_align(window, seq);
            if(align.score > best_score) {
                best_score = align.score;
                best_align = align;
                best_pos_global = c + align.end_ref - (*align.aligned_ref).size() + 1;
            }
        }

        // Create sam read from results
        sam_read* sread = new sam_read{
            {},                         // No tags
            read->get_id(),
            ref_id,
            (best_score <= 0) ? "*" : make_cigar(best_align),
            "*",
            seq,
            read->get_quality(),
            static_cast<uint16_t>((best_score <= 0) ? 0x4 : 0),
            (best_score <= 0) ? 0 : best_pos_global,
            0,
            0,
            static_cast<char>(255)
        };
        mapped_reads.push_back(sread);
    }
    return mapped_reads;
}

void sam_to_vcf(std::string file_path, std::string& ref, std::string ref_id, std::vector<sam_read*>& reads, int CHROMO) {
    // File handling
    std::ofstream file(file_path);
    if(!file.is_open()) {
        throw std::invalid_argument("FILE DOESN'T EXIST OR UNABLE TO OPEN FILE");
    }

    // Stream headers to file
    file << "##fileformat=VCFv4.2\n";
    file << "##source=cubio\n";
    file << "##reference=" << ref_id << "\n";
    file << "##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    // Find (based on reference seq) and stream vcf to file
    std::unordered_set<std::string> vis;
    for(auto& read : reads) {
        if(read->flags & 0x4) {
            continue;
        }
        const std::string& seq = read->seq;
        size_t pos = read->pos;

        // Iterate over all read nucleotides
        for(size_t i = 0; i < seq.size(); ++i) {
            size_t ref_pos = pos + i - 1;
            if(ref_pos >= ref.size()) {
                break;
            }
            char base_ref = ref[ref_pos];
            char base_alt = seq[i];

            // Stream mismatch to file
            if(base_ref != base_alt && base_ref != 'N' && base_alt != 'N') {
                std::string key = std::to_string(CHROMO) + ":" + std::to_string(ref_pos + 1);
                if(vis.count(key)) {
                    continue;
                }
                file << CHROMO << "\t"
                     << ref_pos + 1 << "\t"
                     << "." << "\t"
                     << base_ref << "\t"
                     << base_alt << "\t"
                     << "255" << "\t"               // Placeholder quality
                     << "PASS" << "\t"
                     << "DP=1\n";
                vis.insert(key);
            }
        }
    }
    file.close();
}