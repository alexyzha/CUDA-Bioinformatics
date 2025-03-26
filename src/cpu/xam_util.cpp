#include "headers/xam_util.h"
#include "headers/util_structs.h"

void sam_read_to_file(std::ofstream& out, sam_read& read) {
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
    for(auto& tag : tags) {
        out << tag << "\t";
    }
    out << "\n";
}

void sam_to_file(sam_container& sam, std::string file_path) {
    std::ofstream file(file_path);
    if(!file.is_open()) {
        throw std::invalid_argument("FILE DOESN'T EXIST OR CANNOT OPEN FILE");
    }
    const auto& headers = sam.get_headers();
    for(auto& [tag, lines] : headers) {
        for(auto& line : lines) {
            file << "@" << tag << "\t" << line << "\n";
        }
    }
    const auto& reads = sam.get_reads();
    for(auto& read : reads) {
        sam_read_to_file(file, *read);
    }
    file.close();
}

std::string make_cigar(alignment& align) {
    std::string cigar;
    char prev = '\0';
    int count = 0;
    for(size_t i = 0; i < (*align.aligned_ref).size(); ++i) {
        char rf = (*align.aligned_ref)[i];
        char rd = (*align.aligned_read)[i];
        char cur;
        if(rf == '-') {
            cur = 'I';
        } else if(rd == '-') {
            cur = 'D';
        } else {
            cur = 'M';
        }
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
    if(count > 0) {
        cigar += std::to_string(count) + prev;
    }
    return cigar;
}

std::vector<sam_read*> map_reads_to_ref(std::string& ref, std::string ref_id, std::vector<fq_read*>& reads, size_t k) {
    if(k > 32 || k <= 0) {
        throw std::invalid_argument("K CANNOT BE GREATER THAN 32 OR LESS THAN 1");
    }
    std::unordered_map<uint64_t, std::vector<size_t>> ref_index;
    uint64_t hash = 0;
    uint64_t mask = (1 << ((k * 2) % 64)) - 1;
    for(size_t i = 0; i < ref.size(); ++i) {
        hash <<= 2;
        hash |= (static_cast<char>(0b11) & base_to_bit(static_cast<char>(ref[i])));
        if(i >= k - 1) {
            ref_index[hash & mask].push_back(i - k + 1);
        }
    }

    std::vector<sam_read*> mapped_reads;
    for(size_t i = 0; i < reads.size(); ++i) {
        fq_read* read = reads[i];
        const std::string& seq = read->get_seq();
        if(seq.size() < k) {
            continue;
        }
        std::unordered_set<size_t> cand;
        hash = 0;
        for(size_t j = 0; j < seq.size(); ++j) {
            hash <<= 2;
            hash |= (static_cast<char>(0b11) & base_to_bit(static_cast<char>(seq[j])));
            if(j >= k - 1) {
                auto it = ref_index.find(hash & mask);
                if(it != ref_index.end()) {
                    for(auto pos : it->second) {
                        cand.insert(pos);
                    }
                }
            }
        }

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
    std::ofstream file(file_path);
    if(!file.is_open()) {
        throw std::invalid_argument("FILE DOESN'T EXIST OR UNABLE TO OPEN FILE");
    }

    // Headers
    file << "##fileformat=VCFv4.2\n";
    file << "##source=cubio";
    file << "##reference=" << ref_id << "\n";
    file << "##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    std::unordered_set<std::string> vis;
    for(auto& read : reads) {
        if(read->flags & 0x4) {
            continue;
        }
        const std::string& seq = read->seq;
        size_t pos = read->pos;
        for(size_t i = 0; i < seq.size(); ++i) {
            size_t ref_pos = pos + i - 1;
            if(ref_pos >= ref.size()) {
                break;
            }
            char base_ref = ref[ref_pos];
            char base_alt = seq[i];
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