#include "headers/xam_util.h"

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

    /*
    
        DO LATER
    
    */

}

std::vector<sam_read*> map_reads_to_ref(std::string ref, std::vector<fq_read*>& reads, size_t k) {
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
                best_pos_global = c + align.end_ref - align.aligned_ref.size() + 1;
            }
        }

        sam_read* sread = new sam_read{
            {},                         // No tags
            read->get_id(),
            ref,
            (best_score <= 0) ? "*" : make_cigar(best_align),
            "*",
            seq,
            read->get_quality(),
            (best_score <= 0) ? 0x4 : 0,
            (best_score <= 0) ? 0 : best_pos_global,
            0,
            0,
            255
        };

        mapped_reads.push_back(sread);
    }
    return mapped_reads;
}