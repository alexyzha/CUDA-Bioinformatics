#include "headers/fx_util.h"

void fq_to_file(const std::vector<fq_read*>& reads, std::string file_path) {
    std::ofstream file(file_path);
    for(auto& read : reads) {
        read->to_file(file);
    }
}

std::vector<fq_read*> filter_fq(const std::vector<fq_read*>& reads, char FILTER_BY, char THRESH, double PERC) {
    // Return object
    std::vector<fq_read*> filtered_reads;
    
    // Main processing loop
    for(auto& read : reads) {
        switch(FILTER_BY) {
            // AVERAGE_DISCARD_WHOLE: average quality < THRESH = discard
            case AVERAGE_DISCARD_WHOLE: {
                double average = static_cast<double>(sum(read->get_quality(), 0)) / read->size();
                if(average >= static_cast<double>(THRESH)) {
                    filtered_reads.push_back(
                        new fq_read(read->get_id(),
                                    read->size(),
                                    read->get_seq(),
                                    read->get_quality(),
                                    read->get_metadata())
                    );
                }
                break;
            }

            // SINGLE_DISCARD_WHOLE: single read quality < THRESH = discard
            case SINGLE_DISCARD_WHOLE: {
                bool discard = false;
                for(size_t i = 0; i < read->size(); ++i) {
                    char quality = static_cast<char>(((*read)[i] >> 8));
                    if(quality < THRESH) {
                        discard = true;
                        break;
                    }
                }
                if(!discard) {
                    filtered_reads.push_back(
                        new fq_read(read->get_id(),
                                    read->size(),
                                    read->get_seq(),
                                    read->get_quality(),
                                    read->get_metadata())
                    );
                }
                break;
            }

            // SLIDING_WINDOW: window of size static_cast<int>(PERC) average quality < THRESH = trim
            case SLIDING_WINDOW: {
                size_t WINDOW_SIZE = std::min(static_cast<size_t>(PERC), read->size());
                int sum = 0;
                int trim_index = read->size();

                // Set up the first window
                for(size_t i = 0; i < WINDOW_SIZE; ++i) {
                    char quality = static_cast<char>(((*read)[i] >> 8));
                    sum += quality;
                }

                // Validate first window, iff valid, try to search for largest trim_index
                if(sum < (THRESH * WINDOW_SIZE)) {
                    trim_index = 0;
                } else {
                    for(size_t i = WINDOW_SIZE; i < read->size(); ++i) {
                        sum -= static_cast<char>(((*read)[i - WINDOW_SIZE] >> 8));
                        sum += static_cast<char>(((*read)[i] >> 8));
                        if(sum < (THRESH * WINDOW_SIZE)) {
                            trim_index = i;
                            break;
                        }
                    }
                }

                // Trim and push new read
                if(trim_index) {
                    filtered_reads.push_back(
                        new fq_read(read->get_id(),
                                    trim_index,
                                    read->get_seq().substr(0, trim_index),
                                    read->get_quality().substr(0, trim_index),
                                    read->get_metadata())
                    );
                }
                break;
            }

            // PROPORTION_DISCARD_WHOLE: (percentage of read quality < THRESH) > PERC = discard
            case PROPORTION_DISCARD_WHOLE: {
                size_t ct = 0;
                for(size_t i = 0; i < read->size(); ++i) {
                    char quality = static_cast<char>(((*read)[i] >> 8));
                    if(quality >= THRESH) {
                        ++ct;
                    }
                }
                if(static_cast<double>(ct) / read->size() >= PERC) {
                    filtered_reads.push_back(
                        new fq_read(read->get_id(),
                                    read->size(),
                                    read->get_seq(),
                                    read->get_quality(),
                                    read->get_metadata())
                    );
                }
                break;
            }
        }
    }
    return filtered_reads;
}

std::vector<double> gc_per_read(const std::vector<fq_read*>& reads) {
    std::vector<double> percs;
    for(auto& read : reads) {
        int count = 0;
        for(size_t i = 0; i < read->size(); ++i) {
            // static_cast<char> preserves lowest 8 bits of the uint16_t bitmask, which is where the nucleotide is
            char base = static_cast<char>((*read)[i]);

            // Assumes no lowercase c/g
            if(base == 'C' || base == 'G') {
                ++count;
            }
        }

        // Cast to double before division to avoid flooring division
        percs.push_back(static_cast<double>(count) / read->size());
    }
    return percs;
}

double gc_global(const std::vector<fq_read*>& reads) {
    // Genomes are large. EX: 3B human bases > INT_MAX, which is 2.1B. Therefore, use 64-bit int
    uint64_t count = 0;
    uint64_t total_bases = 0;
    for(auto& read : reads) {
        total_bases += read->size();
        for(size_t i = 0; i < read->size(); ++i) {
            char base = static_cast<char>((*read)[i]);

            // Case sensitive, doesn't capture lowercase c/g
            if(base == 'C' || base == 'G') {
                ++count;
            }
        }
    }
    
    // Cast to double before division to prevent flooring
    return static_cast<double>(count) / total_bases;
}

std::unordered_map<uint64_t, uint64_t> count_kmer(const std::vector<fq_read*>& reads, size_t k) {
    // Using 2 bits per nucleotide, max nucleotides in a 64-bit int = 64/2 = 32
    if(k > 32 || k <= 0) {
        throw std::invalid_argument("K CANNOT BE GREATER THAN 32 OR LESS THAN 1");
    }

    // Rolling hash for O(1) kmer key creation
    uint64_t hash = 0;
    
    // Mask is a string of 1 bits from 0 (lsb) to the k*2th bit. EX: k = 2, mask = 0xf
    uint64_t mask = (1 << ((k * 2) % 64)) - 1;
    std::unordered_map<uint64_t, uint64_t> kmers;
    for(auto& read : reads) {
        if(read->size() < k) {
            continue;
        }
        for(size_t i = 0; i < read->size(); ++i) {
            hash <<= 2;

            // base_to_bit returns 0x80 in some cases. For ACGT bases, mask out non-2-bit values
            hash |= (static_cast<char>(0b11) & base_to_bit(static_cast<char>((*read)[i])));

            // From index 0 to k - 1, hash is not a yet a full kmer
            if(i >= k - 1) {
                ++kmers[hash & mask];
            }
        }
    }
    return kmers;
}

std::unordered_map<uint64_t, std::unordered_set<int>> index_kmer(const std::vector<fq_read*>& reads, size_t k) {
    // Using 2 bits per nucleotide, max nucleotides in a 64-bit int = 64/2 = 32
    if(k > 32 || k <= 0) {
        throw std::invalid_argument("K CANNOT BE GREATER THAN 32 OR LESS THAN 1");
    }

    // Rolling hash for O(1) kmer key creation
    uint64_t hash = 0;

    // Mask is a string of 1 bits from 0 (lsb) to the k*2th bit. EX: k = 2, mask = 0xf
    uint64_t mask = (1 << ((k * 2) % 64)) - 1;
    std::unordered_map<uint64_t, std::unordered_set<int>> kmer_map;
    for(size_t i = 0; i < reads.size(); ++i) {
        if(reads[i]->size() < k) {
            continue;
        }
        for(size_t j = 0; j < reads[i]->size(); ++j) {
            hash <<= 2;

            // base_to_bit returns 0x80 in some cases. For ACGT bases, mask out non-2-bit values
            hash |= (static_cast<char>(0b11) & base_to_bit(static_cast<char>((*reads[i])[j])));
            
            // From index 0 to k - 1, hash is not a yet a full kmer
            if(j >= k - 1) {
                kmer_map[hash & mask].insert(i);
            }
        }
    }
    return kmer_map;
}

alignment local_align(const std::string& ref, const std::string& read) {
    // 2D dp array
    std::vector<std::vector<int>> dp(ref.size() + 1, std::vector<int>(read.size() + 1, 0));
    int max_score = 0;
    int max_i = 0;
    int max_j = 0;
    
    // Run DP
    for(size_t i = 1; i <= ref.size(); ++i) {
        for(size_t j = 1; j <= read.size(); ++j) {
            int dgnl = dp[i - 1][j - 1] + (ref[i - 1] == read[j - 1] ? ALN_MATCH : ALN_MISMATCH);
            int up = dp[i - 1][j] + ALN_GAP;
            int left = dp[i][j - 1] + ALN_GAP;
            dp[i][j] = std::max({0, dgnl, up, left});
            if(dp[i][j] > max_score) {
                max_score = dp[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    // Traceback from DP results
    std::string* aligned_ref = new std::string("");
    std::string* aligned_read = new std::string("");
    int i = max_i;
    int j = max_j;

    // Main traceback loop
    while(i > 0 && j > 0 && dp[i][j] > 0) {
        int score = dp[i][j];
        int match = ref[i - 1] == read[j - 1] ? ALN_MATCH : ALN_MISMATCH;
        if(score == dp[i - 1][j - 1] + match) {
            aligned_ref->push_back(ref[i - 1]);
            aligned_read->push_back(read[j - 1]);
            --i;
            --j;
        } else if(score == dp[i - 1][j] + ALN_GAP) {
            aligned_ref->push_back(ref[i - 1]);
            aligned_read->push_back('-');
            --i;
        } else if(score == dp[i][j - 1] + ALN_GAP) {
            aligned_ref->push_back('-');
            aligned_read->push_back(read[j - 1]);
            --j;
        } else {
            // erm, what the sigma?
            break;
        }
    }

    // Traced sequences are built from high to low indices (backwards). Need to reverse before returning
    std::reverse(aligned_ref->begin(), aligned_ref->end());
    std::reverse(aligned_read->begin(), aligned_read->end());

    // Return custom alignment struct (headers/util_structs.h)
    return {max_score, max_i - 1, max_j - 1, aligned_ref, aligned_read};
}

std::vector<std::unordered_set<int>*> cluster_by_kmer(std::unordered_map<uint64_t, std::unordered_set<int>>& kmer_map, int READS, int THRESH) {
    // 2D "matrix" between all reads that have overlapping kmers
    std::unordered_map<int, std::unordered_map<int, int>> overlaps;
    
    // Iterate through all pairs of reads that are indexed with kmer x
    for(auto& [kmer, reads] : kmer_map) {
        for(auto i = reads.begin(); i != reads.end(); ++i) {
            for(auto j = std::next(i); j != reads.end(); ++j) {
                // Only count once to save space and time. Union find joining duplicate pairs is redundant
                ++overlaps[*i][*j];
            }
        }
    }

    // Use union find to join reads that have >= THRESH similar kmers
    union_find* uf = new union_find();
    for(auto& [i, row] : overlaps) {
        for(auto& [j, shared] : row) {
            if(shared >= THRESH) {
                uf->join(i, j);
            }
        }
    }

    // Convert to clusters using union find output
    std::vector<std::unordered_set<int>*> islands(READS, nullptr);
    for(int i = 0; i < READS; ++i) {
        int root = uf->find(i);
        if(root == i) {
            continue;
        }
        if(!islands[root]) {
            islands[root] = new std::unordered_set<int>();
        }
        islands[root]->insert(i);
    }
    return islands;
}