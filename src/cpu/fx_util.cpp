#include "headers/fx_util.h"

std::vector<fq_read*> filter_fa(const std::vector<fq_read*>& reads, char FILTER_BY, char THRESH, double PERC) {
    std::vector<fq_read*> filtered_reads;
    for(auto& read : reads) {
        switch(FILTER_BY) {
            case AVERAGE_DISCARD_WHOLE: {
                double average = static_cast<double>(sum(read->get_quality(), PHRED_BEGIN)) / read->size();
                if(average >= static_cast<double>(FILTER_BY - PHRED_BEGIN)) {
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
            case SLIDING_WINDOW: {
                size_t WINDOW_SIZE = std::min(static_cast<size_t>(PERC), read->size());
                int sum = 0;
                int trim_index = read->size();
                for(size_t i = 0; i < WINDOW_SIZE; ++i) {
                    char quality = static_cast<char>(((*read)[i] >> 8));
                    sum += quality - PHRED_BEGIN;
                }
                if(static_cast<double>(sum / read->size()) < (FILTER_BY - PHRED_BEGIN)) {
                    trim_index = 0;
                } else {
                    for(size_t i = WINDOW_SIZE; i < read->size(); ++i) {
                        sum -= static_cast<char>(((*read)[i - WINDOW_SIZE] >> 8));
                        sum += static_cast<char>(((*read)[i] >> 8));
                        if(static_cast<double>(sum / read->size()) < (FILTER_BY - PHRED_BEGIN)) {
                            trim_index = i - WINDOW_SIZE + 1;
                            break;
                        }
                    }
                }
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
            case PROPORTION_DISCARD_WHOLE: {
                size_t ct = 0;
                for(size_t i = 0; i < read->size(); ++i) {
                    char quality = static_cast<char>(((*read)[i] >> 8));
                    if(quality < THRESH) {
                        ++ct;
                    }
                }
                if(static_cast<double>(ct / read->size()) >= PERC) {
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
}

std::vector<double> gc_per_read(const std::vector<fa_read*>& reads) {
    std::vector<double> percs;
    for(auto& read : reads) {
        int count = 0;
        for(size_t i = 0; i < read->size(); ++i) {
            char base = static_cast<char>((*read)[i]);
            if(base == 'C' || base == 'G') {
                ++count;
            }
        }
        percs.push_back(static_cast<double>(count) / read->size());
    }
    return percs;
}

double gc_global(const std::vector<fa_read*>& reads) {
    uint64_t count = 0;                             // 3bil human bases > INT_MAX
    uint64_t total_bases = 0;
    for(auto& read : reads) {
        total_bases += read->size();
        for(size_t i = 0; i < read->size(); ++i) {
            char base = static_cast<char>((*read)[i]);
            if(base == 'C' || base == 'G') {
                ++count;
            }
        }
    }
    return static_cast<double>(count / total_bases);
}

std::unordered_map<uint64_t, uint64_t> count_kmer(const std::vector<fa_read*>& reads, size_t k) {
    if(k > 32 || k <= 0) {
        throw std::invalid_argument("K CANNOT BE GREATER THAN 32 OR LESS THAN 1");
    }
    uint64_t hash = 0;
    uint64_t mask = (1 << ((k * 2) % 64)) - 1;
    std::unordered_map<uint64_t, uint64_t> kmers;
    for(auto& read : reads) {
        if(read->size() < k) {
            continue;
        }
        for(size_t i = 0; i < read->size(); ++i) {
            hash <<= 2;
            hash |= (static_cast<char>(0b11) & base_to_bit(static_cast<char>((*read)[i])));
            if(i >= k - 1) {
                ++kmers[hash & mask];
            }
        }
    }
    return kmers;
}

std::unordered_map<uint64_t, std::unordered_set<int>> index_kmer(const std::vector<fa_read*>& reads, size_t k) {
    if(k > 32 || k <= 0) {
        throw std::invalid_argument("K CANNOT BE GREATER THAN 32 OR LESS THAN 1");
    }
    uint64_t hash = 0;
    uint64_t mask = (1 << ((k * 2) % 64)) - 1;
    std::unordered_map<uint64_t, std::unordered_set<int>> kmer_map;
    for(size_t i = 0; i < reads.size(); ++i) {
        if(reads[i]->size() < k) {
            continue;
        }
        for(size_t j = 0; j < reads[i]->size(); ++j) {
            hash <<= 2;
            hash |= (static_cast<char>(0b11) & base_to_bit(static_cast<char>((*reads[i])[j])));
            if(i >= k - 1) {
                kmer_map[hash & mask].insert(i);
            }
        }
    }
    return kmer_map;
}

alignment local_align(const std::string& ref, const std::string& read) {
    std::vector<std::vector<int>> dp(ref.size() + 1, std::vector<int>(read.size() + 1, 0));
    int max_score = 0;
    int max_i = 0;
    int max_j = 0;
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
    std::string aligned_ref = "";
    std::string aligned_read = "";
    int i = max_i;
    int j = max_j;
    while(i > 0 && j > 0 && dp[i][j] > 0) {
        int score = dp[i][j];
        int match = ref[i - 1] == read[j - 1] ? ALN_MATCH : ALN_MISMATCH;
        if(score == dp[i - 1][j - 1] + match) {
            aligned_ref.push_back(ref[i - 1]);
            aligned_read.push_back(read[j - 1]);
            --i;
            --j;
        } else if(score == dp[i - 1][j] + ALN_GAP) {
            aligned_ref.push_back(ref[i - 1]);
            aligned_read.push_back('-');
            --i;
        } else if(score == dp[i][j - 1] + ALN_GAP) {
            aligned_ref.push_back('-');
            aligned_read.push_back(read[j - 1]);
            --j;
        } else {
            // erm, what the sigma?
            break;
        }
    }
    std::reverse(aligned_ref.begin(), aligned_ref.end());
    std::reverse(aligned_read.begin(), aligned_read.end());
    return {max_score, max_i - 1, max_j - 1, aligned_ref, aligned_read};
}

std::vector<std::unordered_set<int>*> cluster_by_kmer(std::unordered_map<uint64_t, std::unordered_set<int>>& kmer_map, int READS, int THRESH) {
    std::unordered_map<int, std::unordered_map<int, int>> overlaps;
    for(auto& [kmer, reads] : kmer_map) {
        for(auto i = reads.begin(); i != reads.end(); ++i) {
            for(auto j = std::next(i); j != reads.end(); ++j) {
                ++overlaps[*i][*j];
                ++overlaps[*j][*i];
            }
        }
    }
    union_find* uf;
    for(auto& [i, row] : overlaps) {
        for(auto& [j, shared] : row) {
            if(shared >= THRESH) {
                uf->join(i, j);
            }
        }
    }
    // Convert to islands
    std::vector<std::unordered_set<int>*> islands(READS);
    for(int i = 0; i < READS; ++i) {
        int root = uf->find(i);
        islands[root]->insert(i);
    }
    return islands;
}