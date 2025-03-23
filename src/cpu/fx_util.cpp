#include "headers/fx_util.h"

std::vector<fq_read*> filter_fa(std::vector<fq_read*> reads, char FILTER_BY, char THRESH, double PERC) {
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

std::vector<double> gc_per_read(std::vector<fa_read*> reads) {
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

double gc_global(std::vector<fa_read*> reads) {
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

std::unordered_map<uint64_t, uint64_t> count_kmer(std::vector<fa_read*> reads, size_t k) {
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