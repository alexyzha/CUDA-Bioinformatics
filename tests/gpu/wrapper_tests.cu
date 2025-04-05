#include "headers/test_util.h"
#include "headers/global_wrappers.cuh"
#include "headers/test_prototypes.cuh"
#include "../../cpu/headers/fx_util.h"

void CU_WRAPPER_TESTS(std::vector<TEST_RESULT*>& RESULTS) {
    srand(0xDEADBEEF);

    /*
     *  Phred: ! = min, K = max, 6/7 ~= mid
     */

    // Filter fastq by average discard whole
    RESULTS.push_back(TEST("CU_WRAPPER_TESTS", "FILTER_FQ_AVERAGE_DISCARD_WHOLE", [](){
        // Create test variables
        std::vector<fq_read*> reads;
        std::unordered_set<std::string> good_ids = {"10"};
        char THRESH = '6';
        for(int i = 0; i < 10; ++i) {
            std::string qual = "";
            std::string seq = "";
            int sum = 0;
            for(int i = 0; i < 100; ++i) {
                char q = (rand() % 42) + 33;
                sum += q;
                qual.push_back(q);
                seq.push_back('*');
            }

            // Check by id
            if(sum >= (THRESH * 100)) {
                good_ids.insert(std::to_string(i));
            }
            reads.push_back(new fq_read(
                std::to_string(i),
                100,
                seq,
                qual
            ));
        }

        // Add one guaranteed to pass
        reads.push_back(new fq_read(
            "10",
            100,
            std::string(100, 'A'),
            std::string(100, '8')
        ));

        // Run kernel wrapper
        std::vector<fq_read*> ret = cu_filter_fq(reads, AVERAGE_DISCARD_WHOLE, THRESH, 0, 0.0);

        // Check results against good_ids
        for(auto& read : ret) {
            EXPECT_TRUE(good_ids.count(read->get_id()));
            good_ids.erase(read->get_id());
        }
        EXPECT_TRUE(good_ids.empty());

        // Clean
        for(auto& read : reads) {
            delete read;
        }
    }));

    // Filter fastq by single discard whole
    RESULTS.push_back(TEST("CU_WRAPPER_TESTS", "FILTER_FQ_SINGLE_DISCARD_WHOLE", [](){
        std::vector<fq_read*> reads;
        std::unordered_set<std::string> good_ids = {"10"};
        char THRESH = '7';
        for(int i = 0; i < 10; ++i) {
            std::string qual = "";
            std::string seq = "";
            unsigned char min = 0;
            for(int i = 0; i < 100; ++i) {
                char q = (rand() % 42) + 33;
                min = (min < q ? min : q);
                qual.push_back(q);
                seq.push_back('*');
            }

            // Check by id
            if(min >= THRESH) {
                good_ids.insert(std::to_string(i));
            }
            reads.push_back(new fq_read(
                std::to_string(i),
                100,
                seq,
                qual
            ));
        }

        // Add one guaranteed to pass
        reads.push_back(new fq_read(
            "10",
            100,
            std::string(100, 'A'),
            std::string(100, '8')
        ));

        // Run kernel wrapper
        std::vector<fq_read*> ret = cu_filter_fq(reads, SINGLE_DISCARD_WHOLE, THRESH, 0, 0.0);

        // Check results against good_ids
        for(auto& read : ret) {
            EXPECT_TRUE(good_ids.count(read->get_id()));
            good_ids.erase(read->get_id());
        }
        EXPECT_TRUE(good_ids.empty());

        // Clean
        for(auto& read : reads) {
            delete read;
        }
    }));

    // Trim fastq by sliding window - trim at first window where avg(window) < thresh
    RESULTS.push_back(TEST("CU_WRAPPER_TESTS", "FILTER_FQ_SLIDING_WINDOW", [](){
        std::vector<fq_read*> reads;
        char THRESH = '6';
        std::vector<std::string> exp_seq = {"ABCDEFG", "ABCDE", "ABCDEFGHIJ"};
        std::vector<std::string> exp_qual = {"6666666", "77556", "7777777777"};

        // Add one with very obvious cut
        reads.push_back(new fq_read(
            "0",
            10,
            "ABCDEFGHIJ",
            "6666666555" 
        ));

        // Add one with less obvious cut
        reads.push_back(new fq_read(
            "1",
            10,
            "ABCDEFGHIJ",
            "7755667777"
        ));

        // Add one with no cut
        reads.push_back(new fq_read(
            "2",
            10,
            "ABCDEFGHIJ",
            "7777777777"
        ));

        // Run kernel wrapper
        std::vector<fq_read*> ret = cu_filter_fq(reads, SLIDING_WINDOW, THRESH, 5, 0.0);

        // Check results
        EXPECT_EQ(ret.size(), static_cast<size_t>(3));
        for(int i = 0; i < 3; ++i) {
            EXPECT_EQ(ret[i]->get_seq(), exp_seq[i]);
            EXPECT_EQ(ret[i]->get_quality(), exp_qual[i]);
        }

        // Clean
        for(auto& read : reads) {
            delete read;
        }
    }));

    // Filter fastq by proportion discard whole
    RESULTS.push_back(TEST("CU_WRAPPER_TESTS", "FILTER_FQ_PROPORTION_DISCARD_WHOLE", [](){
        std::vector<fq_read*> reads;
        char THRESH = '6';
        
        // Add one guaranteed to fail
        reads.push_back(new fq_read(
            "0",
            100,
            std::string(100, 'A'),
            std::string(100, '5')
        ));

        // Add one guaranteed to pass
        reads.push_back(new fq_read(
            "1",
            100,
            std::string(100, 'A'),
            std::string(100, '7')
        ));

        // Run kernel wrapper
        std::vector<fq_read*> ret = cu_filter_fq(reads, PROPORTION_DISCARD_WHOLE, THRESH, 0, 0.5);

        // Check results
        EXPECT_EQ(ret.size(), static_cast<size_t>(1));
        EXPECT_EQ(ret[0]->get_id(), std::string("1"));

        // Clean
        for(auto& read : reads) {
            delete read;
        }
    }));

    // Count kmer kernel test
    RESULTS.push_back(TEST("CU_WRAPPER_TESTS", "COUNT_KMERS", [](){
        // Create reads
        std::vector<fq_read*> reads;
        char bases[4] = {'A', 'C', 'G', 'T'};
        for(int i = 0; i < 100; ++i) {
            std::string seq(1000, '\0');
            std::string qual(1000, '?');
            for(int j = 0; j < 1000; ++j) {
                seq[j] = bases[rand() % 4];
            }
            reads.push_back(new fq_read(
                std::to_string(i),
                1000,
                seq,
                qual
            ));
        }

        // Expected from cpu code
        std::unordered_map<uint64_t, uint64_t> exp = count_kmer(reads, 7);

        // Run kernel
        std::unordered_map<uint64_t, uint64_t> ret = cu_count_kmers(reads, 7);
        
        // Check results
        EXPECT_EQ(ret.size(), exp.size());
        for(auto& [key, value] : ret) {
            // EXPECT_NE(exp.find(key), exp.end());
            // EXPECT_EQ(exp[key], value);
        }

        // Clean
        for(auto& read : reads) {
            delete read;
        }
    }));

    // Index kmer kernel test
    RESULTS.push_back(TEST("CU_WRAPPER_TESTS", "INDEX_KMERS", [](){
        // Create testing variables
        std::vector<fq_read*> reads;
        std::vector<std::string> seqs = {"ACG", "GAT", "CAT", "CAG", "ACG"};
        std::unordered_map<uint64_t, std::unordered_set<int>> exp = {
            {0b0001, {0, 4}}, {0b0110, {0, 4}}, {0b1000, {1}},
            {0b0011, {1, 2}}, {0b0100, {2, 3}}, {0b0010, {3}}
        };
        for(int i = 0; i < 5; ++i) {
            reads.push_back(new fq_read(
                std::to_string(i),
                3,
                seqs[i],
                "???"
            ));
        }

        // Run kernel
        std::unordered_map<uint64_t, std::unordered_set<int>> ret = cu_index_kmers(reads, 2);

        // Check results
        EXPECT_EQ(ret.size(), exp.size());
        for(auto& [key, value] : ret) {
            EXPECT_NE(exp.find(key), exp.end());
            EXPECT_EQ(exp[key].size(), value.size());
            for(auto& index : value) {
                EXPECT_TRUE(exp[key].count(index));
            }
        }

        // Clean
        for(auto& read : reads) {
            delete read;
        }
    }));

    // Kmer clustering kernel test
    RESULTS.push_back(TEST("CU_WRAPPER_TESTS", "CLUSTER_KMERS", [](){
        // Create testing variables
        std::vector<fq_read*> reads;
        std::vector<std::string> seqs = {"ACG", "GAT", "CAT", "CAG", "ACG"};
        std::unordered_map<uint64_t, std::unordered_set<int>> kmer_map;
        for(int i = 0; i < 5; ++i) {
            reads.push_back(new fq_read(
                "ID" + std::to_string(i),
                3,
                seqs[i],
                "???"
            ));
        }
        
        // Run clustering kernel
        std::vector<std::unordered_set<int>*> ret = cu_cluster_by_kmer(reads, 2, 2);

        // Check results
        for(int i = 0; i < ret.size(); ++i) {
            if(ret[i] != nullptr) {
                ret[i]->insert(i);
            }
        }
        std::sort(ret.begin(), ret.end(), [&](std::unordered_set<int>* a, std::unordered_set<int>* b){
            if(a && b) {
                return a->size() < b->size();
            }
            return (a ? true : false);
        });
        EXPECT_TRUE(ret[0] != nullptr);
        EXPECT_TRUE(ret[1] != nullptr);
        std::vector<std::vector<int>> exp_in = {{0, 4}, {1, 2, 3}};
        for(int i = 0; i < exp_in.size(); ++i) {
            EXPECT_EQ(ret[i]->size(), exp_in[i].size());
            for(auto& j : exp_in[i]) {
                EXPECT_TRUE(ret[i]->count(j));
            }
        }

        // Clean
        for(auto& read : reads) {
            delete read;
        }
    }));

    // Local align -> cigar test
    RESULTS.push_back(TEST("CU_WRAPPER_TESTS", "LOCAL_ALIGNMENT", [](){
        // Create testing variables;
        std::vector<fq_read*> reads;
        std::string ref = "GATCATG";
        std::vector<std::string> seqs = {
            "GATC",             // Full match
            "ATCA",             // Full match middle
            "ATCGT",            // Mismatch
            "GATATG",           // Gap in read
            "GATTCATG"          // Gap in ref
        };
        std::vector<std::string> exp_cig = {
            "4M",
            "4M",
            "3M1X1M",
            "3M1D3M",
            "2M1I5M"
        };
        std::vector<std::vector<int>> exp_nums {
            {3, 3, 8},          // {end ref, end read, score}
            {4, 3, 8},
            {5, 4, 7},
            {6, 5, 10},
            {6, 7, 12}
        };
        for(auto& seq : seqs) {
            reads.push_back(new fq_read(
                "[" + seq + "]",
                seq.size(),
                seq,
                std::string(seq.size(), '?')
            ));
        }

        // Run alignment kernel
        std::vector<cu_alignment*> ret = cu_local_align(ref, reads);

        // Check results
        EXPECT_EQ(ret.size(), reads.size());
        for(int i = 0; i < ret.size(); ++i) {
            EXPECT_EQ(*ret[i]->cigar, exp_cig[i]);
            EXPECT_EQ(ret[i]->end_ref, exp_nums[i][0]);
            EXPECT_EQ(ret[i]->end_read, exp_nums[i][1]);
            EXPECT_EQ(ret[i]->score, exp_nums[i][2]);
        }

        // Clean up
        for(auto& read : reads) {
            delete read;
        }
        for(auto& aln : ret) {
            delete aln->cigar;
            delete aln;
        }
    }));

}