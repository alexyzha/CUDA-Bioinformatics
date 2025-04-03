#include "headers/wrappers.cuh"

std::vector<fq_read*> cu_filter_fq(const std::vector<fq_read*>& READS, char FILTER_MODE, char THRESH, size_t K, double PROPORTION) {
    // Size variables
    int READ_LEN = READS.size();
    int FILTER_LENGTH = (READ_LEN / 64) + 1;
    int THREADS = MAX_THREADS;
    int BLOCKS = (READ_LEN / THREADS) + 1;
    
    // Host variables
    std::string temp = "";
    std::vector<uint32_t> h_offsets(READ_LEN);
    for(int i = 0; i < READ_LEN; ++i) {
        temp += READS[i]->get_quality();
        h_offsets[i] = temp.size() - 1;
    }
    int SEQ_LEN = temp.size() + 1;                  // temp + '\0'
    const char* h_allseq = temp.c_str();
    
    // Device variables
    char* d_allseq;
    uint32_t* d_offsets;    
    uint64_t* d_filter;

    // Allocate mem for device variables
    CUDA_CHECK(cudaMalloc(&d_allseq, sizeof(char) * SEQ_LEN));
    CUDA_CHECK(cudaMalloc(&d_offsets, sizeof(uint32_t) * READ_LEN));
    CUDA_CHECK(cudaMalloc(&d_filter, sizeof(uint64_t) * FILTER_LENGTH));
    
    // Copy mem host -> device/set mem
    CUDA_CHECK(cudaMemset(d_filter, 0, sizeof(uint64_t) * FILTER_LENGTH));
    CUDA_CHECK(cudaMemcpy(d_allseq, h_allseq, sizeof(char) * SEQ_LEN, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_offsets, h_offsets.data(), sizeof(uint32_t) * READ_LEN, cudaMemcpyHostToDevice));

    // Run kernels
    if(FILTER_MODE == SLIDING_WINDOW) {
        cu_filter_reads_sw <<<BLOCKS, THREADS>>> (
            d_allseq, 
            d_offsets, 
            READ_LEN, 
            K, 
            THRESH, 
            PROPORTION
        );
    } else {
        cu_filter_reads <<<BLOCKS, THREADS>>> (
            d_allseq, 
            d_offsets, 
            READ_LEN, 
            FILTER_MODE, 
            THRESH, 
            d_filter, 
            PROPORTION
        );
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy mem back
    std::vector<fq_read*> filtered_reads;
    if(FILTER_MODE != SLIDING_WINDOW) {
        // Copy mem back
        std::vector<uint64_t> ret_filter(FILTER_LENGTH, 0);
        CUDA_CHECK(cudaMemcpy(ret_filter.data(), d_filter, sizeof(uint64_t) * FILTER_LENGTH, cudaMemcpyDeviceToHost));

        // Filter and return
        for(int i = 0; i < READ_LEN; ++i) {
            int f_index = i / 64;
            int f_bit = i % 64;
            if(ret_filter[f_index] & ((ULL)1 << f_bit)) {
                filtered_reads.push_back(READS[i]);
            }
        }
    } else {
        // Copy mem back
        std::vector<char> temp_allseq(SEQ_LEN);
        CUDA_CHECK(cudaMemcpy(temp_allseq.data(), d_allseq, sizeof(char) * SEQ_LEN, cudaMemcpyDeviceToHost));
        const char* ret_allseq = temp_allseq.data();

        // Filter, trim, and return
        int seq_begin = 0;
        for(int i = 0; i < READ_LEN; ++i) {
            int og_size = READS[i]->size();
            
            // Get seq substring from seq buffer
            std::string cur_seq(ret_allseq + seq_begin);
            if(cur_seq.size() >= og_size) {
                filtered_reads.push_back(READS[i]);
            } else if(cur_seq.size()) {
                filtered_reads.push_back(new fq_read(
                    READS[i]->get_id(),
                    cur_seq.size(),
                    READS[i]->get_seq().substr(0, cur_seq.size()),
                    cur_seq,
                    READS[i]->get_metadata()
                ));
            }

            // Move to next seq begin
            seq_begin = h_offsets[i] + 1;
        }
    }
    
    // Cleanup and return
    CUDA_CHECK(cudaFree(d_filter));
    CUDA_CHECK(cudaFree(d_offsets));
    CUDA_CHECK(cudaFree(d_allseq));
    return filtered_reads;
}

std::unordered_map<uint64_t, uint64_t> cu_count_kmers(const std::vector<fq_read*>& READS, size_t K) {
    /*
     *  Given: 16 bytes per kh_pair, MAP_SIZE = 2 * total kmers, total bytes = 4^k * 2 * 16
     *  k = 13 requires INT_MAX bytes = 2gb vram/normal ram
     *  My GPU has 8gb vram...
     *  Note to self: consider setting K_MAX as 13 or 14
     *  -> Map needs size to be 2x the number of elements it's holding for better nocollide
     */
    if(K > 13) {
        throw std::runtime_error("K is too large, possible memory issues.");
    }
    
    // Size variables 
    int READ_LEN = READS.size(); 
    uint64_t MAP_SIZE = (1ULL << (2 * K)) * 2;
    int THREADS = MAX_THREADS;
    int BLOCKS = (READ_LEN / THREADS) + 1;
    
    // Host variables
    std::string temp = "";
    std::vector<uint32_t> h_offsets(READ_LEN);
    for(int i = 0; i < READ_LEN; ++i) {
        temp += READS[i]->get_seq();
        h_offsets[i] = temp.size() - 1;
    }
    int SEQ_LEN = temp.size() + 1;                  // temp + '\0'
    const char* h_allseq = temp.c_str();
    
    // Device variables
    char* d_allseq;
    uint32_t* d_offsets;
    kh_pair<uint64_t>* d_map = kh_construct<uint64_t>(MAP_SIZE);

    // Allocate mem for device variables
    CUDA_CHECK(cudaMalloc(&d_allseq, sizeof(char) * SEQ_LEN));
    CUDA_CHECK(cudaMalloc(&d_offsets, sizeof(uint32_t) * READ_LEN));

    // Copy mem host -> device
    CUDA_CHECK(cudaMemcpy(d_allseq, h_allseq, sizeof(char) * SEQ_LEN, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_offsets, h_offsets.data(), sizeof(uint32_t) * READ_LEN, cudaMemcpyHostToDevice));

    // Run kernel
    cu_kmer_count <<<BLOCKS, THREADS>>> (
        d_map, 
        d_allseq, 
        d_offsets, 
        K, 
        READ_LEN, 
        MAP_SIZE
    );
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy mem back
    std::vector<kh_pair<uint64_t>> h_map(MAP_SIZE);
    CUDA_CHECK(cudaMemcpy(h_map.data(), d_map, sizeof(kh_pair<uint64_t>) * MAP_SIZE, cudaMemcpyDeviceToHost));

    // Format data
    std::unordered_map<uint64_t, uint64_t> ret;
    ret.reserve(MAP_SIZE / 8);
    for(int i = 0; i < MAP_SIZE; ++i) {
        if(h_map[i].key == (ULL)EMPTY) {
            continue;
        }
        ret[h_map[i].key] = h_map[i].value;
    }

    // Clean and return
    CUDA_CHECK(cudaFree(d_offsets));
    CUDA_CHECK(cudaFree(d_allseq));
    CUDA_CHECK(cudaFree(d_map));
    return ret;
}

std::unordered_map<uint64_t, std::unordered_set<int>> cu_index_kmers(const std::vector<fq_read*>& READS, size_t K) {
    
    // Size variables
    int READ_LEN = READS.size(); 
    uint64_t MAP_SIZE = (1ULL << (2 * K)) * 2;
    int THREADS = MAX_THREADS;
    int BLOCKS = (READ_LEN / THREADS) + 1;

    // Host variables
    std::string temp = "";
    std::vector<uint32_t> h_offsets(READ_LEN);
    for(int i = 0; i < READ_LEN; ++i) {
        temp += READS[i]->get_seq();
        h_offsets[i] = temp.size() - 1;
    }
    int SEQ_LEN = temp.size() + 1;                  // temp + '\0'
    const char* h_allseq = temp.c_str();
    
    // Device variables
    char* d_allseq;
    uint32_t* d_offsets;
    kh_pair<uint32_t[MAP_MAX_INDICES + 1]>* d_map = kh_construct<uint32_t[MAP_MAX_INDICES + 1]>(MAP_SIZE);

    // Allocate mem for device variables
    CUDA_CHECK(cudaMalloc(&d_allseq, sizeof(char) * SEQ_LEN));
    CUDA_CHECK(cudaMalloc(&d_offsets, sizeof(uint32_t) * READ_LEN));

    // Copy mem host -> device
    CUDA_CHECK(cudaMemcpy(d_allseq, h_allseq, sizeof(char) * SEQ_LEN, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_offsets, h_offsets.data(), sizeof(uint32_t) * READ_LEN, cudaMemcpyHostToDevice));

    // Run kernel
    cu_kmer_index <<<BLOCKS, THREADS>>> (
        d_map,
        d_allseq,
        d_offsets,
        K,
        READ_LEN,
        MAP_SIZE
    );
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy mem back
    std::vector<kh_pair<uint32_t[MAP_MAX_INDICES + 1]>> h_map(MAP_SIZE);
    CUDA_CHECK(cudaMemcpy(h_map.data(), d_map, sizeof(kh_pair<uint32_t[MAP_MAX_INDICES + 1]>) * MAP_SIZE, cudaMemcpyDeviceToHost));

    // Format data, reserve space
    std::unordered_map<uint64_t, std::unordered_set<int>> ret;
    ret.reserve(MAP_SIZE / 4);
    for(int i = 0; i < MAP_SIZE; ++i) {
        if(h_map[i].key == (ULL)EMPTY) {
            continue;
        }

        // Insert up until MAP_MAX_INDICES matching indices
        for(int j = 1; j <= (h_map[i].value[0] - EMPTY); ++j) {
            ret[h_map[i].key].insert(h_map[i].value[j]);
        }
    }

    // Clean and return
    CUDA_CHECK(cudaFree(d_offsets));
    CUDA_CHECK(cudaFree(d_allseq));
    CUDA_CHECK(cudaFree(d_map));
    return ret;
}

std::vector<std::unordered_set<int>*> cu_cluster_by_kmer(const std::vector<fq_read*>& READS, size_t K, size_t THRESH) {
    /*
     *  Given: 130 * 8 = 1040 bytes per kh_pair, MAP_SIZE = 2 * total kmers, total bytes = 4^k * 2 * 1040
     *  k = 10 requires 2.18B bytes = 2gb vram/normal ram
     *  My GPU has 8gb vram...
     *  Note to self: consider setting K_MAX as 10 or 11
     *  -> Map needs size to be 2x the number of elements it's holding for better nocollide
     */
    if(K > 10) {
        throw std::runtime_error("K is too large, possible memory issues.");
    }

    // Size variables
    int READ_LEN = READS.size(); 
    uint64_t MAP_SIZE = (1ULL << (2 * K)) * 2;
    uint64_t MAX_EDGES = (ULL)READ_LEN * MAP_MAX_INDICES;
    uint64_t UF_SIZE = (ULL)READ_LEN;
    int THREADS = MAX_THREADS;
    int CLUSTER_MAP_SIZE = READ_LEN * 2;

    // Host variables
    std::string temp = "";
    std::vector<uint32_t> h_offsets(READ_LEN);
    for(int i = 0; i < READ_LEN; ++i) {
        temp += READS[i]->get_seq();
        h_offsets[i] = temp.size() - 1;
    }
    int SEQ_LEN = temp.size() + 1;                  // temp + '\0'
    const char* h_allseq = temp.c_str();
    
    // Device variables for kmer indexing
    char* d_allseq;
    uint32_t* d_offsets;
    kh_pair<uint32_t[MAP_MAX_INDICES + 1]>* d_map = kh_construct<uint32_t[MAP_MAX_INDICES + 1]>(MAP_SIZE);

    // Allocate mem for device variables
    CUDA_CHECK(cudaMalloc(&d_allseq, sizeof(char) * SEQ_LEN));
    CUDA_CHECK(cudaMalloc(&d_offsets, sizeof(uint32_t) * READ_LEN));

    // Copy mem host -> device
    CUDA_CHECK(cudaMemcpy(d_allseq, h_allseq, sizeof(char) * SEQ_LEN, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_offsets, h_offsets.data(), sizeof(uint32_t) * READ_LEN, cudaMemcpyHostToDevice));

    // Run kmer indexing kernel
    cu_kmer_index <<<(READ_LEN / THREADS) + 1, THREADS>>> (
        d_map,
        d_allseq,
        d_offsets,
        K,
        READ_LEN,
        MAP_SIZE
    );
    CUDA_CHECK(cudaDeviceSynchronize());

    // Device variables for kmer overlaps
    uint32_t* d_edgelist;
    uint32_t* d_edgecount_ptr;

    // Allocate mem for device variables
    CUDA_CHECK(cudaMalloc(&d_edgelist, sizeof(uint32_t) * MAX_EDGES));
    CUDA_CHECK(cudaMalloc(&d_edgecount_ptr, sizeof(uint32_t)));

    // Set values for device variables
    CUDA_CHECK(cudaMemset(d_edgelist, 0, sizeof(uint32_t) * MAX_EDGES));
    CUDA_CHECK(cudaMemset(d_edgecount_ptr, 0, sizeof(uint32_t)));

    // Run kmer overlap kernel
    cu_get_kmer_overlaps <<<(MAP_SIZE / THREADS) + 1, THREADS>>> (
        d_map,
        MAP_SIZE,
        d_edgelist,
        d_edgecount_ptr,
        MAX_EDGES
    );
    CUDA_CHECK(cudaDeviceSynchronize());

    // Device/host variables for union find
    cu_union_find* d_uf = cu_uf_construct(UF_SIZE);
    uint32_t h_edgecount = 0;
    CUDA_CHECK(cudaMemcpy(&h_edgecount, d_edgecount_ptr, sizeof(uint32_t), cudaMemcpyDeviceToHost));

    // Run UF kernel
    cu_get_uf <<<((h_edgecount / 2) / THREADS) + 1, THREADS>>> (
        d_uf,
        UF_SIZE,
        READ_LEN,
        d_edgelist,
        h_edgecount
    );
    CUDA_CHECK(cudaDeviceSynchronize());

    // Device variables for cluster kernel
    kh_pair<uint32_t[MAP_MAX_INDICES + 1]>* d_clusters = kh_construct<uint32_t[MAP_MAX_INDICES + 1]>(CLUSTER_MAP_SIZE);

    // Run cluster kernel
    cu_get_clusters <<<(UF_SIZE / THREADS) + 1, THREADS>>> (
        d_uf,
        d_clusters,
        UF_SIZE,
        CLUSTER_MAP_SIZE,
        THRESH
    );
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy mem back
    std::vector<kh_pair<uint32_t[MAP_MAX_INDICES + 1]>> h_ret_map(CLUSTER_MAP_SIZE);
    CUDA_CHECK(cudaMemcpy(h_ret_map.data(), d_clusters, sizeof(kh_pair<uint32_t[MAP_MAX_INDICES + 1]>) * CLUSTER_MAP_SIZE, cudaMemcpyDeviceToHost));

    // Format data
    std::vector<std::unordered_set<int>*> ret(READ_LEN, nullptr);

    // Enter data
    for(int i = 0; i < CLUSTER_MAP_SIZE; ++i) {
        if(h_ret_map[i].key == (ULL)EMPTY) {
            continue;
        }

        // Insert up until MAP_MAX_INDICES indices
        int root = h_ret_map[i].key;
        if(!ret[root]) {
            ret[root] = new std::unordered_set<int>();
        }
        for(int j = 1; j <= (h_ret_map[i].value[0] - EMPTY); ++j) {
            ret[root]->insert(h_ret_map[i].value[j]);
        }
    }

    // Clean and return
    CUDA_CHECK(cudaFree(d_clusters));
    cu_uf_destruct(d_uf);
    CUDA_CHECK(cudaFree(d_edgecount_ptr));
    CUDA_CHECK(cudaFree(d_edgelist));
    CUDA_CHECK(cudaFree(d_offsets));
    CUDA_CHECK(cudaFree(d_allseq));
    return ret;
}

std::vector<cu_alignment*> cu_local_align(const std::string& REF, const std::vector<fq_read*>& READS) {
    // Size consts
    int READ_LEN = READS.size();
    int THREADS = MAX_THREADS;
    int BLOCKS = (READ_LEN / MAX_THREADS) + 1;
    int DP_SIZE = (MAX_REF_LEN + 1) * (MAX_READ_LEN + 1) * READ_LEN;                // 200 * 170 * 4b * 25kseq ~= 3.4gb
    int CIGAR_SIZE = (MAX_CIGAR_LEN + 1) * READ_LEN;                                // 340 * 4b * 25kseq ~= 34m
    int ALIGN_LEN = READ_LEN * 3;

    // Host variables
    std::string temp = "";
    std::vector<uint32_t> h_offsets(READ_LEN);
    for(int i = 0; i < READ_LEN; ++i) {
        temp += READS[i]->get_seq();
        h_offsets[i] = temp.size() - 1;
    }
    int SEQ_LEN = temp.size() + 1;
    const char* h_allseq = temp.c_str();

    // Device variables
    char* d_allseq;
    char* d_cigarbuf;
    uint32_t* d_offsets;
    int* d_cache;
    int* d_align;

    // Allocate mem for device variables
    CUDA_CHECK(cudaMalloc(&d_allseq, sizeof(char) * SEQ_LEN));
    CUDA_CHECK(cudaMalloc(&d_cigarbuf, sizeof(char) * CIGAR_SIZE));
    CUDA_CHECK(cudaMalloc(&d_offsets, sizeof(uint32_t) * READ_LEN));
    CUDA_CHECK(cudaMalloc(&d_cache, sizeof(int) * DP_SIZE));
    CUDA_CHECK(cudaMalloc(&d_align, sizeof(int) * ALIGN_LEN));
    
    // Copy mem host -> device/set mem
    CUDA_CHECK(cudaMemset(d_cigarbuf, 0, sizeof(char) * CIGAR_SIZE));
    CUDA_CHECK(cudaMemset(d_cache, 0, sizeof(int) * DP_SIZE));
    CUDA_CHECK(cudaMemset(d_align, 0, sizeof(int) * ALIGN_LEN));
    CUDA_CHECK(cudaMemcpy(d_allseq, h_allseq, sizeof(char) * SEQ_LEN, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_offsets, h_offsets.data(), sizeof(uint32_t) * READ_LEN, cudaMemcpyHostToDevice));
    
    // Run local align kernel
    cu_local_alignment <<<BLOCKS, THREADS>>> (
        d_allseq,
        d_cigarbuf,
        d_cache,
        d_offsets,
        READ_LEN,
        d_align
    );
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy mem back
    std::vector<int> h_align(ALIGN_LEN);
    std::vector<char> h_cigarbuf(CIGAR_SIZE);
    CUDA_CHECK(cudaMemcpy(h_cigarbuf.data(), d_cigarbuf, sizeof(char) * CIGAR_SIZE, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_align.data(), d_align, sizeof(int) * ALIGN_LEN, cudaMemcpyDeviceToHost));

    // Format data
    std::vector<cu_alignment*> ret;
    for(int i = 0; i < READ_LEN; ++i) {
        int H_ALIGN_BEGIN = i * 3;
        int CIGAR_BEGIN = (MAX_CIGAR_LEN + 1) * i;
        h_cigarbuf[CIGAR_BEGIN + MAX_CIGAR_LEN] = '\0';
        ret.push_back(new cu_alignment{
            new std::string(h_cigarbuf.data() + CIGAR_BEGIN),
            h_align[H_ALIGN_BEGIN],
            h_align[H_ALIGN_BEGIN + 1],
            h_align[H_ALIGN_BEGIN + 2]
        });
    }

    // Clean and return
    CUDA_CHECK(cudaFree(d_align));
    CUDA_CHECK(cudaFree(d_cache));
    CUDA_CHECK(cudaFree(d_offsets));
    CUDA_CHECK(cudaFree(d_cigarbuf));
    CUDA_CHECK(cudaFree(d_allseq));
    return ret;
}