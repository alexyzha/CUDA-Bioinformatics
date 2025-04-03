#include "headers/test_prototypes.cuh"

void CU_UTIL_STRUCTS_TESTS(std::vector<TEST_RESULT*>& RESULTS) {
    srand(0x1EE7B00F);

    // Kmer hashtable constuctor test basic
    RESULTS.push_back(TEST("CU_UTIL_STRUCTS", "KMER_HASHTABLE_CONSTRUCTOR_BASIC", [](){
        // Create test variables
        kh_pair<uint64_t>* d_map = kh_construct<uint64_t>(100);
        kh_pair<uint64_t> h_map[100];
    
        // Copy device -> host
        CUDA_CHECK(cudaMemcpy(h_map, d_map, sizeof(h_map), cudaMemcpyDeviceToHost));
        for(int i = 0; i < 100; ++i) {
            EXPECT_EQ((uint64_t)EMPTY, h_map[i].key);
            EXPECT_EQ((uint64_t)EMPTY, h_map[i].value);
        }

        // Clean up
        CUDA_CHECK(cudaFree(d_map));
    }));
    
    // Kmer hashtable constuctor test array
    RESULTS.push_back(TEST("CU_UTIL_STRUCTS", "KMER_HASHTABLE_CONSTRUCTOR_ARRAY", [](){
        // Create test variables
        kh_pair<uint64_t[100]>* d_map = kh_construct<uint64_t[100]>(10);
        kh_pair<uint64_t[100]> h_map[10];
    
        // Copy device -> host
        CUDA_CHECK(cudaMemcpy(h_map, d_map, sizeof(h_map), cudaMemcpyDeviceToHost));
        for(int i = 0; i < 10; ++i) {
            EXPECT_EQ((uint64_t)EMPTY, h_map[i].key);
            for(int j = 0; j < 100; ++j) {
                EXPECT_EQ((uint64_t)EMPTY, h_map[i].value[j]);
            }
        }

        // Clean up
        CUDA_CHECK(cudaFree(d_map));
    }));
    
    // Hashing test
    RESULTS.push_back(TEST("CU_UTIL_STRUCTS", "KMER_HASHING", [](){
        // Create test variables
        uint64_t* d_key;
        uint64_t* d_value;
        std::vector<uint64_t> h_key(1024);
        for(auto& i : h_key) {
            i = static_cast<uint64_t>(rand()) ^ (static_cast<uint64_t>(rand()) << 32);
        }
        auto validate_hash = [](uint64_t key, uint64_t hash) -> bool {
            key = (~key) + (key << 21);
            key = key ^ (key >> 24);
            key = (key + (key << 3)) + (key << 8);
            key = key ^ (key >> 14);
            key = (key + (key << 2)) + (key << 4);
            key = key ^ (key >> 28);
            key = key + (key << 31);
            return (key == hash);
        };

        // Allocate mem, set/copy mem
        CUDA_CHECK(cudaMalloc(&d_key, sizeof(uint64_t) * 1024));
        CUDA_CHECK(cudaMalloc(&d_value, sizeof(uint64_t) * 1024));
        CUDA_CHECK(cudaMemset(d_value, 0, sizeof(uint64_t) * 1024));
        CUDA_CHECK(cudaMemcpy(d_key, h_key.data(), sizeof(uint64_t) * 1024, cudaMemcpyHostToDevice));

        // Run kernel
        cu_hashing <<<32, 32>>> (d_value, d_key, 1024);
        CUDA_CHECK(cudaDeviceSynchronize());

        // Copy mem back & validate results
        std::vector<uint64_t> h_value(1024);
        CUDA_CHECK(cudaMemcpy(h_value.data(), d_value, sizeof(uint64_t) * 1024, cudaMemcpyDeviceToHost));
        for(int i = 0; i < 1024; ++i) {
            EXPECT_TRUE(validate_hash(h_key[i], h_value[i]));
        }

        // Clean up
        CUDA_CHECK(cudaFree(d_value));
        CUDA_CHECK(cudaFree(d_key));
    }));

    // Union find find test
    RESULTS.push_back(TEST("CU_UTIL_STRUCTS", "UF_FIND", [](){
        // Consts
        int NODES = 200;
        
        // Create union find and index list
        v_union_find val_uf(NODES);
        for(int i = 0; i < 256; ++i) {
            val_uf.join(rand() % NODES, rand() % NODES);
        }

        std::vector<uint32_t> h_indexlist(NODES);
        std::iota(h_indexlist.begin(), h_indexlist.end(), 0);
        
        // Device variables
        cu_union_find* d_uf;
        uint32_t* d_roots;
        uint32_t* d_indexlist;
        
        // Allocate/copy/set mem
        CUDA_CHECK(cudaMalloc(&d_uf, sizeof(cu_union_find)));
        CUDA_CHECK(cudaMalloc(&d_roots, sizeof(uint32_t) * NODES));
        CUDA_CHECK(cudaMalloc(&d_indexlist, sizeof(uint32_t) * NODES));
        CUDA_CHECK(cudaMemset(d_roots, 0, sizeof(uint32_t) * NODES));
        CUDA_CHECK(cudaMemcpy(d_indexlist, h_indexlist.data(), sizeof(uint32_t) * NODES, cudaMemcpyHostToDevice));
        
        // Special mem handle for struct
        cu_union_find* temp_uf = new cu_union_find;
        CUDA_CHECK(cudaMalloc(&temp_uf->h, sizeof(uint32_t) * NODES));
        CUDA_CHECK(cudaMalloc(&temp_uf->p, sizeof(uint32_t) * NODES));
        CUDA_CHECK(cudaMemcpy(temp_uf->h, val_uf.h, sizeof(uint32_t) * NODES, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(temp_uf->p, val_uf.p, sizeof(uint32_t) * NODES, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_uf, temp_uf, sizeof(cu_union_find), cudaMemcpyHostToDevice));

        // Run kernel
        cu_uf_find <<<16, 16>>> (d_uf, d_indexlist, d_roots, NODES);
        CUDA_CHECK(cudaDeviceSynchronize());

        // Validate results
        std::vector<uint32_t> h_roots(NODES);
        CUDA_CHECK(cudaMemcpy(h_roots.data(), d_roots, sizeof(uint32_t) * NODES, cudaMemcpyDeviceToHost));
        for(int i = 0; i < NODES; ++i) {
            EXPECT_EQ(h_roots[i], static_cast<uint32_t>(val_uf.find(h_indexlist[i])));
        }
        
        // Clean up
        CUDA_CHECK(cudaFree(d_indexlist));
        CUDA_CHECK(cudaFree(d_roots));
        cu_uf_destruct(d_uf);
        free(temp_uf);
    }));

    // Union find con test
    RESULTS.push_back(TEST("CU_UTIL_STRUCTS", "UF_CON", [](){
        // Consts
        int NODES = 200;
        int EDGES = 256 * 2;
        int BOOLS = 256;

        // Create union find and edge list
        v_union_find val_uf(NODES);
        std::vector<uint32_t> h_edgelist;
        for(int i = 0; i < 256; ++i) {
            int first = rand() % NODES;
            int second = rand() % NODES;
            h_edgelist.push_back(first);
            h_edgelist.push_back(second);
            val_uf.join(first, second);
        }

        // Device variables
        cu_union_find* d_uf;
        uint32_t* d_edgelist;
        char* d_ret;

        // Allocate/copy/set mem
        CUDA_CHECK(cudaMalloc(&d_uf, sizeof(cu_union_find)));
        CUDA_CHECK(cudaMalloc(&d_ret, sizeof(char) * BOOLS));
        CUDA_CHECK(cudaMalloc(&d_edgelist, sizeof(uint32_t) * EDGES));
        CUDA_CHECK(cudaMemset(d_ret, 0, sizeof(char) * BOOLS));
        CUDA_CHECK(cudaMemcpy(d_edgelist, h_edgelist.data(), sizeof(uint32_t) * EDGES, cudaMemcpyHostToDevice));
        
        // Special mem handle for struct
        cu_union_find* temp_uf = new cu_union_find;
        CUDA_CHECK(cudaMalloc(&temp_uf->h, sizeof(uint32_t) * NODES));
        CUDA_CHECK(cudaMalloc(&temp_uf->p, sizeof(uint32_t) * NODES));
        CUDA_CHECK(cudaMemcpy(temp_uf->h, val_uf.h, sizeof(uint32_t) * NODES, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(temp_uf->p, val_uf.p, sizeof(uint32_t) * NODES, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_uf, temp_uf, sizeof(cu_union_find), cudaMemcpyHostToDevice));

        // Run kernel
        cu_uf_con <<<16, 16>>> (d_uf, d_edgelist, d_ret, EDGES, BOOLS);
        CUDA_CHECK(cudaDeviceSynchronize());

        // Validate results
        std::vector<char> h_ret(BOOLS);
        CUDA_CHECK(cudaMemcpy(h_ret.data(), d_ret, sizeof(char) * BOOLS, cudaMemcpyDeviceToHost));
        for(int i = 0; i < NODES; ++i) {
            int first = h_edgelist[i * 2];
            int second = h_edgelist[i * 2 + 1];
            EXPECT_EQ(static_cast<bool>(h_ret[i]), val_uf.con(first, second));
        }

        // Clean up
        CUDA_CHECK(cudaFree(d_edgelist));
        CUDA_CHECK(cudaFree(d_ret));
        cu_uf_destruct(d_uf);
        free(temp_uf);
    }));

    // Union find join test
    RESULTS.push_back(TEST("CU_UTIL_STRUCTS", "UF_JOIN", [](){
        // Consts
        int NODES = 200;
        int EDGES = 256 * 2;
        int NUM_RUNS = 10;
        
        // Create edge list & validation
        v_union_find val_uf(NODES);     
        std::vector<uint32_t> h_edgelist;
        for(int i = 0; i < 256; ++i) {
            int first = rand() % NODES;
            int second = rand() % NODES;
            h_edgelist.push_back(first);
            h_edgelist.push_back(second);
            val_uf.join(first, second);
        } 

        // Device variables
        cu_union_find* d_uf = cu_uf_construct(NODES);
        uint32_t* d_edgelist;
        
        // Allocate mem
        CUDA_CHECK(cudaMalloc(&d_edgelist, sizeof(uint32_t) * EDGES));
        CUDA_CHECK(cudaMemcpy(d_edgelist, h_edgelist.data(), sizeof(uint32_t) * EDGES, cudaMemcpyHostToDevice));

        // Run kernel
        for(int i = 0; i < NUM_RUNS; ++i) {
            cu_uf_join <<<16, 16>>> (d_uf, d_edgelist, EDGES);
            CUDA_CHECK(cudaDeviceSynchronize());
        }
        
        // Copy mem over
        cu_union_find h_uf;
        v_union_find ret_uf(NODES);
        CUDA_CHECK(cudaMemcpy(&h_uf, d_uf, sizeof(cu_union_find), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(ret_uf.h, h_uf.h, sizeof(uint32_t) * NODES, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(ret_uf.p, h_uf.p, sizeof(uint32_t) * NODES, cudaMemcpyDeviceToHost));

        // Validate
        for(int i = 0; i < NODES; ++i) {
            // Expect size of joined trees to be equal
            EXPECT_EQ(
                ret_uf.h[ret_uf.find(i)],
                val_uf.h[val_uf.find(i)]
            );
        }

        // Clean
        CUDA_CHECK(cudaFree(d_edgelist));
        cu_uf_destruct(d_uf);
    }));

}