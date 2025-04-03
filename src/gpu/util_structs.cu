#include "headers/util_structs.cuh"

#ifdef C_CU_UF

cu_union_find* cu_uf_construct(int n) {
    // Alloc host mem & device mem to copy over
    uint32_t* host_h = new uint32_t[n];
    uint32_t* host_p = new uint32_t[n];
    uint32_t* dev_h;
    uint32_t* dev_p;
    CUDA_CHECK(cudaMalloc(&dev_h, sizeof(uint32_t) * n));
    CUDA_CHECK(cudaMalloc(&dev_p, sizeof(uint32_t) * n));

    // Init values for union find
    for(int i = 0; i < n; ++i) {
        host_h[i] = 1;
        host_p[i] = i;
    }

    // Copy mem over
    CUDA_CHECK(cudaMemcpy(dev_h, host_h, sizeof(uint32_t) * n, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dev_p, host_p, sizeof(uint32_t) * n, cudaMemcpyHostToDevice));

    // Create host side struct
    cu_union_find host_uf;
    host_uf.h = dev_h;
    host_uf.p = dev_p;
    host_uf.LEN = n;

    // Copy host uf data over to device uf
    cu_union_find* dev_uf;
    CUDA_CHECK(cudaMalloc(&dev_uf, sizeof(cu_union_find)));
    CUDA_CHECK(cudaMemcpy(dev_uf, &host_uf, sizeof(cu_union_find), cudaMemcpyHostToDevice));

    // Cleanup and return
    delete[] host_h;
    delete[] host_p; 
    return dev_uf;
}

void cu_uf_destruct(cu_union_find* d_uf) {
    // Copy to host to access p/h pointers
    cu_union_find h_uf;
    CUDA_CHECK(cudaMemcpy(&h_uf, d_uf, sizeof(cu_union_find), cudaMemcpyDeviceToHost));

    // Free device arrays
    CUDA_CHECK(cudaFree(h_uf.h));
    CUDA_CHECK(cudaFree(h_uf.p));

    // Free device uf struct
    CUDA_CHECK(cudaFree(d_uf));
}

__device__ int __cu_uf_find(cu_union_find* UF, int x) {
    // Avoid recursion, do iteration even if it kinda sucks
    int total_h = 0;
    while(true) {
        int px = UF->p[x];
        if(px == x) {
            break;
        }
        int gpx = UF->p[px];
        int local_h = atomicExch(&UF->h[x], 0);
        total_h += local_h;

        // Compress
        atomicCAS(&UF->p[x], px, gpx);
        x = px;
    }

    // Path compressed + aggregate
    atomicAdd(&UF->h[x], total_h);
    return x;
}

/*
 *  Modified join-by-height -> h[root] = count of all nodes in tree
 *  @param UF `cu_union_find*`
 *  @param x `int`
 *  @param y `int`
 *  @return `void`
 */
__device__ void __cu_uf_join(cu_union_find* UF, int x, int y) {
    if(__cu_uf_con(UF, x, y)) {
        return;
    }
    
    // Need to use atomics for thread safety
    while(true) {
        int px = __cu_uf_find(UF, x);
        int py = __cu_uf_find(UF, y);
        
        // Joined
        if(px == py) {
            break;
        }
        
        // Favor smaller indexed root
        if(px > py) {
            int temp = px;
            px = py;
            py = temp;
        }

        // Link py to px
        if(atomicCAS(&UF->p[py], py, px) == py) {
            int moved_h = atomicExch(&UF->h[py], 0);
            atomicAdd(&UF->h[px], moved_h);
            __threadfence();
            break;
        }

    }
}

__device__ bool __cu_uf_con(cu_union_find* UF, int x, int y) {
    return __cu_uf_find(UF, x) == __cu_uf_find(UF, y);
}

#endif

#ifdef C_KMER_HASH_TABLE

__device__ uint64_t __kh_hash(uint64_t key) {
    key = (~key) + (key << 21);
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8);
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4);
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

#endif