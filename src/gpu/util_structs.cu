#include "headers/util_structs.cuh"

#ifdef C_CU_UF

cu_union_find* cu_uf_construct(int n) {
    // Alloc host mem & device mem to copy over
    uint32_t* host_h = new uint32_t[n];
    uint32_t* host_p = new uint32_t[n];
    uint32_t* dev_h;
    uint32_t* dev_p;
    cudaMalloc(&dev_h, sizeof(uint32_t) * n);
    cudaMalloc(&dev_p, sizeof(uint32_t) * n);

    // Init values for union find
    for(int i = 0; i < n; ++i) {
        host_h[i] = 1;
        host_p[i] = i;
    }

    // Copy mem over
    cudaMemcpy(dev_h, host_h, sizeof(uint32_t) * n, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_p, host_p, sizeof(uint32_t) * n, cudaMemcpyHostToDevice);

    // Create host side struct
    cu_union_find host_uf;
    host_uf.h = dev_h;
    host_uf.p = dev_p;
    host_uf.LEN = n;

    // Copy host uf data over to device uf
    cu_union_find* dev_uf;
    cudaMalloc(&dev_uf, sizeof(cu_union_find));
    cudaMemcpy(dev_uf, &host_uf, sizeof(cu_union_find), cudaMemcpyHostToDevice);

    // Cleanup and return
    delete[] host_h;
    delete[] host_p; 
    return dev_uf;
}

__device__ int cu_uf_find(cu_union_find* UF, int x) {
    // Avoid recursion, do iteration even if it kinda sucks
    int root = x;
    while(true) {
        int parent = UF->p[root];
        if(parent == root) {
            break;
        }
        root = parent;
    }

    // Thread safety
    while(UF->p[x] != root) {
        int parent = UF->p[x];
        atomicCAS(&UF->p[x], parent, root);
        x = parent;
    }

    // Path compressed
    return root;
}

/*
 *  Modified join-by-height -> h[root] = count of all nodes in tree `NEED TO FIX`
 *  @param UF `cu_union_find*`
 *  @param x `int`
 *  @param y `int`
 *  @return `void`
 */
__device__ void cu_uf_join(cu_union_find* UF, int x, int y) {
    int px = cu_uf_find(UF, x);
    int py = cu_uf_find(UF, y);
    if(px == py) {
        return;
    }
    
    // Need to use atomics for thread safety
    while(true) {
        px = cu_uf_find(UF, px);
        py = cu_uf_find(UF, py);
        
        // Joined
        if(px == py) {
            return;
        }
        
        // Join logic, favor x's parent as root in tiebreaker; favor lower index as root
        if(px < py) {
            if(atomicCAS(&UF->p[py], py, px) == py) {
                // Equal height
                if(UF->h[px] == UF->h[py]) {
                    atomicAdd(&UF->h[px], 1);
                }
                return;
            }
        } else {
            if(atomicCAS(&UF->p[px], px, py) == px) {
                // Equal height
                if(UF->h[px] == UF->h[py]) {
                    atomicAdd(&UF->h[py], 1);
                }
                return;
            }
        }

    }
}

__device__ bool cu_uf_con(cu_union_find* UF, int x, int y) {
    return cu_uf_find(UF, x) == cu_uf_find(UF, y);
}

#endif

#ifdef C_KMER_HASH_TABLE

__device__ uint64_t kh_hash(uint64_t key) {
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