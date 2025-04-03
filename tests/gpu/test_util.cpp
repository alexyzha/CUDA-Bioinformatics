#include "headers/test_util.h"

#ifdef C_V_UF

v_union_find::v_union_find(int n) {
    h = new uint32_t[n];
    p = new uint32_t[n];
    for(int i = 0; i < n; ++i) {
        h[i] = 1;
        p[i] = i;
    }
}

v_union_find::~v_union_find() {
    delete[] h;
    delete[] p;
}

int v_union_find::find(int x) {
    if(x == p[x]) {
        return x;
    }
    h[p[x]] += h[x];
    h[x] = 0;
    return p[x] = find(p[x]);
}

void v_union_find::join(int x, int y) {
    int px = find(x);
    int py = find(y);
    if(px == py) {
        return;
    }
    if(px > py) {
        std::swap(px, py);
    }
    p[py] = px;
    h[px] += h[py];
    h[py] = 0;
}

bool v_union_find::con(int x, int y) {
    return find(x) == find(y);
}

#endif

void EXPECT_TRUE(bool ACT) {
    if(!ACT) {
        throw std::runtime_error("EXPECT_TRUE failed");
    }
}

void EXPECT_FALSE(bool ACT) {
    if(ACT) {
        throw std::runtime_error("EXPECT_FALSE failed");
    }
}