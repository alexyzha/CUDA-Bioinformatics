#include "headers/util_structs.h"

alignment::alignment() : score(0), end_ref(0), end_read(0), aligned_ref(nullptr), aligned_read(nullptr) {};

alignment::alignment(int s, int erf, int erd, std::string* arf, std::string* ard) : score(s), end_ref(erf), end_read(erd), aligned_ref(arf), aligned_read(ard) {};

int union_find::find(int x) {
    if(!p.count(x)) {
        p[x] = x;
        h[x] = 1;
    }
    if(x != p[x]) {
        return p[x] = find(p[x]);
    }
    return p[x];
}

void union_find::join(int x, int y) {
    int px = find(x);
    int py = find(y);
    if(px == py) {
        return;
    }
    if(h[px] >= h[py]) {
        h[px] += (h[px] == h[py]);
        p[py] = px;
    } else {
        p[px] = py;
    }
}

bool union_find::con(int x, int y) {
    return find(x) == find(y);
}