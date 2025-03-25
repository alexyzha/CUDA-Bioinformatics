#include "headers/util_structs.h"

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