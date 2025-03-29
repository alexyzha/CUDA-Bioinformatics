#include "headers/util_structs.h"

alignment::alignment() : score(0), end_ref(0), end_read(0), aligned_ref(nullptr), aligned_read(nullptr) {};

alignment::alignment(int s, int erf, int erd, std::string* arf, std::string* ard) : score(s), end_ref(erf), end_read(erd), aligned_ref(arf), aligned_read(ard) {};

int union_find::find(int x) {
    // Need to check null case of map: no parent means parent needs to be set to x
    if(!p.count(x)) {
        p[x] = x;

        // Height can technically be any default value, but is set to 1 for conventions
        h[x] = 1;
    }

    // Recursively do path compression
    if(x != p[x]) {
        return p[x] = find(p[x]);
    }
    return p[x];
}

void union_find::join(int x, int y) {
    int px = find(x);
    int py = find(y);

    // Pair already joined case
    if(px == py) {
        return;
    }

    // Prefer x as root to y in tiebreaker case
    if(h[px] >= h[py]) {
        h[px] += (h[px] == h[py]);                      // "++h[px]" and "p[py] = px" for tiebreaker case
        p[py] = px;                                     // (h[px] == h[py]) is a shorthand for "add 1 only when there is a height tie"
    } else {
        p[px] = py;
    }
}

bool union_find::con(int x, int y) {
    return find(x) == find(y);
}