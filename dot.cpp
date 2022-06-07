#include <cstdio>

int main() {
    int v, e, x, y, z;
    freopen("input/g4.txt", "r", stdin);
    freopen("g4.dot", "w", stdout);
        printf("graph {\nnode [shape=point]\n");
    scanf("%d %d", &v, &e);
    for (int i = 0; i < e; i++) {
        scanf("%d %d %d", &x, &y, &z);
        if (x == 413 || y == 413 || x == 842 || y == 842)
            continue;
        printf("%d -- %d [color=%s]\n", x, y, z >= 0 ? "red" : "blue");
    }
    printf("}\n");
    return 0;
}
