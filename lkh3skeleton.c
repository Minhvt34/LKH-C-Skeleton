#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef struct Edge {
    int to; /* destination node index*/
    int alpha; /* placeholder for alpha-value*/
    double length; /* Euclidean distance */
} Edge;

typedef struct Node {
    int id; /* 1-based TSPLIBn id */
    double x, y; /* coordinates */
    int cand_count; /* current # of candidates */
    Edge *cands; /* dynamically sized array */
    int prev, next; /* prev/next node in tour */
    
} Node;

/* Global isntance data */
static Node *nodes = NULL;
static int N = 0; /* # of nodes */

/* Helper functions */
static double dist(int i, int j) {
    double dx = nodes[i].x - nodes[j].x;
    double dy = nodes[i].y - nodes[j].y;
    return round(sqrt(dx * dx + dy * dy)); /* TSPLIB-style Euclidean distance rounding*/
}

static void fatal(const char *msg) {
    fprintf(stderr, "Error: %s\n", msg);
    exit(EXIT_FAILURE);
}

/* TSPLIB (very small subset)*/
static void read_tsplib(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fatal("Failed to open .tsp file");
    }

    char buf[256];
    
    while (fgets(buf, sizeof(buf), file)) {
        if (sscanf(buf, "DIMENSION : %d", &N) == 1) {
            nodes = (Node *)malloc(N * sizeof(Node));
            if (!nodes) {
                fatal("Failed to allocate memory for nodes");
            }
        }
        if (strcmp(buf, "NODE_COORD_SECTION\n", 18) == 0) break;
    }

    if (!nodes) fatal("Dimension missing before NODE_COORD_SECTION");

    for(int i = 0; i < N; ++i) {
        int id;
        double x, y;
        if (fscanf(file, "%d %lf %lf", &id, &x, &y) != 3) {
            fatal("Error reading node coordinates");
        }

        nodes[i].id = id;
        nodes[i].x = x;
        nodes[i].y = y;
    }

    fclose(file);
}

/* Candidate set (k nearest by distance)*/

#define K_CANDIDATES 20

static int cmp_edges(const void *a, const void *b) {
    Edge *e1 = (Edge *)a;
    Edge *e2 = (Edge *)b;
    return e1->length - e2->length;
}

static void build_candidate_set(void) {
    for(int i = 0; i < N; ++i) {
        /*allocate maximum K_CANDIDATES edges*/
        Edge *cands = (Edge *)malloc(K_CANDIDATES * sizeof(Edge));
        if (!cands) {
            fatal("Failed to allocate memory for candidate edges");
        }

        /*fill cands with k nearest neighbors*/
        for(int j = 0; j < N; ++j) {
            if (j != i) {
                double length = dist(i, j);
                if (k < K_CANDIDATES) {
                    cands[k++] = (Edge){.to = j, .alpha = 0, .length = length};
                    if (k == K_CANDIDATES){
                        qsort(cands, k, sizeof(Edge), cmp_edges);
                    }
                } else if (length < cands[K_CANDIDATES - 1].length) {
                    cands[K_CANDIDATES - 1] = (Edge){.to = j, .alpha = 0, .length = length};
                    qsort(cands, K_CANDIDATES, sizeof(Edge), cmp_edges);
                }
            }
        }

        /*store sorted cands*/
        nodes[i].cand_count = K_CANDIDATES;
        nodes[i].cands = cands;
    }
}

/* Greedy nearest-neighbor tour*/
static double build_initial_tour(int *order){
    bool *visited = (bool *)calloc(N, sizeof(bool));
    if (!visited) {
        fatal("Failed to allocate memory for visited array");
    }

    int current = 0; visited[0] = true; order[0] = 0;
    double cost = 0.0;

    for(int k = 1; k < N; ++k) {
        double best_d = 1e100;
        int best_j = -1;

        for(int j = 0; j < N; ++j) {
            if (!visited[j]){
                double d = dist(current, j);
                if (d < best_d) {
                    best_d = d;
                    best_j = j;
                }
            }
        }

        visited[best_j] = true;
        order[k] = best_j;
        cost += best_d;
        current = best_j;
    }

    cost += dist(order[N-1], order[0]); /*close*/

    free(visited);
    return cost;
}

static void tour_from_order(int *order) {
    for(int i = 0; i < N; ++i) {
        int j = (i+1)%N;
        nodes[order[i]].next = order[j];
        nodes[order[j]].prev = order[i];
    }
}

/* Simple 2-opt improvement*/
static double tour_length(void){
    double length = 0.0;
    int start = 0, curr = start;
    do{
        int next = nodes[curr].next;
        length += dist(curr, next);
        curr = next;
    } while (curr != start);
    return length;
}

static bool two_opt_pass(void){
    bool improved = false;
    for(int i = 0; i < N;++i){
        int a = i;
        int b = nodes[a].next;
        for(int k = 0; k < nodes[a].cand_count; ++k){
            int c = nodes[a].cands[k].to;
            int d = nodes[c].next;
            if(c == a || c == b || d == a) continue;

            double delta = dist(a, b) + dist(c, d) - dist(a, c) - dist(b, d);
            if(delta > 1e-9){
                /*reverse segment [b .. c]*/
                int p = b, q;
                while(p != c){
                    q = nodes[p].next;
                    int tmp = nodes[p].next;
                    nodes[p].next = nodes[p].prev;
                    nodes[p].prev = tmp;
                    p = tmp;
                }

                /*reconnect*/
                nodes[a].next = c;
                nodes[b].prev = d;
                nodes[c].prev = a;
                nodes[d].next = b;
                improved = true;
            }
        }
    }

    return improved;
}

static void two_opt_local_search(void) {
    while (two_opt_pass()) {
        /*repeat until no improvement*/
    }
}

static void print_tour(int *order) {
    for(int i = 0; i < N; ++i) {
        printf("%d ", order[i]);
    }
    printf("\n");
}


int main() {
    
    if(argc < 2){
        fprintf(stderr, "Usage: %s instance.tsp [time] \n", argv[0]);
        return EXIT_FAILURE;
    }

    read_tsplib(argv[1]);
    build_candidate_set();

    int *order = (int *)malloc(N * sizeof(int));
    if (!order) {
        fatal("Failed to allocate memory for order array");
    }

    double length = build_initial_tour(order);
    tour_from_order(order);
    free(order);
    printf("Initial tour length: %.2f\n", length);

    two_opt_local_search();
    length = tour_length();
    printf("Optimized tour length: %.2f\n", length);
    print_tour(order);
    
    for(int i = 0; i < N; ++i){
        free(nodes[i].cands);
    }
    free(nodes);

    return 0;
}



