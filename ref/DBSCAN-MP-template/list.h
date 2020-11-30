#include "point.h"
typedef struct node_s node_t;
struct node_s {
    unsigned int index;
    double radius;
    node_t *next;
};

typedef struct epsilon_neighbours_s index_list;
struct epsilon_neighbours_s {
    unsigned int num_members;
    node_t *head, *tail;
};

node_t *create_node(unsigned int index);
int append_at_end(
     unsigned int index,
     index_list *en);

void print_list(
    point_t *points,
    index_list *en);
