#include <stdlib.h>
#include <stdio.h>
#include "list.h"
#include "constants.h"

node_t *create_node(unsigned int index)
{
    node_t *n = (node_t *)calloc(1, sizeof(node_t));
    if (n == NULL)
        perror("Failed to allocate node.");
    else
    {
        n->index = index;
        n->next = NULL;
    }
    return n;
}

int append_at_end(
    unsigned int index,
    index_list *en)
{
    node_t *n = create_node(index);
    if (n == NULL)
    {
        free(en);
        return FAILURE;
    }
    if (en->head == NULL)
    {
        en->head = n;
        en->tail = n;
    }
    else
    {
        en->tail->next = n;
        en->tail = n;
    }
    ++(en->num_members);
    return SUCCESS;
}

void print_list(
    point_t *points,
    index_list *en)
{
    if (en)
    {
        node_t *h = en->head;
        while (h)
        {
            printf("(%lf, %lf, %lf) - rad: %f\n",
                   points[h->index].x,
                   points[h->index].y,
                   points[h->index].z,
                   h ->radius);
            h = h->next;
        }
    }
}
