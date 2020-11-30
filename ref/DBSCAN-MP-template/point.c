#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "point.h"
#include "constants.h"


double dist(point_t a, point_t b)
{
    return sqrt(pow(a.x - b.x, 2) +
            pow(a.y - b.y, 2) +
            pow(a.z - b.z, 2));
}

double dist2(double coX, double coY, double coZ,
             double a, double b, double c)
{
    return sqrt(pow(coX - a, 2) +
            pow(coY - b, 2) +
            pow(coZ - c, 2));
}

unsigned int parse_input(
    FILE *file,
    point_t **points,
    double *epsilon,
    unsigned int *minpts)
{
    unsigned int num_points, i = 0;
    fscanf(file, "%lf %u %u\n",
           epsilon, minpts, &num_points);
    point_t *p = (point_t *)
        calloc(num_points, sizeof(point_t));
    if (p == NULL) {
        perror("Failed to allocate points array");
        return 0;
    }
    while (i < num_points) {
          fscanf(file, "%lf %lf %lf\n",
                 &(p[i].x), &(p[i].y), &(p[i].z));
          p[i].cluster_id = UNCLASSIFIED;
          ++i;
    }
    *points = p;
    return num_points;
}

void print_points(
    point_t *points,
    unsigned int num_points)
{
    unsigned int i = 0;
    printf("Number of points: %u\n"
        " x     y     z     cluster_id\n"
        "-----------------------------\n"
        , num_points);
    while (i < num_points) {
          printf("%5.2lf %5.2lf %5.2lf: %d\n",
                 points[i].x,
                 points[i].y, points[i].z,
                 points[i].cluster_id);
          ++i;
    }
}
