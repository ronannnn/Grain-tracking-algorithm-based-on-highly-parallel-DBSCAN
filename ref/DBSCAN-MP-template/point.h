#pragma once
typedef struct point_s point_t;
struct point_s {
    double x, y, z;
    int cluster_id;
};

typedef struct point_r point_rad;
struct point_r {
    double radius;
    int index;
};

double dist(point_t a, point_t b);
double dist2(double coX, double coY, double coZ, double a, double b, double c);
unsigned int parse_input(
    FILE *file,
    point_t **points,
    double *epsilon,
    unsigned int *minpts);
void print_points(
    point_t *points,
    unsigned int num_points);
