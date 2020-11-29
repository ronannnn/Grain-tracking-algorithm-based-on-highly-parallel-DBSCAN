/**
 * Author:   Ruonan Chen (ruonanch@usc.edu)
 * Date:     10/29/20
 * Filename: dbscan_pbc.c
 */
#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "math.h"

#define MAX_LENGTH_OF_LINE 1024 // max length of each line for file reader
#define SKIP_LINE 2             // skip line number when reading file
#define MIN_POINTS 4            // predefined minimum number of points(neighbors)
#define EPS 2.5                 // predefined radius
#define POINT_NUM 12800         // total number of points

#define DEBUG 1 // debug mode

typedef struct tagAtom {
  int id;
  int cluster_id; // when it is -1, it means the atom is labeled as NOISE
  double coord[3];
} Atom;

static char *delim = " "; // one or more space

static Atom atoms[POINT_NUM]; // atom array
static double dist[POINT_NUM][POINT_NUM];

static char *substring(const char *string, int position, int length) {
  char *pointer = calloc(length + 1, sizeof(char));
  for (int i = 0; i < length; i++) pointer[i] = string[i + position];
  return pointer;
}

// parse dimension of each atom
static Atom *parse(char *line, int id) {
  Atom *atom = calloc(1, sizeof(Atom));
  char *ptr;
  char *token = strtok_r(line, delim, &ptr);
  for (int i = 0; i < 3; i++) {
    if (token == NULL) return NULL; // exclude error lines
    atom->coord[i] = strtod(token, NULL);
    token = strtok_r(NULL, delim, &ptr);
  }
  atom->id = id;
  return atom;
}

static void fetch_data(FILE *stream) {
  char line[MAX_LENGTH_OF_LINE];
  // skip lines
  for (int i = 0; i < SKIP_LINE; i++) {
    fgets(line, MAX_LENGTH_OF_LINE, stream);
  }
  // read line by line
  for (int i = 0; i < POINT_NUM; i++) {
    fgets(line, MAX_LENGTH_OF_LINE, stream);
    Atom *atom = parse(substring(line, 0, (int)strlen(line) - 2), i + 1);
    if (atom != NULL) {
      atoms[i] = *atom;
    } else {
      i--;
    }
  }
}

static void print_atoms() {
  for (int i = 0; i < POINT_NUM; i++) {
    printf(
      "id: %d, cluster id: %d\n"
      "coords: %f, %f, %f\n",
      atoms[i].id, atoms[i].cluster_id,
      atoms[i].coord[0], atoms[i].coord[1], atoms[i].coord[2]
    );
  }
}

static void export_atoms() {
  FILE *output = fopen("atoms.txt", "w");
  fprintf(output, "%d", POINT_NUM);
  for (int i = 0; i < POINT_NUM; i++) {
    fprintf(output, "%f %f %f %d\n", atoms[i].coord[0], atoms[i].coord[1], atoms[i].coord[2], atoms[i].cluster_id);
  }
}

static void calculate_distances() {
  double xmax, xmin, ymax, ymin, zmax, zmin;
  // find the box size (man & min for each direction)
  xmax = atoms[0].coord[0];
  xmin = atoms[0].coord[0];
  ymax = atoms[0].coord[1];
  ymin = atoms[0].coord[1];
  zmax = atoms[0].coord[2];
  zmin = atoms[0].coord[2];

  for (int i = 1; i < POINT_NUM; i++) {
    if (xmax < atoms[i].coord[0]) xmax = atoms[i].coord[0];
    if (xmin > atoms[i].coord[0]) xmin = atoms[i].coord[0];
    if (ymax < atoms[i].coord[1]) ymax = atoms[i].coord[1];
    if (ymin > atoms[i].coord[1]) ymin = atoms[i].coord[1];
    if (zmax < atoms[i].coord[2]) zmax = atoms[i].coord[2];
    if (zmin > atoms[i].coord[2]) zmin = atoms[i].coord[2];
  }
  double xhalf = (xmax - xmin) / 2.0;
  double yhalf = (ymax - ymin) / 2.0;
  double zhalf = (zmax - zmin) / 2.0;
  double xlength = (xmax - xmin);
  double ylength = (ymax - ymin);
  double zlength = (zmax - zmin);

  for (int i = 0; i < POINT_NUM; i++) {
    for (int j = i + 1; j < POINT_NUM; j++) {
      double xdistance = abs(atoms[i].coord[0] - atoms[j].coord[0]);
      double ydistance = abs(atoms[i].coord[1] - atoms[j].coord[1]);
      double zdistance = abs(atoms[i].coord[2] - atoms[j].coord[2]);
      if (xdistance > xhalf) xdistance = xlength - xdistance;
      if (ydistance > yhalf) ydistance = ylength - ydistance;
      if (zdistance > zhalf) zdistance = zlength - zdistance;

      double distance = sqrt(pow(xdistance, 2) + pow(ydistance, 2) + pow(zdistance, 2));
      dist[i][j] = distance;
      dist[j][i] = distance;
    }
  }
}

/**
 * @return int[]{neighbor size, unclustered neighbor size}
 */
static void get_neighbors(int atom_id, int *neighbors, int *size, int *unclustered_size) {
  for (int i = 0; i < POINT_NUM; i++) {
    if (atom_id == i) continue;
    if (dist[atom_id][i] <= EPS) {
      neighbors[(*size)++] = i;
      if (atoms[i].cluster_id == 0) {
        (*unclustered_size)++;
      }
    }
  }
  /*for (int j = 0; j < *size; j++) {
    printf("%d, ", neighbors[j] + 1);
  }
  printf("\n");*/
}

static int dbscan() {
  int cid = 1;
  int atom_in_c[POINT_NUM];
  int visited[POINT_NUM] = {0};
  for (int i = 0; i < POINT_NUM; i++) {
    if (visited[i] == 1) continue;
    visited[i] = 1;
    int *neighbors = calloc(POINT_NUM, sizeof(int)), size = 0, unclustered_size = 0;
    get_neighbors(i, neighbors, &size, &unclustered_size);
    if (DEBUG) printf("[DEBUG] id: %d, size: %d, unclustered size: %d\n", i + 1, size, unclustered_size);
    if (unclustered_size < MIN_POINTS) {
      atoms[i].cluster_id = -1;
    } else {
      atoms[i].cluster_id = cid;
      atom_in_c[cid]++;
      // expand cluster
      for (int j = 0; j < size; j++) {
        int nid = neighbors[j];
        if (visited[nid] == 1) continue;
        visited[nid] = 1;
        int *n_neighbors = calloc(POINT_NUM, sizeof(int)), n_size = 0, n_unclustered_size = 0;
        get_neighbors(nid, n_neighbors, &n_size, &n_unclustered_size);
        if (DEBUG) printf("[DEBUG] id: %d, n size: %d, n unclustered size: %d\n", nid + 1, n_size, n_unclustered_size);
        if (n_size >= MIN_POINTS) {
          for (int k = 0; k < n_size; k++) {
            int n_nid = n_neighbors[k];
            // check if n_nid is already in neighbor array to avoid infinite loop
            int add_new_neighbor = 1;
            for (int l = 0; l < size; l++) {
              if (neighbors[l] == n_nid) {
                add_new_neighbor = 0;
                break;
              }
            }
            if (add_new_neighbor) neighbors[size++] = n_nid;
          }
        }
        free(n_neighbors);
        if (atoms[nid].cluster_id == 0) { // not yet a member of any cluster
          atoms[nid].cluster_id = cid;
          atom_in_c[cid]++;
        }
      }
      // expansion over, increase the cluster id
      cid++;
    }
    free(neighbors);
  }
  if (DEBUG) { // print the number of atoms in clusters
    for (int m = 1; m < cid; m++) {
      printf("%d, ", atom_in_c[m]);
    }
    printf("The number of clusters: %d\n", cid - 1);
  }
  return cid;
}

int main() {
  fetch_data(fopen("x.txt", "r"));
  calculate_distances();
  dbscan();
  export_atoms();
}
