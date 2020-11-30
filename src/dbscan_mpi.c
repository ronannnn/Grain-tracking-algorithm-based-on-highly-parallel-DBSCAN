/**
 * Author:   Ruonan Chen (ruonanch@usc.edu)
 * Date:     11/29/20
 * Filename: dbscan_mpi.c
 */

#include "mpi.h"
#include "dbscan_serial.c"

int main(int argc, char *argv[]) {
  int rank, size;
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  //Read file of points coordinates into directory
  FILE *in = fopen("input.txt", "r");

  point_t *points;
  double eps = EPS;
  int minPts = MIN_POINTS;
  Atom atoms[POINT_NUM] = fetch_data(in)
  fclose(in);

  // dividing points set into equal size subsets
  MPI_Scatter(points, num, mpi_point_t_type, myLocal, num, mpi_point_t_type, 0, MPI_COMM_WORLD);
}