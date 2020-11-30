#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

int main(int argc, char *argv[])
{
	int		epsilon, minpts, num_points, error, i, rand_seed;
	double		low, high, x, y, z;
	struct timeval	rand_init_timestamp;
	char		filename[256], output[256];
	FILE		*file;


	printf("Enter value for epsilon: ");
	scanf("%d", &epsilon);

	if (epsilon < 1) {
		printf("Invalid value for \"epsilon\".\n");
		exit(1);
	}

	printf("Enter value for minimum number of points that consist a cluster: ");
	scanf("%d", &minpts);

	if (minpts < 1) {
		printf("Invalid value for \"minpts\".\n");
		exit(1);
	}

	printf("Enter number of points to be clustered: ");
	scanf("%d", &num_points);

	if (num_points < 1) {
		printf("Invalid value for \"num_points\".\n");
		exit(1);
	}

	printf("Enter lowest value for point coordinates : ");
	scanf("%lf", &low);

	printf("Enter highest value for point coordinates: ");
	scanf("%lf", &high);

	printf("Enter output file name: ");
	scanf("%s", filename);

	file = fopen(filename, "w");
	if (file == NULL) {
		printf("Could not open file \"%s\".\n", filename);
		exit(1);
	}

	error = gettimeofday(&rand_init_timestamp, NULL);
        if (error != 0) {
                printf("gettimeofday() returned an error.\n");
                exit(1);
        }

        rand_seed = rand_init_timestamp.tv_sec * 1000000 + rand_init_timestamp.tv_usec;

	srand48(rand_seed);

	sprintf(output, "%d %d %d\n", epsilon, minpts, num_points);
	fwrite(output, sizeof(char), strlen(output), file);
	for (i = 0; i < num_points - 1; i++) {
		x = low + (high - low) * drand48();
		y = low + (high - low) * drand48();
		z = low + (high - low) * drand48();
		sprintf(output, "%.2f %.2f %.2f\n", x, y, z);
		fwrite(output, sizeof(char), strlen(output), file);
	}

	x = low + (high - low) * drand48();
	y = low + (high - low) * drand48();
	z = low + (high - low) * drand48();
	sprintf(output, "%.2f %.2f %.2f", x, y, z);
	fwrite(output, sizeof(char), strlen(output), file);

	fclose(file);

	return(0);
}
