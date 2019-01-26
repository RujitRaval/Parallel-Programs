#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

double random(double a, double b) {
	double r;
	r = ((b - a) * ((double)rand() / (double)RAND_MAX)) + a;
	return r;
}

int main(int argc, char*argv[]) {
	long long int n = 10000000000;
	long long int npernode;
	int rank, size, count = 0, i;
	double pi_current, pi_sum, x, y, start, end;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0) {
		printf("  Number of processes: %d \n", size);
		npernode = n * (rank + 256);
		MPI_Bcast(&npernode, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		start = MPI_Wtime();
	}

	for (i = 1; i <= n; i++) {
		x = random(-1, 1);
		y = random(-1, 1);
		if (((x*x) + (y*y)) <= 1.0) {
			count++;
		}
	}

	pi_current = 4.0 * (double)count / (double)n;
	MPI_Reduce(&pi_current, &pi_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		pi_sum = pi_sum / size;
		end = MPI_Wtime();
		printf("Estimated Value of PI  : %f\n", pi_sum);
		printf("Time    : %f\n\n", end - start);
	}
	MPI_Finalize();
	return 0;
}