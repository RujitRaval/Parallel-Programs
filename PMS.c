#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

void merge(int *, int *, int, int, int);
void mergeSort(int *, int *, int, int);

int main(int argc, char** argv) {
	long long int n = 100000000;
	double start, end;
	long long int *original_array = (long long int*)malloc(n * sizeof(long long int));
	long int i;

	int comm_rank;
	int comm_size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	start = MPI_Wtime();
	
	for (i = 0; i < n; i++) {
		original_array[i] = rand() % n;
	}
	long int size = n / comm_size;
	int *sub_array = (int*)malloc(size * sizeof(int));
	MPI_Scatter(original_array, size, MPI_INT, sub_array, size, MPI_INT, 0, MPI_COMM_WORLD);

	int *tmp_array = (int*)malloc(size * sizeof(int));
	mergeSort(sub_array, tmp_array, 0, (size - 1));

	int *sorted = NULL;
	if (comm_rank == 0) {
		sorted = (int*)malloc(n * sizeof(int));
	}

	MPI_Gather(sub_array, size, MPI_INT, sorted, size, MPI_INT, 0, MPI_COMM_WORLD);

	if (comm_rank == 0) {
		int *other_array = (int*)malloc(n * sizeof(int));
		mergeSort(sorted, other_array, 0, (n - 1));
		end = MPI_Wtime();
		printf("Time is %f \n\n", end-start);
		free(sorted);
		free(other_array);
	}
	free(original_array);
	free(sub_array);
	free(tmp_array);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	system("pause");
	return 0;
}

void merge(int *a, int *b, int l, int m, int r) {
	int h, i, j, k;
	h = l;
	i = l;
	j = m + 1;
	while ((h <= m) && (j <= r)) {
		if (a[h] <= a[j]) {
			b[i] = a[h];
			h++;
		}
		else {
			b[i] = a[j];
			j++;
		}
		i++;
	}

	if (m < h) {
		for (k = j; k <= r; k++) {
			b[i] = a[k];
			i++;
		}
	}
	else {
		for (k = h; k <= m; k++) {
			b[i] = a[k];
			i++;
		}
	}

	for (k = l; k <= r; k++) {
		a[k] = b[k];
	}
}

void mergeSort(int *a, int *b, int l, int r) {
	int m;
	if (l < r) {
		m = (l + r) / 2;
		mergeSort(a, b, l, m);
		mergeSort(a, b, (m + 1), r);
		merge(a, b, l, m, r);
	}
}