
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

long count = 0;
long per_thread;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

double randomnum(double a, double b) {
	double r;
	r = ((b - a) * ((double)rand() / (double)RAND_MAX)) + a;
	return r;
}

void *runner() {
	long count_thread = 0;
	unsigned int rand_state = rand();
    	long i;
    	for (i = 0; i < per_thread; i++) {
		double x = randomnum(-1,1);
		double y = randomnum(-1,1);
        	if (x * x + y * y < 1) {
            		count_thread++;
        	}
    	}
    	pthread_mutex_lock(&mutex);
    	count += count_thread;
    	pthread_mutex_unlock(&mutex);
}


int main(int argc, const char *argv[])
{
    	long long int n = atol(argv[1]);
    	int nthread = atoi(argv[2]);
    	per_thread = n / nthread;

	time_t start = time(NULL);
	srand((unsigned)time(NULL));

    	pthread_t *threads = malloc(nthread * sizeof(pthread_t));

	pthread_attr_t attr;
        pthread_attr_init(&attr);

	int i;
    	for (i = 0; i < nthread; i++) {
        	pthread_create(&threads[i], &attr, runner, (void *) NULL);
    	}

    	for (i = 0; i < nthread; i++) {
        	pthread_join(threads[i], NULL);
    	}

    	pthread_mutex_destroy(&mutex);
    	free(threads);
	
	printf("TOTAL NUMBER OF THREADS: %d \n", nthread);
	printf("DARTS INSIDE THE CIRCLES: %ld \n",count);
	printf("VALUE OF PI IS: %f\n", (4. * (double)count) / ((double)per_thread * nthread));
	printf("TOTAL TIME IS: %d SEC\n", (unsigned int)(time(NULL) - start));

    return 0;
}
