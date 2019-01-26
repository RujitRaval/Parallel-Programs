#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

void Usage(char* prog_name);

int thread_count;

int main(int argc, char* argv[]){
 	 
 	 double x,y;
	 int i;
 	 int toss;
 	 double count = 0;
 	 double pi;
 	 double runtime;
	 double distance;

 	if (argc != 3) Usage(argv[0]);
 	thread_count = strtol(argv[1], NULL, 10);
   	toss = strtoll(argv[2], NULL, 10);
   	if (thread_count < 1 || toss < 1) Usage(argv[0]);
  
	
   	srandom(0);

	clock_t begin = clock();
#  pragma omp parallel for num_threads(thread_count) \
      reduction(+: count) private(x, y, distance)
  for(i = 0; i < toss; i++){
    x = (double)drand48();
    y = (double)drand48();
    distance = (x*x)+(y*y);
    if(distance<1){
      ++count;
    }
  }
  clock_t end = clock();
  runtime = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Darts inside the circle are: %f \n", count);
  pi = (double)4 * (count / toss);
  printf("Total time spent = %f \n", runtime);
  printf("pi = %f \n", pi);
return 0;
}

void Usage(char* prog_name) {
   fprintf(stderr, "usage: %s <thread_count> <number_tosses>\n", prog_name);  /* Change */
   exit(0);
}  /* Usage */
