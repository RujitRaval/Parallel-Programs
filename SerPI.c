#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void Usage(char* prog_name);



int main(int argc, char* argv[]){
 	 double size = 10000000000;
 	 double x,y;
 	 int toss;
 	 double count = 0;
 	 double pi;
 	 double runtime;
	 double distance;

 	if (argc != 3) Usage(argv[0]);
 	thread_count = strtol(argv[1], NULL, 10);
   	n = strtoll(argv[2], NULL, 10);
   	if (thread_count < 1 || n < 1) Usage(argv[0]);
  
	
   	srandom(0);
#  pragma omp parallel for num_threads(thread_count) \
      reduction(+: count) private(x, y, distance)

	clock_t begin = clock();
  for(toss = 0; toss < size; toss++){
    x = (double)drand48();
    y = (double)drand48();
    distance = (x*x)+(y*y)
    if(distance<1){
      ++count;
    }
  }
  clock_t end = clock();
  runtime = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Darts inside the circle are: %f \n", count);
  pi = (double)4 * (count / size);
  printf("Total time spent = %f \n", runtime);
  printf("pi = %f", pi);
return 0;
}

void Usage(char* prog_name) {
   fprintf(stderr, "usage: %s <thread_count> <number_tosses>\n", prog_name);  /* Change */
   exit(0);
}  /* Usage */
