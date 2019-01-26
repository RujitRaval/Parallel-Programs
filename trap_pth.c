#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <time.h>

const long max_thread = 512;
double a, b, h, Mutex_Trap = 0, Sem_Trap = 0, Busy_Trap = 0;
pthread_mutex_t mutex;
sem_t semaphore;
int flag = 0;
unsigned long long int n;  

double f(double x)
{
  double return_val = x*x;
  return return_val;
}     

long long my_first(long rank,long count,long long int n)
{
	return ((rank)*(n)/(count));
}

long long my_last(long rank, long count, long long int n)
{
	return (my_first((rank)+1,count,n)-1);
}

void* mutex_trapezoid(void* rank)
{ 
  long thread_rank = (long)rank;
  double sum = 0.0;
  long long i;
  int first = (int)thread_rank;
  long long my_first_i = my_first(thread_rank, max_thread, n);
  long long my_last_i = my_last(thread_rank, max_thread, n);
  if( first == 1)
    sum += (f(a)+f(b))/2.0;
  for( i= my_first_i; i <= my_last_i; i++)
    sum += f(a+(i*h));
  sum = sum * h;
  pthread_mutex_lock(&mutex);
  Mutex_Trap += sum;
  pthread_mutex_unlock(&mutex);
  return NULL;
}
void* sem_trapezoid(void *rank)
{
  long thread_rank = (long)rank;
  double sum = 0.0;
  long long i;
  int first = (int)thread_rank;
  long long my_first_i = my_first(thread_rank, max_thread, n);
  long long my_last_i = my_last(thread_rank, max_thread, n);

  if( first == 1)
    sum += (f(a)+f(b))/2.0;
  for( i= my_first_i; i <= my_last_i; i++)
    sum += f(a+(i*h));
  sum = sum * h;
  sem_wait(&semaphore);
  Sem_Trap += sum;
  sem_post(&semaphore);

  return NULL;
}

void* busy_wait_trapezoid(void* rank)
{
  long thread_rank = (long)rank;
  double sum = 0.0;	
  long long i;
  int first = (int)thread_rank;

  long long my_first_i = my_first(thread_rank, max_thread, n);
  long long my_last_i = my_last(thread_rank, max_thread, n);

  if( first == 1)
    sum += (f(a)+f(b))/2.0;

  for( i= my_first_i; i <= my_last_i; i++)
    sum += f(a+(i*h));
  sum = sum * h;

  while(flag != thread_rank );
  Busy_Trap += sum;
  flag = (flag+1) % max_thread;
  return NULL;
}

void mutex_trap()
{
  time_t start = time(NULL);
  pthread_t* thread_handles;
  pthread_mutex_init(&mutex, NULL);
  thread_handles = malloc( max_thread * sizeof(pthread_t));
  long i;
  for( i=0; i < max_thread; i++)
    pthread_create( &thread_handles[i], NULL, mutex_trapezoid, (void*) i );
  for( i=0; i < max_thread; i++)
    pthread_join( thread_handles[i], NULL);
  free(thread_handles);
  pthread_mutex_destroy(&mutex);
  time_t end = time(NULL);
  printf("TRAPEZOID WITH MUTEX\n");
  printf("Number of threads = %ld \n", max_thread);
  printf("Number of trapezoids = %llu \n", n);
  printf("Integral from %f to %f = %.15f\n",a, b, Mutex_Trap); 
  printf("Time is : %d Seconds \n\n", (unsigned int)( end- start));  
}


void sem_trap()
{
	time_t start = time(NULL);
    pthread_t* thread_handles;
  sem_init(&semaphore, 0, 1);
  thread_handles = malloc( max_thread * sizeof(pthread_t));
  long i;
  for( i=0; i < max_thread; i++)
    pthread_create( &thread_handles[i], NULL, sem_trapezoid, (void*) i);
  for( i=0; i < max_thread; i++)
    pthread_join( thread_handles[i], NULL);
  free(thread_handles);
  sem_destroy(&semaphore);
   time_t end = time(NULL);
  printf("TRAPEZOID WITH SEMAPHORE\n");
  printf("Number of threads = %ld \n", max_thread);
  printf("Number of trapezoids = %llu \n", n);
  printf("Integral from %f to %f = %.15f\n",a, b, Sem_Trap);  
  printf("Time is : %d Seconds \n\n", (unsigned int)( end- start));  
}

void busy_wait_trap()
{
	time_t start = time(NULL);
    pthread_t* thread_handles;
  thread_handles = malloc( max_thread * sizeof(pthread_t));
  long i;
  for( i=0; i < max_thread; i++)
    pthread_create( &thread_handles[i], NULL, busy_wait_trapezoid, (void*) i);
  for( i=0; i < max_thread; i++)
    pthread_join( thread_handles[i], NULL);
  free(thread_handles);
   time_t end = time(NULL);
  printf("TRAPEZOID new WITH BUSY_WAIT\n");
  printf("Number of threads = %ld \n", max_thread);
  printf("Number of trapezoids = %llu \n", n);
  printf("Integral from %f to %f = %.15f\n", a, b, Busy_Trap); 
printf("Time is : %d Seconds \n\n", (int)( end- start));    
}



int main( int argc, char **argv )
{
  printf("Enter a, b, and n\n");
  scanf("%lf",&a);
  scanf("%lf",&b);
  scanf("%llu",&n);
  h = (b-a)/n;
  mutex_trap();
  sem_trap();
  busy_wait_trap();
  return 0;
}
