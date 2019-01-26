// RUJIT RAVAL
// Student ID # 1512338
// HIGH PERFORMANCE COMPUTING
// ASSIGNMENT - 0
// PROGRAM - 4
// Estimating speed of *, /, sqrt & sin

#include<stdio.h>
#include<time.h>
#include<math.h>

//Multiplication function
float multiplication(long int x){
	return x*5;
}

//Division function
float division(long int x){
	return x/2;
	
}

//Function for sqrt()
float sqrroot(long int x){
	return sqrt(x); 
}

//Function for sin()
float sin_value(long int x){
	return sin(x);
}


int main(int argc, char **argv){
	clock_t t_mul, t_div, t_sqrt, t_sin;
	long int data[1000000];
	long int i = 1;
	
	//Calling multiplication
	t_mul =clock();
	for (i=1; i<=1000000; i++){
		multiplication(data[i]);
	}
	t_mul = clock()-t_mul;

	//Calling division
	t_div=clock();
	for (i=1; i<=1000000; i++){
		division(data[i]);
			
	}
	t_div = clock()-t_div;
	
	//Calling function for sqrt
	t_sqrt=clock();
	for (i=1; i<=1000000; i++){
		sqrroot(data[i]);
			
	}
	t_sqrt = clock()-t_sqrt;

	//Calling function for sin
	t_sin=clock();
	for (i=1; i<=1000000; i++){
		sin_value(data[i]);
			
	}
	t_sin = clock()-t_sin;
	
	//Calculating time
	double timetaken_mul = ((double)t_mul)/(CLOCKS_PER_SEC/1000);
	double timetaken_div = ((double)t_div)/(CLOCKS_PER_SEC/1000);
	double timetaken_sqrt = ((double)t_sqrt)/(CLOCKS_PER_SEC/1000);
	double timetaken_sin = ((double)t_sin)/(CLOCKS_PER_SEC/1000);

	//Printing time
	printf("Time taken for multiplication = %f milliseconds. \n", timetaken_mul);
	printf("Time taken for division = %f milliseconds. \n", timetaken_div);
	printf("Time taken for sqrt function = %f milliseconds. \n", timetaken_sqrt);
	printf("Time taken for sin function = %f milliseconds. \n", timetaken_sin);
	return 0;
}

// For compile: gcc -Wall -g prog4.c -o prog4 -lm

