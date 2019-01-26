// RUJIT RAVAL
// Student ID # 1512338
// HIGH PERFORMANCE COMPUTING
// ASSIGNMENT - 0
// PROGRAM - 2
// Archimedes algorithm to estimate pi

#include<stdio.h>
#include<math.h>
#include<time.h>

int main(int argc, char **argv){
	clock_t t_cal;
	printf("Archimedes Algorithm: \n");
	
	//Assigning values 1 to 100
	int n[100];
	for(int i=1; i<=100; i++){
		n[i]= i;
	}		
	
	//Inscribe Polygon
	float ins[100];
	
	//Circumscribe Polygon
	float cir[100];
	
	//Calculating value of PI
	t_cal =clock();	
	for(int i=1; i<=100; i++){
		double rad = (double)((n[i]-2)* 180)/(2*n[i]);
		double deg = (double)(rad * M_PI) / 180.0;
		ins[i] = n[i]* sin(deg)* cos(deg) ; //Lowerbound Value
		cir[i] = n[i] / tan (deg); //Upperbound Value
	}
	t_cal = clock()-t_cal;

	//Printing Values
	printf("SIDES \t INSCRIBED \t CIRCUMSCRIBED \n ------------------------------------------ \n");
	for(int i=1; i<=100; i++){
		printf("%d \t %f \t %f \n", n[i], ins[i], cir[i]);	
	}
	
	//Calculating time
	double timetaken_cal = ((double)t_cal)/(CLOCKS_PER_SEC/1000);
	//Printing time
	printf("Time taken for calculation is %f milliseconds. \n", timetaken_cal);
	printf("\n");
	return 0;
}
