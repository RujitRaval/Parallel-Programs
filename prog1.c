// RUJIT RAVAL
// Student ID # 1512338
// HIGH PERFORMANCE COMPUTING
// ASSIGNMENT - 0
// PROGRAM - 1
// Print "Hello, <NAME>" where NAME is input from keyboard.


#include<stdio.h>

int main(int argc, char **argv){
	if (argc == 1){
		printf("No argument applied ! \n");
		return 1;
	}
	
	//Print 1st argument
	printf("Hello, %s",argv[1]);
	
	//Print multiple arguments (if given)	
	for(int i=2; i<argc; i++){
		printf(" %s", argv[i]);
	}	

	//New line after the end of all arguments
	printf("\n");
	return 0;
}
