// RUJIT RAVAL
// Student ID # 1512338
// HIGH PERFORMANCE COMPUTING
// ASSIGNMENT - 0
// PROGRAM - 3
// Matrix - vactor multiplication

#include<stdio.h>
#include<stdlib.h>

int main(int argc, char **argv){
	FILE *inputfile; //Pointer to file
	inputfile = fopen("mv.txt","r"); //Open a file
	
	int rows,cols;
	fscanf(inputfile, "%d %d", &rows, &cols); //Getting value of #of rows and columns
	
	int matrix[rows][cols];
	
	//initializing matrix
	for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++){
			fscanf(inputfile, "%d", &matrix[i][j]);
		}
	} 
	
	int v_size;
	fscanf(inputfile, "%d", &v_size); //Getting size of vector
	int v[v_size];	
	
	//initializing vector
	for(int i=0; i<v_size; i++){
		fscanf(inputfile, "%d", &v[i]);
	}
	
	//Function for multiplication
	int ans[rows][cols];
	for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++){
			ans[i][j] = matrix[i][j] * v[i];
		}
	} 
	
	//Printing the value
	for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++){
			printf("%d \t", ans[i][j]);
		}
		printf("\n");
	} 
	printf("\n");
	return 0;
}
