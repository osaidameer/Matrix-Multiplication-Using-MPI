#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<limits.h>
#include<mpi.h>

#define MAX_LEN 256

struct matrixStruct {
    char id;
    int matrix[50][50];
};

struct listStruct {
    struct matrixStruct listObj[7];
};

struct listStruct matrices;

int matrixRes[50][50] = {0};
int resRows = 0;
int resCols = 0;
int line_counter= 0;
int bracket[50][50] = {0};
int matrixA[50][50], matrixB[50][50], matrixC[50][50];
int rank, size;

void file_reading(int rows[], int columns[], FILE *fp, char lines[][10]){
    char buffer[MAX_LEN];
    int i;
    while (fgets(buffer, MAX_LEN, fp)) {
        //printf("%s\n", buffer);
        for(i=0; i<strlen(buffer); i++){
            lines[line_counter][i] = buffer[i];
        }
        rows[line_counter] = ((buffer[0] - '0') * 10) + ((buffer[1] - '0') * 1);
        columns[line_counter] = ((buffer[5] - '0') * 10) + ((buffer[6] - '0') * 1);
        if (rows[line_counter] < 20 || rows[line_counter] > 50 || columns[line_counter] < 20 || columns[line_counter] > 50){
            printf("\nMatrix at line %d has invalid dimensions\n", line_counter);
            exit(1);
        }
        line_counter++;
    }
}

void check_dimensions(const int dims[], int size, const int rows[], const int columns[]){
    // check uniqueness
    // check if matrices are square
    // check if there are valid number of matrices
    if (line_counter < 5 || line_counter > 7){
        printf("\nERROR! Number of Matrices must be between 5 and 7, number is %d\n", line_counter);
        exit(1);
    }
    int i;
    int j;
    for(i=0; i<size; i++){
        if (i < size - 1){
            if (rows[i] == columns[i]){
                printf("\nERROR! Matrix %d is a square Matrix\n", i+1);
                exit(1);
            }
        }
        if (i < size -2){
            if (columns[i] != rows[i+1]){
                printf("\nERROR! Matrices %d and %d are not suitable for multiplication\n", i+1, i+2);
                exit(1);
            }
        }
        for(j=i + 1; j<size; j++){
            if (dims[j] == dims[i]) {
                printf("\nERROR! Matrices %d and %d have the same dimensions\n", i+1, j+1);
                exit(1);
            }
        }
    }
}

void set_dimensions(int dimensions[], const int matrix_rows[], const int matrix_columns[]){
    //printf("Dimensions: \n");
    int i;
    for (i=0; i<line_counter; i++){
        dimensions[i] = matrix_rows[i];
        //printf("%d ", dimensions[i]);
    }
    dimensions[line_counter] = matrix_columns[line_counter - 1];
    //printf("%d\n", dimensions[line_counter]);
}

void initialise_matrix(int index, int rows, int col, char id){
    if(rank==0){
	    matrices.listObj[index].id = id;
	    printf("\nMatrix %c\n", id);
	    int i;
	    int j;
	    for (i=0; i<rows; i++){
		for (j=0; j<col; j++){
		    matrices.listObj[index].matrix[i][j] = rand() % 5 + 1;
		    printf("%d ", matrices.listObj[index].matrix[i][j]);
		}
		printf("\n");
	    }
	    printf("\n\n");
    }
}

void populate_matrices(int matrix_rows[], int matrix_columns[], int size){
    char id = 'A';
    int i;
    for(i=0; i<size; i++, id++){
        initialise_matrix(i, matrix_rows[i], matrix_columns[i], id);
    }
}


int MatrixChainOrder(int p[], int i, int j)
{
    if (i == j)
        return 0;
    int k;
    int min = INT_MAX;
    int count;

    for (k = i; k < j; k++)
    {
        count = MatrixChainOrder(p, i, k)
                + MatrixChainOrder(p, k + 1, j)
                + p[i - 1] * p[k] * p[j];
//		printf("%d,", count);
        if (count < min)
//        	printf("\nHERE %d \n", count);
            min = count;
//            printf("%d\n", k);
//            bracket[i][j] = k;
    }

    return min;
}

int main(int argc, char **argv){
    //srand(time(0));
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Request request = MPI_REQUEST_NULL;    
    
    int slave_count = size - 1;
    double t1, t2;
    
    FILE* fp;
    char lines[7][10] = {0};
    int matrix_rows[7] = {0};
    int matrix_columns[7] = {0};

    fp = fopen("matrices.txt", "r");
    if (fp == NULL) {
        perror("Failed: ");
        return 1;
    }
    
   if (rank==0){ 
   	printf("\n---------------------------------------------------------------------BLOCKING---------------------------------------------------------------------\n");    
   }    
   sleep(1);

    file_reading(matrix_rows, matrix_columns, fp, lines);

    int *dimensions;
    dimensions = (int *)malloc(sizeof(int) * line_counter + 1);
    set_dimensions(dimensions, matrix_rows, matrix_columns);
    check_dimensions(dimensions, line_counter + 1, matrix_rows, matrix_columns);
    populate_matrices(matrix_rows, matrix_columns, line_counter);
    
    if (rank == 0){
	    printf("\nMatrices Populated!\n");
	    int cost = MatrixChainOrder(dimensions, 1, line_counter);
	    printf("\nMinimum Cost of Multiplying the Matrices: %d\n", cost);
	    printf("\nDimensions: [ ");
	    for(int i=0; i<line_counter + 1; i++){
	    	printf("%d ", dimensions[i]);
	    }
	    printf("]\n");
    }

    
    int rows = 0;
    int offset = 0;
    
    
    if (rank == 0){
    	t1 = MPI_Wtime();
    	//int resR, resC;
    	rows = matrix_rows[0]/slave_count;
	for(int i=0; i<matrix_rows[0]; i++){
		for(int j=0; j<matrix_columns[0]; j++){
			matrixC[i][j] = matrices.listObj[0].matrix[i][j];
		}
	}
    	for(int i=0; i<line_counter - 1; i++){
    		offset = 0;

    		for(int j=1; j<=slave_count; j++){
			MPI_Send(&offset, 1, MPI_INT, j, 1, MPI_COMM_WORLD);
			MPI_Send(&rows, 1, MPI_INT, j, 1, MPI_COMM_WORLD);
			MPI_Send(&matrix_columns[i], 1, MPI_INT, j, 1, MPI_COMM_WORLD);	
			MPI_Send(&matrix_columns[i+1], 1, MPI_INT, j, 1, MPI_COMM_WORLD);						
			int *mA = (int *)malloc(rows * matrix_columns[i] * sizeof(int));
			
			for(int k=offset, m=0; k<offset+rows; k++){
				for(int l=0; l<matrix_columns[i]; l++, m++){
					mA[m] = matrixC[k][l];
				}
			}
			//send previous received result, which would have columns of the previous matrix as well				
			MPI_Send(mA, rows * matrix_columns[i], MPI_INT, j, 1, MPI_COMM_WORLD);	

			//send B
			int *mB = (int *)malloc(matrix_columns[i+1] * matrix_rows[i+1] * sizeof(int));
			
			for(int k=0, m=0; k<matrix_rows[i+1]; k++){
				for(int l=0; l<matrix_columns[i+1]; l++, m++){
					mB[m] = matrices.listObj[i+1].matrix[k][l];
				}
			}			
			
			MPI_Send(mB, matrix_rows[i+1] * matrix_columns[i+1], MPI_INT, j, 1, MPI_COMM_WORLD);	

			printf("\n\nSending from row [%d] till [%d] to process %d", offset, offset + rows, j);
			offset += rows;
    		}
    		//clearing C before receiving from all processes 
		//memset(matrixC, 0, sizeof(matrixC));
		
    		for(int j=1; j<=slave_count; j++){
    			int *mC = (int *)malloc(rows * matrix_columns[i+1] * sizeof(int));
    			MPI_Recv(&offset, 1, MPI_INT, j, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    			MPI_Recv(&rows, 1, MPI_INT, j, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    			MPI_Recv(mC, rows * matrix_columns[i+1], MPI_INT, j, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    			
    			for(int k=offset, m=0; k<offset+rows; k++){
    				for(int l=0; l<matrix_columns[i+1]; l++, m++){
    					matrixC[k][l] = mC[m];
    				}
    			}
    		}	    		
    		
    	}
    	printf("\n\nMatrix C: \n");
    	for(int i=0; i<matrix_rows[0]; i++){
    		for(int j=0; j<matrix_columns[line_counter - 1]; j++){
    			printf("%d ", matrixC[i][j]);
    		}
    		printf("\n");
    	}
    }
   
   else if (rank > 0){
   	//int colA, colB, colC;
   	int c1, c2;
   	

   	for (int i=0; i<line_counter - 1; i++){
   		int resMatrix[50][50] = {0};	
   		MPI_Recv(&offset, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   		MPI_Recv(&rows, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   		MPI_Recv(&c1, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
     		MPI_Recv(&c2, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 	
     		
     		int *mA = (int *)malloc(rows * c1 * sizeof(int));
     		int *mB = (int *)malloc(c1 * c2 * sizeof(int));
     		int *mR = (int *)malloc(rows * c2 * sizeof(int));
     		
   		MPI_Recv(mA, rows * c1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   		MPI_Recv(mB, c1 * c2, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
   		
   		for(int j=0, l=0; j<rows; j++){
   			for(int k=0; k<c1; k++, l++){
   				matrixA[j][k] = mA[l];
   			}
   		}
   		
   		for(int j=0, l=0; j<c1; j++){
   			for(int k=0; k<c2; k++, l++){
   				matrixB[j][k] = mB[l];
   			}
   		}
   		
   		/*
   		if(rank==1){
	   		printf("\n\nSlave 1 Received: \n");
	   		for(int i=0; i<rows; i++){
	   			for(int j=0; j<matrix_columns[0]; j++){
	   				printf("% d", matrixA[i][j]);
	   			}
	   			printf("\n");
	   		}
   		}*/
   		
   		//MULTIPLICATION LOOP 
   		for(int j=0; j<rows; j++){
   			for(int k=0; k<c2; k++){
   				for(int l=0; l<c1; l++){
   					resMatrix[j][k] += matrixA[j][l] * matrixB[l][k];
   				}
   			}
   		}
   		
   		for(int j=0, l=0; j<rows; j++){
   			for(int k=0; k<c2; k++, l++){
   				mR[l] = resMatrix[j][k];
   			}
   		}
   		
   		//memset(matrixA, 0, sizeof(matrixA));
   		//memset(matrixB, 0, sizeof(matrixB));
   		
   		MPI_Send(&offset, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
     		MPI_Send(&rows, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
     		MPI_Send(mR, rows * c2, MPI_INT, 0, 2, MPI_COMM_WORLD); 		
   	}
   }
   
    if (rank==0){
    	t2 = MPI_Wtime();
   	FILE *fp1;
   	   
   	fp1 = fopen("blocking.txt", "a");
	printf("\n\nTotal Execution time for blocking Multiplicaiton %f\n", t2 - t1);
	fprintf(fp1, "Number of Processes: %d, Execution Time %f\n", slave_count, t2-t1);
	fclose(fp1);    	
    }   

    //printf("\n%d\n", matrix_columns[line_counter - 1]);
    MPI_Finalize();
    fclose(fp);
    return 0;
}
