#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void iterate(float** hotOld, float** hotNew, int start, int end);

bool isStable(float** hot, int start, int end);

void initialize(float** rowsOld,float** rowsNew, int start, int end);

float **initArray(int rows, int cols);

double When();

//#pragma omp for private(row,col)

//using namespace std;

const int NUM_OF_ROWS = 16384;
const int NUM_OF_COLS = 16384;
const float EPSILON = 0.1;
int numCells;

struct arg_struct {
   float** hotOld;
   float** hotNew;
   int start;
   int end;
};

struct ret_struct {
   bool isStable;
   int numCells;
};

int main(int argc, char *argv[]) {

   int iproc, nproc;
   char host[255], message[55];
   MPI_Status status;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nproc);
   MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

   gethostname(host,253);
   printf("I am proc %d of %d running on %s\n", iproc+1, nproc,host);


   int i, j, proc;
   numCells = 0;
   int iterations = 0;
   
   //Master Process
   if(iproc == 0) {

      double begin = When();  

      int range = NUM_OF_ROWS/(nproc-1);
      
      //initialize
      for(proc=1;proc<nproc;proc++) {
    	  int start = range * (proc - 1) - 1;
		  if(start < 0) {
			  start = 0;
		  }
		  int end = range * proc + 1;
		  if(end > NUM_OF_COLS) {
			  end = NUM_OF_COLS;
		  }
		  printf("%d sending to %d\n",iproc,proc);
    	  MPI_Send(&start,1,MPI_INT,proc,0,MPI_COMM_WORLD);
    	  MPI_Send(&end,1,MPI_INT,proc,0,MPI_COMM_WORLD);
      }
      
      bool stable = true;
      do{
    	  float** rows[(nproc-2) * 2];
    	  int counter = 0;
    	  for(proc=1;proc<nproc;proc++) {
    		  if(proc!=1 && proc != nproc-1) {
    			  float** temp = initArray(1,NUM_OF_COLS);
    			  MPI_Recv(&(temp[0][0]),NUM_OF_COLS,MPI_INT,proc,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    			  rows[counter] = temp;
    			  counter++;
    			  free(temp);
    		  }
			  float** temp = initArray(1,NUM_OF_COLS);
			  MPI_Recv(&(temp[0][0]),NUM_OF_COLS,MPI_INT,proc,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			  rows[counter] = temp;
			  counter++;
			  free(temp);
    	  }
    	  counter = 0;
    	  for(proc=1;proc<nproc;proc++) {
    		  if(proc == 1) {
    			  float** temp = rows[counter+1];
				  MPI_Send(&(temp[0][0]),NUM_OF_COLS,MPI_INT,proc,0,MPI_COMM_WORLD);
				  counter++;
				  free(temp);
    		  }
    		  if(proc == nproc-1) {
    			  float** temp = rows[counter-1];
				  MPI_Send(&(temp[0][0]),NUM_OF_COLS,MPI_INT,proc,0,MPI_COMM_WORLD);
				  counter++;
				  free(temp);
    		  }
    		  else {
    			  float** temp = rows[counter-1];
    			  MPI_Send(&(temp[0][0]),NUM_OF_COLS,MPI_INT,proc,0,MPI_COMM_WORLD);
    			  counter++;
    			  temp = rows[counter+1];
				  MPI_Send(&(temp[0][0]),NUM_OF_COLS,MPI_INT,proc,0,MPI_COMM_WORLD);
				  counter++;
				  free(temp);
    		  }
    	  }
    	  for(proc=1;proc<nproc;proc++) {
    		  bool tempBool = false;
    		  int tempCells = 0;
    		  MPI_Recv(&tempBool,1,MPI_INT,proc,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    		  MPI_Recv(&tempCells,1,MPI_INT,proc,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    		  numCells += tempCells;
    		  if(!tempBool) {
    			  stable = tempBool;
    			  numCells = 0;
    		  }
    	  }
    	  iterations++;
      }while(!stable);
      MPI_Finalize();
      double end = When(); 
   }
   
   //Slave Process
   else {
	   
	  //initialize
	  int start, end, i, j;
      MPI_Recv(&start, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      printf("%d: start=%d\n", iproc,start);
      MPI_Recv(&end, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      printf("%d: end=%d\n", iproc,end);
      float** rowsOld = initArray(end-start,NUM_OF_COLS);
      float** rowsNew = initArray(end-start,NUM_OF_COLS);
      float** temp = initArray(end-start,NUM_OF_COLS);
      initialize(rowsNew,rowsOld,start,end);
      printf("%d: initialized\n", iproc);
      
      bool stable = true;
      do{
    	  iterate(rowsOld,rowsNew,start,end);
    	  if(iproc != 1) {
    		  float** topRow = initArray(1,NUM_OF_COLS);
    		  for(j=0;j<NUM_OF_COLS;j++) {
    			  topRow[0][j] = rowsNew[1][j];
    		  }
    		  MPI_Send(&(topRow[0][0]),NUM_OF_COLS,MPI_INT,0,0,MPI_COMM_WORLD);
    		  MPI_Recv(&(topRow[0][0]),NUM_OF_COLS,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    		  for(j=0;j<NUM_OF_COLS;j++) {
    			  rowsNew[0][j] = topRow[0][j];
    		  }
    		  free(topRow);
    	  }
    	  if(iproc != nproc-1) {
    		  float** bottomRow = initArray(1,NUM_OF_COLS);
    		  for(j=0;j<NUM_OF_COLS;j++) {
				  bottomRow[0][j] = rowsNew[1][j];
			  }
			  MPI_Send(&(bottomRow[0][0]),NUM_OF_COLS,MPI_INT,0,0,MPI_COMM_WORLD);
			  MPI_Recv(&(bottomRow[0][0]),NUM_OF_COLS,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			  for(j=0;j<NUM_OF_COLS;j++) {
				  rowsNew[end-start][j] = bottomRow[end-start][j];
			  }
			  free(bottomRow);
    	  }
    	  bool tempBool = isStable(rowsNew,start,end);
    	  MPI_Send(&tempBool,1,MPI_INT,0,0,MPI_COMM_WORLD);
    	  MPI_Send(&numCells,1,MPI_INT,0,0,MPI_COMM_WORLD);
    	  MPI_Recv(&tempBool,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    	  if(!tempBool) {
    		  stable = tempBool;
    		  temp = rowsOld;
    		  rowsOld = rowsNew;
    		  rowsNew = temp;
    	  }
      }while(!stable);

      free(rowsOld);
      free(temp);
      free(rowsNew);
   }

//   double stop = When();
//   printf("%d number of cells\n", numCells);
//   printf("%d iterations\n", numOfIterations);
//   //printf("%f seconds\n", difftime(end, start));
//   printf("%f seconds\n", stop-begin);
   return 0;
}

float **initArray(int rows, int cols) {
	int i;
	float** array = (float**)malloc(rows*sizeof(float*));
	      for (i = 0; i < rows; i++) {
	         array[i] = (float*)malloc(cols * sizeof(float));
	      }
	 return array;
}

void initialize(float** rowsOld,float** rowsNew, int start, int end) {
   int i, j;
   //#pragma omp parallel for private(i,j)
   for (i = 0; i < end-start; i++) {
      for (j = 0; j < NUM_OF_COLS; j++) {
		if (j == NUM_OF_COLS - 1) {
		   rowsOld[i][j] = 0;
		   rowsNew[i][j] = 0;
		}
		else if (i+start == NUM_OF_ROWS - 1) {
		   rowsOld[i][j] = 100;
		   rowsNew[i][j] = 100;
		}
		else if (((i+start)%20) == 0) {
			rowsOld[i][j] = 100;
			rowsNew[i][j] = 100;
		}
		else if((j%20) == 0) {
		   rowsOld[i][j] = 0;
		   rowsNew[i][j] = 0;
		}
		else {
		   rowsOld[i][j] = 50;
		   rowsNew[i][j] = 50;
		}
     }
  }
}

void iterate(float** hotOld,float** hotNew,int start, int end) {
   int i, j;
   //#pragma omp parallel for private(i,j)
   for (i = 1; i < end-start-1; i++) {
      for (j = 1; j < NUM_OF_COLS - 1; j++) {
		if(((i+start)%20) == 0) {
		   hotNew[i][j] = 100;
		}
		else if((j%20) == 0) {
		   hotNew[i][j] = 0;
		}
		else {
		   float temp = (hotOld[i - 1][j] + hotOld[i][j+1]+hotOld[i+1][j]+hotOld[i][j-1]+4*hotOld[i][j])/8;
		   hotNew[i][j] = temp;
			}
		  }
   }
}

bool isStable(float** hot,int start, int end) {
   bool isStable = true;
   numCells = 0;
   int i, j;
   for (i = 1; i < end-start-1; i++) {
      for (j = 1; j < NUM_OF_COLS - 1; j++) {
		if (hot[i][j] > 50) {
		   numCells++;
		}
		float test = fabsf(hot[i][j] - (hot[i - 1][j] + hot[i + 1][j] + hot[i][j - 1] + hot[i][j + 1]) / 4);
		if (test > EPSILON) {
		   if ((i % 20) == 0) {}
		   else if ((j % 20) == 0) {}
		   else {
			  //printf("i: %d\n", i);
			  //printf("j: %d\n", j);
			  //printf("%f\n", test);
			  isStable = false;
		   }
	}
      }
   }
   return isStable;
}

/*return the current time in seconds, using a double precision number.*/

double When()
{
   struct timeval tp;
   gettimeofday(&tp, NULL);
   return ((double)tp.tv_sec + (double)tp.tv_usec * 1e-6);
}

