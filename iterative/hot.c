#include <stdbool.h>
#include <time.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void initialize(float** hotOld, float** hotNew);

void iterate(float** hotOld, float** hotNew);

bool isStable(float** hotOld);

double When();

//#pragma omp for private(row,col)

//using namespace std;

const int NUM_OF_ROWS = 4096;
const int NUM_OF_COLS = 4096;
const float EPSILON = 0.1;

int main(int argc, char *argv[]) {
   int i, j;
   //printf("argv[1]: %s\n", argv[1]);
   //int numThreads = atoi(argv[1]);
   //printf("num of threads %d\n", numThreads);
   //omp_set_num_threads(numThreads);

   time_t start, end;
   clock_t ticks;
   long count;

   time(&start);

   double begin = When();  

   float** hotplateOld = (float**)malloc(NUM_OF_ROWS*sizeof(float*));
   for (i = 0; i < NUM_OF_ROWS; i++) {
      hotplateOld[i] = (float*)malloc(NUM_OF_COLS * sizeof(float));
   }
   float** hotplateNew = (float**)malloc(NUM_OF_ROWS * sizeof(float*));
   for (i = 0; i < NUM_OF_ROWS; i++) {
      hotplateNew[i] = (float*)malloc(NUM_OF_COLS * sizeof(float));
   }
   float** temp = (float**)malloc(NUM_OF_ROWS * sizeof(float*));
   for (i = 0; i < NUM_OF_ROWS; i++) {
      temp[i] = (float*)malloc(NUM_OF_COLS * sizeof(float));
   }
   initialize(hotplateOld,hotplateNew);
   //printf("%f", hotplateNew[0][0]);
   
   int numOfIterations = 0;
   do{
	numOfIterations++;
	iterate(hotplateOld, hotplateNew);
	//printf("iterations %d\n", numOfIterations);
	//printf("New %f\n", hotplateNew[1][1]);
	//printf("Old %f\n", hotplateOld[1][1]);

	temp = hotplateOld;
	hotplateOld = hotplateNew;
	hotplateNew = temp;
   } while (!isStable(hotplateOld));
   time(&end);

   double stop = When();

   printf("%d iterations\n", numOfIterations);
   //printf("%f seconds\n", difftime(end, start));
   printf("%f seconds\n", stop-begin);
   return 0;
}

void initialize(float** hotOld, float** hotNew) {
   int i, j;
   double start = When();
   #pragma omp parallel for private(i,j)
   for (i = 0; i < NUM_OF_ROWS; i++) {
      for (j = 0; j < NUM_OF_COLS; j++) {
	if (i == 0 || j == 0 || j == 4095) {
   	   hotOld[i][j] = 0;
	   hotNew[i][j] = 0;
	}
	else if (i == 4095) {
	   hotOld[i][j] = 100;
	   hotNew[i][j] = 100;
	}
	else if (i == 400 && j > 0 && j < 331) {
	   hotOld[i][j] = 100;
	   //hotNew[i][j] = 100;
	}
	else if (i == 200 && j == 500) {
	   hotOld[i][j] = 100;
	   //hotNew[i][j] = 100;
	}
	else {
	   hotOld[i][j] = 50;
	   //hotNew[i][j] = 50;
	}
     }
  }
  double end = When();
  printf("Initialization time: %f seconds\n",end-start);
}

void iterate(float** hotOld, float** hotNew) {
   int i, j;
   #pragma omp parallel for private(i,j)
   for (i = 1; i < NUM_OF_ROWS - 1; i++) {
      //printf("thread num; %d,%d\n", omp_get_thread_num(), i); 
      for (j = 1; j < NUM_OF_COLS - 1; j++) {
	float temp = (hotOld[i - 1][j] + hotOld[i][j+1]+hotOld[i+1][j]+hotOld[i][j-1]+4*hotOld[i][j])/8;
	hotNew[i][j] = temp;
      }
   }
   #pragma omp parallel for private(j)
   for (j = 1; j < 331; j++) {
      hotNew[400][j] = 100;
   }
   hotNew[200][500] = 100;
}

bool isStable(float** hot) {
   int numOfCells = 0;
   int i, j;
   for (i = 1; i < NUM_OF_ROWS - 1; i++) {
      for (j = 1; j < NUM_OF_COLS - 1; j++) {
	if (hot[i][j] > 50) {
	   numOfCells++;
	}
	float test = fabsf(hot[i][j] - (hot[i - 1][j] + hot[i + 1][j] + hot[i][j - 1] + hot[i][j + 1]) / 4);
	if (test > EPSILON) {
	   if (i == 200 && j == 500) {}
	   else if (i == 400 && j > 0 && j < 331) {}
	   else {
	      //printf("i: %d\n", i);
	      //printf("j: %d\n", j);
	      //printf("%f\n", test);
	      numOfCells = 0;
	      return false;
	   }
	}
      }
   }
   printf("numberOfCells: %d\n", numOfCells);
   return true;
}

/*return the current time in seconds, using a double precision number.*/

double When()
{
   struct timeval tp;
   gettimeofday(&tp, NULL);
   return ((double)tp.tv_sec + (double)tp.tv_usec * 1e-6);
}

