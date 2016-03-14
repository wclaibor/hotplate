#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pthread.h"

void *initialize(void* argument);

void *iterate(void* argument);

void *isStable(void* argument);

double When();

//#pragma omp for private(row,col)

//using namespace std;

const int NUM_OF_ROWS = 16384;
const int NUM_OF_COLS = 16384;
const float EPSILON = 0.1;

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
   int i, j;
   //printf("argv[1]: %s\n", argv[1]);
   int numThreads = atoi(argv[1]);
   //printf("num of threads %d\n", numThreads);
   //omp_set_num_threads(numThreads);

   time_t start, end;
   clock_t ticks;
   long count;

   //create threads
   int nThreads = numThreads;
   pthread_t threads[nThreads];
   long thread;
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
   int range = NUM_OF_ROWS/nThreads;
   double beginInit = When();
   for(thread=0;thread<nThreads;thread++) {
      int start = range * thread;
      int end = range * (thread+1);
      struct arg_struct *args = malloc(sizeof *args);
      args->hotOld = hotplateOld;
      args->hotNew = hotplateNew;
      args->start = start;
      args->end = end;
      pthread_create(&threads[thread], NULL, &initialize, args);
   }
   for(thread=0;thread<nThreads;thread++) {
      pthread_join(threads[thread], NULL);
   }
   double stopInit = When();
   printf("Initialization time: %f seconds\n",stopInit-beginInit);

   //printf("%f", hotplateNew[0][0]);
   int numOfIterations = 0;
   bool stable;
   int numCells;
   //for(i=0;i<98;i++){
   do {
	stable = true;
	numCells = 0;
	numOfIterations++;
	for(thread=0;thread<nThreads;thread++) {
	   int start = range * thread;
	   int end = range * (thread+1);
	   struct arg_struct *args=malloc(sizeof *args);
	   args->hotOld = hotplateOld;
	   args->hotNew = hotplateNew;
	   args->start = start;
	   args->end = end;
	   pthread_create(&threads[thread], NULL, &iterate,args); 
        }
	for(thread=0;thread<nThreads;thread++) {
	   pthread_join(threads[thread], NULL);
        }
	printf("iterations %d\n", numOfIterations);
	temp = hotplateOld;
	hotplateOld = hotplateNew;
	hotplateNew = temp;
	for(thread=0;thread<nThreads;thread++) {
	   int start = range * thread;
	   int end = range * (thread+1);
	   //printf("start: %d\nend: %d\n", start,end);
	   struct arg_struct *args=malloc(sizeof *args);
	   args->hotOld = hotplateOld;
	   args->hotNew = hotplateNew;
	   args->start = start;
	   args->end = end;
	   pthread_create(&threads[thread], NULL, &isStable,args); 
	}
	for(thread=0;thread<nThreads;thread++) {
	   struct ret_struct *stat;
	   pthread_join(threads[thread], (void*)&stat);
	   
	   numCells += stat->numCells;
	   if(stat->isStable == false) {
	      stable = false;
	   }
           free(stat);
        }
   } while (!stable);
   time(&end);

   double stop = When();
   printf("%d number of cells\n", numCells);
   printf("%d iterations\n", numOfIterations);
   //printf("%f seconds\n", difftime(end, start));
   printf("%f seconds\n", stop-begin);
   return 0;
}

void *initialize(void * argument) {
   int i, j;
   //#pragma omp parallel for private(i,j)
   struct arg_struct *args = argument;
   float** hotNew = args->hotNew;
   float** hotOld = args->hotOld;
   int start = args->start;
   int end = args->end;
   for (i = start; i < end; i++) {
      for (j = 0; j < NUM_OF_COLS; j++) {
	if (j == NUM_OF_COLS - 1) {
   	   hotOld[i][j] = 0;
	   hotNew[i][j] = 0;
	}
	else if (i == NUM_OF_ROWS - 1) {
	   hotOld[i][j] = 100;
	   hotNew[i][j] = 100;
	}
	else if ((i%20) == 0) {
	   hotOld[i][j] = 100;
	   hotNew[i][j] = 100;
	}
	else if((j%20) == 0) {
	   hotOld[i][j] = 0;
	   hotNew[i][j] = 0;
	}
	else {
	   hotOld[i][j] = 50;
	   //hotNew[i][j] = 50;
	}
     }
  }
  free(args);
}

void *iterate(void * argument) {
   struct arg_struct *args = argument;
   float** hotNew = args->hotNew;
   float** hotOld = args->hotOld;
   int start = args->start;
   if(start == 0) start = 1;
   int end = args->end;
   if(end >= NUM_OF_ROWS) end = NUM_OF_ROWS-1;
   int i, j;
   //#pragma omp parallel for private(i,j)
   for (i = start; i < end; i++) {
      for (j = 1; j < NUM_OF_COLS - 1; j++) {
	if((i%20) == 0) {
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
   free(args);
}

void *isStable(void* argument) {
   struct arg_struct *args = argument;
   float** hot = args->hotNew;
   int start = args->start;
   if(start == 0) start = 1;
   int end = args->end;
   if(end >= NUM_OF_ROWS) end = NUM_OF_ROWS-1;
   bool isStable = true;
   int numOfCells = 0;
   int i, j;
   for (i = start; i < end; i++) {
      for (j = 1; j < NUM_OF_COLS - 1; j++) {
	if (hot[i][j] > 50) {
	   numOfCells++;
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
   struct ret_struct *ret = malloc(sizeof *ret);
   ret->numCells = numOfCells;
   ret->isStable = isStable;
   pthread_exit(ret);
}

/*return the current time in seconds, using a double precision number.*/

double When()
{
   struct timeval tp;
   gettimeofday(&tp, NULL);
   return ((double)tp.tv_sec + (double)tp.tv_usec * 1e-6);
}

