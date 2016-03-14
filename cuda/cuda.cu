#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define BLOCKSIZE 1024
#define MAXIT 1
#define TOTROWS		(BLOCKSIZE*8)
#define TOTCOLS		(BLOCKSIZE*8)
#define NOTSETLOC       -1 // for cells that are not fixed

#define QMAX(x,y) (((x) > (y))? (x): (y))


int *lkeepgoing;
float *iplate;
float *oplate;
float *fixed;
float *tmp;
int ncols, nrows;

double When();
void Compute();


int main(int argc, char *argv[])
{
	double t0, tottime;
	ncols = TOTCOLS;
	nrows = TOTROWS;

	cudaMalloc((void **) &lkeepgoing, nrows * ncols * sizeof(int));
	cudaMalloc((void **) &iplate, nrows * ncols * sizeof(float));
	cudaMalloc((void **) &oplate, nrows * ncols * sizeof(float));
	cudaMalloc((void **) &fixed,  nrows * ncols * sizeof(float));
	fprintf(stderr,"Memory allocated\n");

	t0 = When();
	/* Now proceed with the Jacobi algorithm */
	Compute();

	tottime = When() - t0;
	printf("Total Time is: %lf sec.\n", tottime);

	return 0;
}

__global__ void InitArrays(float *ip, float *op, float *fp, int *kp, int ncols)
{
	int i;
	float *fppos, *oppos, *ippos;
        int *kppos;
        int blockOffset;
        int rowStartPos;
        int colsPerThread;
	
        // Each block gets a row, each thread will fill part of a row

	// Calculate the offset of the row
        blockOffset = blockIdx.x * ncols;
        // Calculate our offset into the row
	rowStartPos = threadIdx.x * (ncols/blockDim.x);
        // The number of cols per thread
        colsPerThread = ncols/blockDim.x;

	ippos = ip + blockOffset+ rowStartPos;
	fppos = fp + blockOffset+ rowStartPos;
	oppos = op + blockOffset+ rowStartPos;
	kppos = kp + blockOffset+ rowStartPos;

	for (i = 0; i < colsPerThread; i++) {
		fppos[i] = NOTSETLOC; // Not Fixed
		ippos[i] = 50;
		oppos[i] = 50;
	        kppos[i] = 1; // Keep Going
	}
        // Insert code to set the rest of the boundary and fixed positions
}
__global__ void doCalc(float *iplate, float *oplate, int ncols)
{
	/* Compute the 5 point stencil for my region */
}

__global__ void doCheck(float *iplate, float *oplate, float *fixed, int *lkeepgoing, int ncols)
{
	// Calculate keepgoing array
}

__global__ void reduceSingle(int *idata, int *single, int nrows)
{
	// Reduce rows to the first element in each row
	int i;
    int rowStartPos;
    int colsPerThread;
    int *mypart;
    extern __shared__ int parts[];
	
    // Each block gets a row, each thread will reduce part of a row

        // Calculate our offset into the row
	rowStartPos = threadIdx.x * (nrows/blockDim.x);
        // The number of cols per thread
    colsPerThread = nrows/blockDim.x;

	//mypart = idata + blockOffset + rowStartPos;

	// sum my part of one dim array
	parts[threadIdx.x] = 0;
	for (i = rowStartPos; i < colsPerThread + start; i++) {
		parts[threadIdx.x] += idata[i];
	}
	int tid = threadIdx.x
	if(tid <512) { parts[tid] += parts[tid+512];}
	__syncthreads();
	if(tid <256) { parts[tid] += parts[tid+256];}
	__syncthreads();
	if(tid <128) { parts[tid] += parts[tid+128];}
	__syncthreads();
	if(tid <64) { parts[tid] += parts[tid+64];}
	__syncthreads();
	if(tid <32) { parts[tid] += parts[tid+32];}
	__syncthreads();
	if(tid == 0) {
		for(i=0;i<32;i++) {
			*single += parts[i];
		}
	}
}
__global__ void reduceSum(int *idata, int *odata, unsigned int ncols)
{
	// Reduce rows to the first element in each row
	int i;
    int blockOffset;
    int rowStartPos;
    int colsPerThread;
    int *mypart;
	
    // Each block gets a row, each thread will reduce part of a row

	// Calculate the offset of the row
        blockOffset = blockIdx.x * ncols;
        // Calculate our offset into the row
	rowStartPos = threadIdx.x * (ncols/blockDim.x);
        // The number of cols per thread
        colsPerThread = ncols/blockDim.x;

	mypart = idata + blockOffset + rowStartPos;

	// Sum all of the elements in my thread block and put them 
        // into the first column spot
	for (i = 1; i < colsPerThread; i++) {
		mypart[0] += mypart[i];
	}
	__syncthreads(); // Wait for everyone to complete
        // Now reduce all of the threads in my block into the first spot for my row
	if(threadIdx.x == 0) {
		odata[blockIdx.x] = 0;
		for(i = 0; i < blockDim.x; i++) {
			odata[blockIdx.x] += mypart[i*colsPerThread];
		}
	}
	// We cant synchronize between blocks, so we will have to start another kernel
}
	
void Compute()
{
	int *keepgoing_single;
	int *keepgoing_sums;
	int keepgoing;
	int blocksize = BLOCKSIZE;
	int iteration;

	ncols = TOTCOLS;
	nrows = TOTROWS;

	// One block per row
	InitArrays<<< nrows, blocksize >>>(iplate, oplate, fixed, lkeepgoing, ncols);
	cudaMalloc((void **)&keepgoing_single, 1 * sizeof(int));
	keepgoing = 1;
	cudaMalloc((void **)&keepgoing_sums, nrows * sizeof(int));
 	int *peek = (int *)malloc(nrows*sizeof(int));

	for (iteration = 0; (iteration < MAXIT) && keepgoing; iteration++)
	{
		doCalc<<< nrows, blocksize >>>(iplate, oplate, ncols);
		doCheck<<< nrows, blocksize >>>(iplate, oplate, fixed, lkeepgoing, ncols);
		//reduce value to first cell in each row`
        	reduceSum<<< nrows, blocksize>>>(lkeepgoing, keepgoing_sums, ncols);
		// cudaMemcpy(peek, keepgoing_sums, nrows*sizeof(int), cudaMemcpyDeviceToHost);
		// fprintf(stderr, "after cudaMemcpy \n");
		// int i;
 	// 	for(i = 0; i < nrows; i++) {
		// 	fprintf(stderr, "%d, ",peek[i]);
		// }
		// Now we have the sum for each row in the first column, 
		//  reduce to one value
		double t0 = When();
		int singleLoop;
		for(singleLoop = 0; singleLoop < 1000; singleLoop++) {}
			reduceSingle<<<1, blocksize, blocksize*sizeof(int)>>>(keepgoing_sums, keepgoing_single);
		}
		printf("reduce single:%f\n", When() - t0);
		keepgoing = 0;
		cudaMemcpy(&keepgoing, keepgoing_single, 1 * sizeof(int), cudaMemcpyDeviceToHost);
		fprintf(stderr, "keepgoing = %d\n", keepgoing);

		/* swap the new value pointer with the old value pointer */
		tmp = oplate;
		oplate = iplate;
		iplate = tmp;
	}
	free(peek);
	cudaFree(keepgoing_single);
	cudaFree(keepgoing_sums);
	fprintf(stderr,"Finished in %d iterations\n", iteration);
}

/* Return the current time in seconds, using a double precision number.       */
double When()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}
