#include "R.h"
#include "Rinternals.h"


// Function createSamplePath
// This function return an array of paired node number 
// of size (size). The array returned is the one 
// specified by the offset
// 
// Arguments:
// total_n        = number of nodes in g1
// offset (1:...) = starting position in pathways
//                   note: the total number of offset is ((total_n*(total_n-1)/2)/size)+1 --undirected
//                  total number of offset is ((total_n*(total_n-1))/size)+1
// size           = size of sampling (returned array)
//
// Return
// And 
void createSamplePath(int *total_n, int *offset, int *size, int* path, int *directed){
int i, j, k, t,stotal,ssize, soffset,sdirected;
ssize=*size;
soffset=*offset;
stotal=*total_n;
sdirected=*directed;
k=0;
t=0;
//Rprintf("%i %i %i", stotal, soffset, ssize);
//PROTECT(z = allocVector(INTSXP, ssize));
 //--Undirected
 if (sdirected==0) {
	 for(j = 0; j < stotal; j++) {
		for(i = j; i < stotal; i++){
			
			if (i>j) {			
				if (((t/ssize)+1)==soffset) {			
					path[k] =i+1;
					path[k+1] =j+1;
					k=k+2;
				}
				t++;
			}
			if (((t/ssize)+1)>soffset) break;
		}
	}
 } else {
	  for(j = 0; j < stotal; j++) {
		for(i = 0; i < stotal; i++){
			
			if (i!=j) {			
				if (((t/ssize)+1)==soffset) {			
					path[k] =i+1;
					path[k+1] =j+1;
					k=k+2;
				}
				t++;
			}
			if (((t/ssize)+1)>soffset) break;
		}
	}	 
 }
//UNPROTECT(1); /*z*/
//return z;
}
