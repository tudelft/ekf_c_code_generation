#include "maths.h"
#include <math.h>

// columns major. upper factor
void chol(float *U, float *A, float *iDiag, int n)
{
    // rosetta code
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < (i+1); j++) {
            float s = 0;
            for (int k = 0; k < j; k++) {
                s += U[i * n + k] * U[j * n + k];
            }
            if (i == j) {
                U[i * n + j] = sqrt(A[i * n + i] - s);
                iDiag[j] = 1.f / U[i * n + j];
            } else {
                U[i * n + j] = iDiag[j] * (A[i * n + j] - s);
            }
        }
    }
}

// column major, upper factor
void chol_solve(float *U, float* iDiag, int n, float *b, float *x) {
    // Antoine Drouin, 2007, modified
	int j,k;
	float t;

    for(j = 0 ; j < n ; j++) { // solve Uty=b
        t = b[j];
        for(k = j - 1 ; k >= 0 ; k--)
            t -= U[k + n*j] * x[k];
        x[j] = t*iDiag[j];
    }
    for(j = n - 1 ; j >= 0 ; j--) { // solve Ux=y
        t = x[j];
        for(k = j + 1 ; k < n ; k++)
            t -= U[j + n*k] * x[k];
        x[j] = t*iDiag[j];
    }
}
