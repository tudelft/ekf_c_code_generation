// columns major. upper factor
void chol(float *U, float *A, float *iDiag, int n);

// Column major upper factor. solve A x = UT U x = b
void chol_solve(float *U, float* iDiag, int n, float *b, float *x);
