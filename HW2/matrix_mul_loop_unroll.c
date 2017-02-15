void seq_cal_ikj_loop(double *C, double *A, double *B, int n, int l){
    if(n%l || n%b){
        printf("l and b should be divisible by n\n");
        exit(-1);
    }
    int i, j, k;
    //double cij;
    for(i = 0; i< n; i++){
        for(k = 0; k< n;k++){
            for(j = 0; j < n; j += l){
                /*
                cij = C[i*n + j];
                cij += A[i*n +k]* B[k*n+ j]; // cij +=A[i][k]+ B[k][j]
                C[i*n +j] = cij;
                */
                if(l == 1){
                    C[i*n +j] += A[i*n +k]* B[k*n+ j];
                }
                else if(l == 2){
                    C[i*n +j] += A[i*n +k]* B[k*n+ j];
                    C[i*n +j+1] += A[i*n +k]* B[k*n+ j+1];
                }
                else if(l == 4){
                    C[i*n +j] += A[i*n +k]* B[k*n+ j];
                    C[i*n +j+1] += A[i*n +k]* B[k*n+ j+1];
                    C[i*n +j+2] += A[i*n +k]* B[k*n+ j+2];
                    C[i*n +j+3] += A[i*n +k]* B[k*n+ j+3];
                }
                else if(i == 8){
                    C[i*n +j] += A[i*n +k]* B[k*n+ j];
                    C[i*n +j+1] += A[i*n +k]* B[k*n+ j+1];
                    C[i*n +j+2] += A[i*n +k]* B[k*n+ j+2];
                    C[i*n +j+3] += A[i*n +k]* B[k*n+ j+3];
                    C[i*n +j+4] += A[i*n +k]* B[k*n+ j+4];
                    C[i*n +j+5] += A[i*n +k]* B[k*n+ j+5];
                    C[i*n +j+6] += A[i*n +k]* B[k*n+ j+6];
                    C[i*n +j+7] += A[i*n +k]* B[k*n+ j+7];
                }
                else{
                    printf("loop size not supported\n");
                    exit(-1);
                }
            }
        }
    }
}
