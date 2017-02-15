#include "run_unroll.h"

// preprocessor for unrolling
// unroll_0 means 2^0 loop size, no unrolling
// unroll_3 means 2^3 loop size, using size-8 loop
#define UNROLL_1(N) \
    C[(i+ii)*n +j+jj +(N) ] += tmp* B[(k+kk)*n + j+ jj + (N)];

#define UNROLL_2(N) \
    UNROLL_1(2*(N)) \
    UNROLL_1(2*(N)+1)

#define UNROLL_4(N) \
    UNROLL_2(2*(N)) \
    UNROLL_2(2*(N) +1)

#define UNROLL_8(N) \
    UNROLL_4(2*(N)) \
    UNROLL_4(2*(N) +1)

#define UNROLL_16(N) \
    UNROLL_8(2*(N)) \
    UNROLL_8(2*(N) +1)

#define UNROLL_32(N) \
    UNROLL_16(2*(N)) \
    UNROLL_16(2*(N) +1)


#define UNROLL_64(N) \
    UNROLL_32(2*(N)) \
    UNROLL_32(2*(N) +1)








// actual 
// you can change block size and loop unrolling length
void seq_cal_ikj_block(double *C, double *A, double *B, int n, int b){
    int i, j, k, ii, jj, kk;

    if(n%b){
        printf("warning: b should be divisible by n\n");
        exit(-1);
    }

    for(i = 0; i < n; i+= b){
        for(k = 0; k < n; k+=b){
            for(j = 0; j < n; j+=b){
                // C(i,j) = C(i,j)+ A(i,k)*B(k,j)

                //start inner block multiplication
                for(ii = 0; ii < b; ii++){
                    for(kk = 0; kk < b; kk++ ){
                        for(jj = 0; jj < b ;jj++){
                            C[(i+ii)*n +j+jj] += A[(i+ii)*n + k + kk]* B[(k+kk)*n + j+ jj];
                        }
                    }
                }
                // end of inner multiplication
            }
        }
    }
}

void seq_cal_ikj_unroll_2(double *C, double *A, double *B, int n, int b){
    int i, j, k, ii, jj, kk;
    if(n%b){
        printf("warning: b should be divisible by n\n");
        exit(-1);
    }
    register float tmp;
    
    for(i = 0; i < n; i+= b){
        for(k = 0; k < n; k+=b){
            for(j = 0; j < n; j+=b){
                // C(i,j) = C(i,j)+ A(i,k)*B(k,j)

                //start inner block multiplication
                for(ii = 0; ii < b; ii++){
                    for(kk = 0; kk < b; kk++ ){
                        for(jj = 0; jj < b ;jj+=2){
                                tmp = A[(i+ii)*n + k + kk];
                                UNROLL_2(0)
                        }
                    }
                }
                // end of inner multiplication
            }
        }
    }
}

void seq_cal_ikj_unroll_4(double *C, double *A, double *B, int n, int b){
    int i, j, k, ii, jj, kk;

    if(n%b){
        printf("warning: b should be divisible by n\n");
        exit(-1);
    }
    register float tmp;

    for(i = 0; i < n; i+= b){
        for(k = 0; k < n; k+=b){
            for(j = 0; j < n; j+=b){
                // C(i,j) = C(i,j)+ A(i,k)*B(k,j)

                //start inner block multiplication
                for(ii = 0; ii < b; ii++){
                    for(kk = 0; kk < b; kk++ ){
                        for(jj = 0; jj < b ;jj+=4){
                                tmp = A[(i+ii)*n + k + kk];
                                UNROLL_4(0)
                                                        }
                    }
                }
                // end of inner multiplication
            }
        }
    }
}

void seq_cal_ikj_unroll_8(double *C, double *A, double *B, int n, int b){
    int i, j, k, ii, jj, kk;

    if(n%b){
        printf("warning: b should be divisible by n\n");
        exit(-1);
    }
    register float tmp;

    for(i = 0; i < n; i+= b){
        for(k = 0; k < n; k+=b){
            for(j = 0; j < n; j+=b){
                // C(i,j) = C(i,j)+ A(i,k)*B(k,j)

                //start inner block multiplication
                for(ii = 0; ii < b; ii++){
                    for(kk = 0; kk < b; kk++ ){
                        for(jj = 0; jj < b ;jj+=8){
                                
                                tmp = A[(i+ii)*n + k + kk];
                                UNROLL_8(0)
                        }
                    }
                }
                // end of inner multiplication
            }
        }
    }
}
void seq_cal_ikj_unroll_16(double *C, double *A, double *B, int n, int b){
    int i, j, k, ii, jj, kk;

    if(n%b){
        printf("warning: b should be divisible by n\n");
        exit(-1);
    }
    register float tmp;

    for(i = 0; i < n; i+= b){
        for(k = 0; k < n; k+=b){
            for(j = 0; j < n; j+=b){
                // C(i,j) = C(i,j)+ A(i,k)*B(k,j)

                //start inner block multiplication
                for(ii = 0; ii < b; ii++){
                    for(kk = 0; kk < b; kk++ ){
                        for(jj = 0; jj < b ;jj+=16){
                                tmp = A[(i+ii)*n + k + kk];
                                UNROLL_16(0)
                        }
                    }
                }
                // end of inner multiplication
            }
        }
    }
}

void seq_cal_ikj_unroll_32(double *C, double *A, double *B, int n, int b){
    int i, j, k, ii, jj, kk;

    if(n%b){
        printf("warning: b should be divisible by n\n");
        exit(-1);
    }
    register float tmp;

    for(i = 0; i < n; i+= b){
        for(k = 0; k < n; k+=b){
            for(j = 0; j < n; j+=b){
                // C(i,j) = C(i,j)+ A(i,k)*B(k,j)

                //start inner block multiplication
                for(ii = 0; ii < b; ii++){
                    for(kk = 0; kk < b; kk++ ){
                        for(jj = 0; jj < b ;jj+=32){
                                tmp = A[(i+ii)*n + k + kk];
                                UNROLL_32(0)
                        }
                    }
                }
                // end of inner multiplication
            }
        }
    }
}

void seq_cal_ikj_unroll_64(double *C, double *A, double *B, int n, int b){
    int i, j, k, ii, jj, kk;

    if(n%b){
        printf("warning: b should be divisible by n\n");
        exit(-1);
    }
    register float tmp;

    for(i = 0; i < n; i+= b){
        for(k = 0; k < n; k+=b){
            for(j = 0; j < n; j+=b){
                // C(i,j) = C(i,j)+ A(i,k)*B(k,j)

                //start inner block multiplication
                for(ii = 0; ii < b; ii++){
                    for(kk = 0; kk < b; kk++ ){
                        for(jj = 0; jj < b ;jj+=64){
                                tmp = A[(i+ii)*n + k + kk];
                                UNROLL_64(0)
                        }
                    }
                }
                // end of inner multiplication
            }
        }
    }
}




/*
 * generate a random matrix with double value -1~1
 * init_type
 *  1: init with random number (-1~1)
 *  0: init with 0
 */
int init_matrix(double ** m, int n, int init_type){

    

    int i, j;
    double tmp;
    srand(time(NULL));
    // allocate space
    *m = (double*)malloc(sizeof(double)*n*n);
    if(*m ==NULL){
        perror("allocate error");
        exit(-1);
    }

    // initiate with random number
    for( i = 0; i< n; i++){
        for(j = 0; j < n ; j++){

            if(init_type == 1){
                tmp = (double)rand()/RAND_MAX*2.0-1.0;
            }
            else{
                tmp = 0;
            }
            (*m)[i*n + j] = tmp;
        }
    }
    return 0;
}

void free_matrix(double *m){
    if(m != NULL){
        free(m);
    }
}

void print_matrix(double *m, int n){
    int i, j;
    for(i = 0 ; i < n ; i++){
        for(j = 0; j< n; j++){
            printf("%4.3lf ", m[i*n +j]);
        }
        printf("\n");
    }
}


// return 0 if two matrix are equal
/*
int compare_matrix(double *A, double *B, int n){
    int i, j;
    int equal =0;
    for(i = 0; i < n; i++){
        for(j = 0; j < n ; j++){
            if(A[i*n +j]!= B[i*n+j]){
                return -1;
            }
        }
    }
    return 0;
}
*/
void blas_cal(double *C, double* A, double *B, int n){
    int m, k;
    m = n;
    k = n;

    double  alpha =1;
    double  beta = 0;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                            m, n, k, alpha, A, k, B, n, beta, C, n);
}

// verify results usig blas
int verify(double *C, double *C_blas, int n){
    int i, j;

    for(i = 0; i < n; i++){
        for(j = 0; j< n; j++){
            if(abs(C[i*n +j] -  C_blas[i*n+j]) > 0.000001){
                printf("%3.12lf != %3.12lf\n", C[i*n +j], C_blas[i*n +j]);
                return -1;
            }
        }
    }
    return 0;
}



/*
 * main routine
 * input:
 *  N:matrix size
 * return:
 *  0 if success
 */
int main(int argc, char * argv[]){
    int N, i, j, retval;
    double *C, *A, *B, *C_blas;

    // block matrix multiplication
    int block_size;

    double btime, etime;


    // configurations
    int unroll_level[] = {1,2,4,8,16,32,64};

    double total_time = 0;
    // l1 data cache access and L1 data cache miss
    long long count_l1_a = 0;
    long long count_l1_m = 0;
    long long count_l2_a = 0;
    long long count_l2_m = 0;
    long long count_tlb_m = 0;

    int Events[NUM_EVENTS]={PAPI_L1_DCA, PAPI_L1_DCM,PAPI_L2_DCA, PAPI_L2_DCM,PAPI_TLB_DM};
    int EventSet = PAPI_NULL;
    long long values[NUM_EVENTS];


    if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) 
        exit(1);

    // create event set
    retval = PAPI_create_eventset (&EventSet);
    /* Add events to the eventset */
    retval = PAPI_add_events (EventSet, Events, NUM_EVENTS);


    // matrix size
    if(argc != 4){
        printf("format: aprun seq_matrix_mul N block_size loopsize, now exit\n");
        exit(-1);
    }


    N = atoi(argv[1]);
    if(N >0){
        printf("matrix all have dimension %d x %d\n", N, N);
    }
    block_size=atoi(argv[2]);
    i = atoi(argv[3]);

    void (*unroll[10])(double*, double*, double*, int, int);                                                                                                                                                                                       
    unroll[0] = &seq_cal_ikj_block;
    unroll[1] = &seq_cal_ikj_unroll_2;
    unroll[2] = &seq_cal_ikj_unroll_4;
    unroll[3] = &seq_cal_ikj_unroll_8;
    unroll[4] = &seq_cal_ikj_unroll_16;
    unroll[5] = &seq_cal_ikj_unroll_32;
    unroll[6] = &seq_cal_ikj_unroll_64;

   
    // generate three input matrix
    C = NULL;
    B = NULL;
    A = NULL;
    C_blas =NULL;

    for(j = 0 ; j < NUM_EXP ; j ++){
        printf("****************************************************\n");
        printf("start No.%d experiment \n", j);

        
        //print_matrix(C_blas, N);
        //printf("\t\tElapsed Time: %16.9f second\n", etime - btime);
        

        // can try different block size here
        //for(i = 0; i < 4; i++){
        //for(i = 3; i >=0; i--){
            if(init_matrix(&B, N, 1) == 0){
                printf("\trandom matrix B is generated \n");
            }
            if(init_matrix(&A, N, 1) == 0){
                printf("\trandom matrix A is generated \n");
            }

            // start blas_version

            if(init_matrix(&C_blas, N, 0) == 0){
                printf("\tstart blas  version...\n");
            }

            blas_cal(C_blas, A, B, N);


            printf("start unrolling level %d... \n", i);


            init_matrix(&C, N, 0);

            // init counter
            btime = get_cur_time(); 
            retval = PAPI_start (EventSet);

            //seq_cal_ikj_unroll(C,A,B,N,block_size, unroll_level[i]);
            //unroll[i](C,A,B,N,block_size);
            unroll[i](C,A,B,N,block_size);
            //seq_cal_ikj_unroll(C,A,B,N,  block_size, i);



            // get counter
            retval = PAPI_stop (EventSet, values);
            etime = get_cur_time(); 
            //print_matrix(C, N);
            
            // sum to get avg
            total_time += etime-btime;
            count_l1_a += values[0];
            count_l1_m += values[1];
            count_l2_a += values[2];
            count_l2_m += values[3];
            count_tlb_m += values[4];



            // use cray libsci(blas to verify)
            if(verify(C, C_blas, N) == 0){

//#ifdef VERBOSE
                //printf("\t\tElapsed Time: %16.9f second, %Ld miss in %Ld acess, mr=%lf\t, tlb = %Ld\t", etime - btime, values[1],values[0], (values[1]+0.0)/values[0], values[2]);
                printf("correct\n");
//#endif
            }else{
                free_matrix(C);
                free_matrix(A);
                free_matrix(B);
                free_matrix(C_blas);

                printf("calculation not correct, now exit\n");
                return -1;
            }
            free_matrix(C);

        free_matrix(A);
        free_matrix(B);
        free_matrix(C_blas);

        printf("\tall matrix are freed\n");
        //}
         
    }
    // print average
    printf("the average cacluation time of %d runs:\n", NUM_EXP);

    //fprintf(stderr, "N = %d", N);

    //fprintf(stderr, "\nblk_size\tunroll_level\tgflops\tl1_a\tl1_m\tL1MR\tl2_a\tl2_m\tL2MR\tTLBM");

    // floating point overflow
    double n = N;
    //for(i = 0; i <  4; i ++){
        double avg_time = total_time/NUM_EXP;
        double scalar = 2*n*n*n*(1E-9);
        printf("scalar = %lf\n", scalar);
        double avg_gflops = scalar/avg_time;
        double avg_l1_a = count_l1_a/NUM_EXP;
        double avg_l1_m = count_l1_m/NUM_EXP;
        double avg_l2_a = count_l2_a/NUM_EXP;
        double avg_l2_m = count_l2_m/NUM_EXP;
        double avg_l1_mr = avg_l1_m/avg_l1_a;
        double avg_l2_mr = avg_l2_m/avg_l2_a;
        double avg_tlb_m = count_tlb_m/NUM_EXP;

        printf("unroll level %d:\t %lf s, %lf gflops, L1 data access %lf , L1 data miss %lf, miss rate %lf, tlb mis %lf\n", unroll_level[i],avg_time, avg_gflops, avg_l1_a, avg_l1_m, avg_l1_mr, avg_tlb_m);
        fprintf(stderr, "\n%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",block_size, unroll_level[i], avg_gflops, avg_l1_a, avg_l1_m, avg_l1_mr, avg_l2_a, avg_l2_m, avg_l2_mr,avg_tlb_m);
    //}

    return 0;

}


    


   

