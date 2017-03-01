#include "summa.h"
#include "matrix_op.h"
#include <mpi.h>
#include <matrix_op.h>

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
    char transfer_A = 'N';
    char transfer_B = 'N';

    dgemm_(&transfer_B, &transfer_A, &n, &m,&k, &alpha, B, &n, A, &k, &beta, C,&n);
}

// verify results usig blas
int verify(double *C, double *C_blas, int n){
    int i, j;

    for(i = 0; i < n; i++){
        for(j = 0; j< n; j++){
            if(abs(C[i*n +j] -  C_blas[i*n+j]) > 0.000001){
                printf("element(%d, %d)%3.12lf != %3.12lf\n",i,j, C[i*n +j], C_blas[i*n +j]);
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
// start  prepare mpi
    int nprocs, rank;
    MPI_Comm gcomm;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);
    gcomm = MPI_COMM_WORLD;

    double *C, *B, *A, *C_blas;

    if(argc !=2){
        printf("mpirun run_summa N\n");
    }


 // block configruation end of MPI prepare

    int N = atoi(argv[1]);
    int procs_per_side = sqrt(nprocs);

    // submatrix size
    int sub_N = N/sqrt(nprocs); 

    // cartisian pos
    int my_row = rank/procs_per_side;
    int my_col = rank%procs_per_side;

    MPI_Comm row_comm, col_comm;

    MPI_Comm_split(MPI_COMM_WORLD, my_row, rank, &row_comm);
    MPI_Comm_split(MPI_COMM_WORLD, my_col, rank, &col_comm);

    int i, j,k, retval;
    j = 0;

    // block matrix multiplication
    int b;


    int exp;

    for(b = 1; b<=sub_N; b*=4){

        double time_sum = 0;
        
        for(exp =0; exp <NUM_EXP; exp++){

            // this part is only avaible in rank 0
            if(rank == 0){
                double *C, *A, *B, *C_blas;
                double btime, etime;
                // configurations
                double total_time = 0;

                if(N >0){
                    printf("matrix all have dimension %d x %d\n", N, N);
                    printf("totally %d procs, submatrix_N = %d, block_size = %d \n", nprocs, sub_N, b);
                }
             // generate three input matrix
                C = NULL;
                B = NULL;
                A = NULL;
                C_blas =NULL;
            }

            //printf("I am procs %d, row is %d col is %d\n", rank, my_row, my_col);


            // prepare local space for blocks
            double *local_C, *local_A, *local_B;
            init_matrix(&local_C, sub_N, 0);
            init_matrix(&local_A, sub_N, 0);
            init_matrix(&local_B, sub_N, 0);

            // also each process will have a row vector and col vector
            double *local_col =(double*)malloc(sizeof(double)*sub_N*b);
            double *local_row =(double*)malloc(sizeof(double)*sub_N*b);

            // Step1
            //  only rank 0 need to prepare all the matrix
            if(rank == 0){
                double *tmp_buffer_A;
                double *tmp_buffer_B;
                int ii, jj;
                init_matrix(&tmp_buffer_A, sub_N, 0);
                init_matrix(&tmp_buffer_B, sub_N, 0);

                printf("****************************************************\n");
                printf("start No.%d experiment \n", exp);
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

                init_matrix(&C, N, 0);
                int target=0;
                for(i = 0; i < N; i+= sub_N){
                    for(j = 0; j < N; j+=sub_N){
                        int cur = 0;
                        if(target == 0){
                        //start inner block multiplication
                            for(ii = 0; ii < sub_N; ii++){
                                    for(jj = 0; jj < sub_N ;jj+=1){
                                             // put into loacal buffer
                                             local_A[cur] = A[(i+ii)*N + j+ jj];
                                             local_B[cur] = B[(i+ii)*N + j+ jj];
                                             cur +=1;
                                    }
                            }
                        }
                        else{

                            //start inner block multiplication
                            for(ii = 0; ii < sub_N; ii++){
                                    for(jj = 0; jj < sub_N ;jj+=1){
                                             // put into loacal buffer
                                             tmp_buffer_A[cur] = A[(i+ii)*N + j+jj];
                                             tmp_buffer_B[cur] = B[(i+ii)*N + j+jj];
                                             cur +=1;
                                    }
                            }
                        // end of inner multiplication
                         MPI_Send(tmp_buffer_A, sub_N*sub_N, MPI_DOUBLE, target, 1, gcomm);
                         MPI_Send(tmp_buffer_B, sub_N*sub_N, MPI_DOUBLE, target, 2, gcomm);
                        }
                    target++;
                    }
                }
            }
            else{

                MPI_Recv(local_A, sub_N*sub_N, MPI_DOUBLE, 0, 1, gcomm,MPI_STATUS_IGNORE);
                MPI_Recv(local_B, sub_N*sub_N, MPI_DOUBLE, 0, 2, gcomm, MPI_STATUS_IGNORE);
            }
            // receive block

            // start barrier here
            MPI_Barrier(gcomm);
            if(rank == 0){
                printf("all blocks are distributed\n");
            }
            double t1 = MPI_Wtime();

            /*
             * start main cal
             */
            for(k = 0; k < N/b; k+=1){
                // check wheter I have the col
                if(my_col == k*b/sub_N){
                    // prepare the col
                    //get_col(local_A, local_col,sub_N, (k*b)%sub_N, b);
                    int local_block_index = k%(sub_N/b);
                    //printf("I am procs %d, colblock of index%d\n", rank, local_block_index);
                    get_col(local_A, local_col,sub_N, local_block_index, b);
                }
                MPI_Bcast(local_col, sub_N*b, MPI_DOUBLE, k*b/sub_N/*root*/, row_comm);

                // check whether I have this row(B) and then send to other proc in the same col
                if(my_row == k*b/sub_N){

                    //printf("I am procs %d, row is %d col is %d\n", rank, my_row, my_col);
                    int local_block_index = k%(sub_N/b);
                    //printf("I am procs %d, rowblock of index%d\n", rank, local_block_index);
                    get_row(local_B, local_row, sub_N, local_block_index, b);
                }
                MPI_Bcast(local_row, sub_N*b, MPI_DOUBLE, k*b/sub_N/*root*/, col_comm);

                // do local calculation
                update_local_C(local_C, local_col, local_row, sub_N, sub_N ,b);

            }

            double t2 = MPI_Wtime();
            double time_comp = t2 -t1;

            // barrier here
            MPI_Barrier(gcomm);
            //double global_time_comm;
            double global_time_comp;        
            //MPI_Reduce(&time_comm, &global_time_comm, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm); 
            MPI_Reduce(&time_comp, &global_time_comp, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm); 

            // Print the result    
            if (rank == 0) {       
              printf("exp_%d: Computation Total %lf avg %lf\n",exp,  global_time_comp , global_time_comp/ (nprocs));
              time_sum+=global_time_comp/(nprocs);
            }


            if(rank == 0){
                printf("calculation completed\n");
            }

            // start gathering
            if(rank == 0){
                double * tmp_buffer_C;
                init_matrix(&tmp_buffer_C, sub_N, 0);
                int ii, jj;
                int source = 0;
                for(i = 0; i < N; i+= sub_N){
                    for(j = 0; j < N; j+=sub_N){
                        int cur = 0;
                        if(source == 0){
                        //start inner block multiplication
                            for(ii = 0; ii < sub_N; ii++){
                                    for(jj = 0; jj < sub_N ;jj+=1){
                                             // put into loacal buffer
                                             C[(i+ii)*N + j + jj] = local_C[cur];
                                             cur +=1;
                                    }
                            }
                        }
                        else{

                         MPI_Recv(tmp_buffer_C, sub_N*sub_N, MPI_DOUBLE, source, 3, gcomm,MPI_STATUS_IGNORE);

                            //start inner block multiplication
                            for(ii = 0; ii < sub_N; ii++){
                                    for(jj = 0; jj < sub_N ;jj+=1){
                                             // put into loacal buffer
                                             C[(i+ii)*N + j+jj] = tmp_buffer_C[cur];
                                             cur +=1;
                                    }
                            }
                        // end of inner multiplication
                        }
                        source++;
                    }
                }
               }
            else{
                MPI_Send(local_C, sub_N*sub_N, MPI_DOUBLE, 0, 3, gcomm);
            }
            // barrier here
            MPI_Barrier(gcomm);

            if(rank == 0){
                printf("results gathered\n");
            }



            // combine all the blocks to C

            // now free local buffer
            free(local_C);
            free(local_A);
            free(local_B);


            // end calculation here, rank 0 will verify add record time
            if(rank ==0){

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
            }
        }// finish current experiment
    double n = N;
    if (rank == 0) {       
          double avg_time = time_sum/NUM_EXP;
          double scalar = 2*n*n*n*(1E-9);
          double avg_gflops = scalar/avg_time;
          fprintf(stderr, "blk_size: %d,  computation time in %d exp:  %lf, gflops %f, gflops perprocess %f\n", b, NUM_EXP, avg_time, avg_gflops, avg_gflops/nprocs);

        }
    }
                
    MPI_Barrier(gcomm);
    MPI_Finalize();
    return 0;

}


    


   

