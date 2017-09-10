#include <fstream>
#include <vector>
#include <iostream>
#include <random>

#include <mpi.h>

#include "utility.hpp"


int main(int argc, char **argv){
    MPI_Init(&argc, &argv);

    int id, ier, num_p;
	
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_p);

	double* matrix;
    double* vector;
    double* prevvector;
    double* out;
	int n, m;
    double diff = 1e10;
	double eps = 1e-7;

	int iter = 0;
    while(diff > eps){
		if (id == 0){
			if(iter == 0){		
				eps = 1e-7;

				//n = 7, m = 7;

				//Read matrix from file name t1
				//std::ifstream in_matrix;
				//std::ifstream in_vector;
	
				//in_matrix.open("t1");
				//
				//in_matrix >> n >> m;
                                
                                std::cin >> n;
                                m = n;

				M = new double[m*n];
				V = new double[n];
				prevV = new double[n];
				out = new double[n];

                                make_sth_matrix(M, n, n);

				//for(int i = 0; i < n; ++i){
				//	for(int j = 0; j < n; ++j){
				//		in_matrix >> M[i*n+j];
				//	}
				//}
				//End read
	

				//make_sth_matrix(M, n, m);
				//print(n, m, M);

    	    	for(int i = 0; i < m; ++i){
    	    		vector[i] = ((double) std::rand())/RAND_MAX;
    	    		
    	    	}

				norming(m, vector);

    	    	for(int i = 0; i < n; ++i){
    	    		prevvector[i] = 1e10;
    	    	}
			}


			int num_my_vectors = n % (num_p-1);
			int num_vectors_per_process = n/(num_p-1);
    	    for(int i = 1; i < num_p; ++i){
    	        MPI_Send(&num_vectors_per_process,
			    	1, MPI_INT, i, 0, MPI_COMM_WORLD);
    	        MPI_Send(&n,
			    	1, MPI_INT, i, 1, MPI_COMM_WORLD);
    	        MPI_Send(matrix + (num_my_vectors + (i-1) * num_vectors_per_process) * n,
			    	n * num_vectors_per_process, MPI_DOUBLE_PRECISION, 
    	            	i, 2, MPI_COMM_WORLD);
    	        MPI_Send(vector, n, MPI_DOUBLE_PRECISION, 
    	                i, 3, MPI_COMM_WORLD);
    	    }

			for(int i = 0; i < num_my_vectors; ++i){
				out[i] = 0;
				for(int j = 0; j < m; ++j){
				    out[i] += matrix[i*n + j]*vector[j]; 
				}
			}

			for(int i = 1; i < num_p; ++i){
				MPI_Recv(out + num_my_vectors + (i-1) * num_vectors_per_process, 
					num_vectors_per_process, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD,
							MPI_STATUS_IGNORE);
			}

    	    norming(n, out);
    	    
    	    std::copy(vector, vector + n, prevvector);
    	    std::copy(out, out + n, vector);
    	    
    	    diff = mean_diff(m, vector, prevvector);
			for(int i = 1; i < num_p; ++i){
    	    	MPI_Send(&diff, 1, MPI_DOUBLE_PRECISION, 
    	    	        i, 4, MPI_COMM_WORLD);
			}
			++iter;	
		}else{
			int num_vectors_per_process = 0;
			int n = 0;
    	    MPI_Recv(&num_vectors_per_process,
				1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    	    MPI_Recv(&n,
				1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    	    double lmatrix[n*num_vectors_per_process];
    	    double lvector[n];
    	    MPI_Recv(lmatrix, n*num_vectors_per_process, 
				MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    	    MPI_Recv(lvector, n, MPI_DOUBLE_PRECISION, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    	    

    	    double out[num_vectors_per_process];
			for(int i = 0; i < num_vectors_per_process; ++i){
				out[i] = 0;
				for(int j = 0; j < n; ++j){
				    out[i] += lmatrix[i*n + j]*lvector[j]; 
				}
			}
			MPI_Send(out, num_vectors_per_process, MPI_DOUBLE_PRECISION, 0, id, MPI_COMM_WORLD);
    	    MPI_Recv(&diff, 1, MPI_DOUBLE_PRECISION, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	
	if (id == 0){
		std::cout << "out" << std::endl;
		print(n, out);

		delete [] matrix;
		delete [] vector;
		delete [] out;
		delete [] prevvector;
	} MPI_Finalize(); 
	
}
