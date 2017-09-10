#include <fstream>
#include <vector>
#include <iostream>
#include <random>
#include <chrono>

#include <mpi.h>

template<class T>
void print(int n, T v){
	for(int i = 0; i < n; ++i){
		////std::cout << v[i] << " ";
	}
}

template<class T>
void print(T v, int n, int m){
    for(int i = 0; i < n; ++i){
	for(int j = 0; j < m; ++j){
		////std::cout << v[i*n+j] << " ";
	}
        ////std::cout << std::endl;
    }
}

template<class T>
void make_sth_matrix(T* M, int n, int m){
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < m; ++j){
			int idx = (int) (n*((double) std::rand())/RAND_MAX);
			M[idx*n + i] += 1.0/n;
		}	
	}
}

template<class T>
void matrix_dot_vector(int n, int m, T out, T M, T V){
	for(int i = 0; i < n; ++i){
		out[i] = 0;
		for(int j = 0; j < m; ++j){
			out[i] += M[i*n+j]*V[j]; 
		}
	}
}

template<class T>
double norm(int n, T V){
	double sum = 0;
	for(int i = 0; i < n; ++i){
		sum+=V[i]*V[i];
	}
	return sqrt(sum);
}

template<class T>
void norming(int n, T V){
	double nrm = norm(n, V);
	for(int i = 0; i < n; ++i){
		V[i] /= nrm;
	}
}

template<class T>
double mean_diff(int n, T V1, T V2){
	double sum = 0;
	for(int i = 0; i < n; ++i){
		sum += fabs(V1[i] - V2[i]);
	}
	return sum/n;
}

template<class T>
void power_method(int n, int m, T out, T M, double eps = 1e-16){
	double V[m];
	for(int i = 0; i < m; ++i){
		V[i] = ((double) std::rand())/RAND_MAX;
		
	}
	norming(m, V);
	double prevV[n];
	for(int i = 0; i < n; ++i){
		prevV[i] = 1e10;
	}
	double diff = 1e10;
	while(diff > eps){
		matrix_dot_vector(n, m, out, M, V);
		norming(m, out);
		
		std::copy(V, V + n, prevV);
		std::copy(out, out + n, V);
		
		diff = mean_diff(m, V, prevV);
		////std::cout << diff << std::endl;
	}	
}

int main(int argc, char **argv){
    
    MPI_Init(&argc, &argv);// init

    std::chrono::system_clock::time_point s;

    //get number of thread
    int id, ier, num_p;
	
    MPI_Comm_rank(MPI_COMM_WORLD, &id);//get number of thread
    MPI_Comm_size(MPI_COMM_WORLD, &num_p);

    double* M;//main matrix 
    double* V;//out
    double* prevV;//prev 
    double* out;
    int n, m;
    double diff = 1e10;
    double eps = 1e-7;
    int num_my_vectors = 0;//vectors, for 1 thread
    int num_vectors_per_process = 0;// for each of others

    int iter = 0;
    while(diff > eps){
        if (id == 0){
            s = std::chrono::system_clock::now();
	    if(iter == 0){		
		eps = 1e-7;

                //Input matrix size
                std::cout << "Input number of columns of matrix!" << std::endl;
                std::cin >> n;
                m = n;

	        num_my_vectors = n % (num_p-1);
	        num_vectors_per_process = n/(num_p-1);
	
                // allocate memory
		M = new double[m*n];
		V = new double[n];
		prevV = new double[n];
		out = new double[n];

		make_sth_matrix(M, n, m);
                
    	    	for(int i = 0; i < m; ++i){
    	    		V[i] = 1;
    	    	}

		norming(m, V);
		
    	    	for(int i = 0; i < n; ++i){
    	    		prevV[i] = 1e10;
    	    	}

    	        for(int i = 1; i < num_p; ++i){
    	            MPI_Send(&num_vectors_per_process, 1, MPI_INT, i, 0, MPI_COMM_WORLD);//количество строк матрицы, которые первый поток отправляет всем остальным
    	            MPI_Send(&n, 1, MPI_INT, i, 1, MPI_COMM_WORLD);//размерность строк матрицы
    	        }
	    }

    	    for(int i = 1; i < num_p; ++i){
    	        MPI_Send(M + (num_my_vectors + (i-1) * num_vectors_per_process) * n, n * num_vectors_per_process, MPI_DOUBLE_PRECISION, i, 2, MPI_COMM_WORLD);//Отправляет те строки матрицы, которые надо умножить на вектор для потока каждого
    	        MPI_Send(V, n, MPI_DOUBLE_PRECISION, i, 3, MPI_COMM_WORLD);
    	    }

            
	    for(int i = 0; i < num_my_vectors; ++i){
	    	out[i] = 0;
	    	for(int j = 0; j < m; ++j){
	    	    out[i] += M[i*n + j]*V[j]; 
	    	}
	    }

	    for(int i = 1; i < num_p; ++i){
	    	MPI_Recv(out + num_my_vectors + (i-1) * num_vectors_per_process, 
	    		num_vectors_per_process, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD,
	    				MPI_STATUS_IGNORE);///получает от всех потоков результат
	    }

    	    norming(n, out);
    	    
    	    std::copy(V, V + n, prevV);
    	    std::copy(out, out + n, V);
    	    
    	    diff = mean_diff(m, V, prevV);

            std::cout << "iter: " << iter << " " << diff << std::endl;

            //рассылаем diff каждому потоку, чтобы каждый поток знал когда ему кончать нахуй
	    for(int i = 1; i < num_p; ++i){
    	    	MPI_Send(&diff, 1, MPI_DOUBLE_PRECISION, i, 4, MPI_COMM_WORLD);
	    }

	    }else{//остальные
                if(iter == 0){
    	            MPI_Recv(&num_vectors_per_process, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    	            MPI_Recv(&n, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

    	        double* lM = new double[n*num_vectors_per_process];//указатель, куда поток получит матрицу
    	        double* lV = new double[n];//указатель куда поток получит вектор

    	        MPI_Recv(lM, n*num_vectors_per_process, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//получает матрицу
    	        MPI_Recv(lV, n, MPI_DOUBLE_PRECISION, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//получает вектор

    	        double *out = new double[num_vectors_per_process];
            	for(int i = 0; i < num_vectors_per_process; ++i){
	            out[i] = 0;
	            for(int j = 0; j < n; ++j){
	                out[i] += lM[i*n + j]*lV[j]; 
	            }
	        }
                delete [] out;

	        MPI_Send(out, num_vectors_per_process, MPI_DOUBLE_PRECISION, 0, id, MPI_COMM_WORLD);
    	        MPI_Recv(&diff, 1, MPI_DOUBLE_PRECISION, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	    ++iter;	
	}
	
    MPI_Barrier(MPI_COMM_WORLD);
	if (id == 0){
            std::chrono::system_clock::time_point e = std::chrono::system_clock::now();

            std::cout << "time: " << std::chrono::duration_cast<std::chrono::microseconds>(e - s).count() << std::endl;
		delete [] M;
		delete [] V;
		delete [] out;
		delete [] prevV;
	}
        MPI_Finalize(); 
}
