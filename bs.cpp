#include <ctime>

#include <mpi.h>

template<class T> 
void get_random_vector(T* v, int n){
    for(int i = 0; i < n; ++i)
        v[i] = (int)(10*((double) std::rand())/RAND_MAX);
}

template<class T>
void print(int n, T v, int id){
        std::cout << "id " << id << " start print" << std::endl;
	for(int i = 0; i < n; ++i){
		std::cout << v[i] << " ";
	}
        std::cout << "id " << id << " end print" << std::endl;
}

template<class T>
void swap_if(T& a, T& b, bool c){
    if(c)
        std::swap(a, b);
}

bool less(int a, int b){
    return a < b;
}

bool greater(int a, int b){
    return a > b;
}

template<class T, class C>
void swaping(T* V, int n, C comp){
    for(int j = 0; j + n/2 < n; ++j)
        swap_if(V[j], V[j+n/2], comp(V[j], V[j+n/2]));
}

int main(int argc, char **argv){
    
    std::srand(std::time(0));
    MPI_Init(&argc, &argv);// init

    int id, ier, num_p;
	
    MPI_Comm_rank(MPI_COMM_WORLD, &id);//get number of thread
    MPI_Comm_size(MPI_COMM_WORLD, &num_p);

    //Master thread only!
    if(id == 0){
        int n = 16;
        int* V = new int[n];
        get_random_vector(V, n);
        std::cout << "Initial sequence:" << std::endl;
        print(n, V, id);
        
        //send size
        for(int i = 1; i < num_p; ++i){
            MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        for(int i = 2; i < n; i*=2){
            for(int j = 1; j < n/i+1; ++j){
                MPI_Send(V + i*(j-1), i, MPI_INT, j, i, MPI_COMM_WORLD);   
            }
            for(int j = 1; j < n/i+1; ++j){
                MPI_Recv(V + i*(j-1), i, MPI_INT, j, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   
            }
            print(n, V, id);
        }
        for(int k = n; k > 1; k/=2){
            for(int ki = 0; ki < n/k; ++ki)
                    swaping(V + ki * k, k, greater);
        }
        print(n, V, id);
    }

    if(id > 0){
        int n = 0;
        MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int* V = new int[n/2];

        for(int i = 2; i < n; i*=2){
            if(id < n / i + 1){
                MPI_Recv(V, i, MPI_INT, 0, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                ///std::cout << "id " << id << "receive: " << std::endl;
                ///print(i, V, id);
                if(id % 2 == 1){
                    for(int k = i; k > 1; k/=2){
                        ///std::cout << "id " << id << "k " << k << std::endl;
                        for(int ki = 0; ki < i/k; ++ki)
                                swaping(V + ki * k, k, less);
                    }
                }else{
                    for(int k = i; k > 1; k/=2)
                        for(int ki = 0; ki < i/k; ++ki)
                                swaping(V + ki * k, k, greater);
                }
                ///std::cout << "id " << id << "after swaping: " << std::endl;
                ///print(i, V, id);
                MPI_Send(V, i, MPI_INT, 0, i, MPI_COMM_WORLD);
            }
        }
    }



    MPI_Finalize(); 
}

