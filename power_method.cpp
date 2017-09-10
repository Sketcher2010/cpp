#include <fstream>
#include <vector>
#include <iostream>
#include <random>
#include <chrono>

template<class T>
void print(T v, int n, int m){
    for(int i = 0; i < n; ++i){
	for(int j = 0; j < m; ++j){
		std::cout << v[i*n+j] << " ";
	}
        std::cout << std::endl;
    }
    std::cout << std::endl;
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
void power_method(int n, int m, T out, T M, double eps = 1e-7){
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
        int i = 0;
	while(diff > eps){
		matrix_dot_vector(n, m, out, M, V);
		norming(m, out);
		
		std::copy(V, V + n, prevV);
		std::copy(out, out + n, V);
		

		diff = mean_diff(m, V, prevV);
		///std::cout << i << " " << diff << std::endl;
                ++i;
	}	
}

int main(int argc, char **argv){
    //std::ifstream in_matrix;
    //std::ifstream in_vector;

    //in_matrix.open("t1");
    
    //int n, m;
    //std::cin >> n;
    //m = n;
    //in_matrix >> n >> m;

    //double M[m*n];
    //double test_copy_M[m*n];
    //double out[n];
    
    std::vector<int> ns = {100, 200, 400, 800, 1600, 3200, 6400, 12800}; 

    for(int i = 0; i < ns.size(); ++i){
        
        int n = ns[i], m = ns[i];

        double* M = new double[m*n];
        double* out = new double[n];

        make_sth_matrix(M, n, n);
        
        ///print(M, n, n);
        ///for(int i = 0; i < n; ++i){
        ///	for(int j = 0; j < m; ++j){
        ///		in_matrix >> M[i*n+j];
        ///	}
        ///}
        //std::copy(M, M + n*m, test_copy_M);

        norming(n*m, M);
        std::chrono::system_clock::time_point s = std::chrono::system_clock::now();
        power_method(n, m, out, M);
        std::chrono::system_clock::time_point e = std::chrono::system_clock::now();
        std::cout << "size: " << ns[i] << " time: " << std::chrono::duration_cast<std::chrono::microseconds>(e - s).count() << std::endl;
    }

}
