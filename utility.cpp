#include <fstream>
#include <vector>
#include <iostream>
#include <random>

void print(int n, double* v){
	for(int i = 0; i < n; ++i){
		std::cout << v[i] << " ";
	}
}

void print(int n, int m, double* matrix){
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < m; ++j){
			std::cout << matrix[i*m + j] << " ";
		}
		std::cout << std::endl;
	}
		std::cout << std::endl;
}

void make_sth_matrix(double* matrix, int n, int m){
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < m; ++j){
			int idx = (int) (n*((double) std::rand())/RAND_MAX);
			matrix[idx*n + i] += 1.0/n;
		}	
	}
}

void matrix_dot_vector(int n, int m, double* out, double* matrix, double* vector){
	for(int i = 0; i < n; ++i){
		out[i] = 0;
		for(int j = 0; j < m; ++j){
			out[i] += matrix[i*n+j]*vector[j]; 
		}
	}
}

double norm(int n, double* vector){
	double sum = 0;
	for(int i = 0; i < n; ++i){
		sum+=vector[i]*vector[i];
	}
	return sqrt(sum);
}

void norming(int n, double* vector){
	double nrm = norm(n, vector);
	for(int i = 0; i < n; ++i){
		vector[i] /= nrm;
	}
}

double mean_diff(int n, double* vector1, double* vector2){
	double sum = 0;
	for(int i = 0; i < n; ++i){
		sum += fabs(vector1[i] - vector2[i]);
	}
	return sum/n;
}

void power_method(int n, int m, double* out, double* matrix, double eps){
	double vector[m];
	for(int i = 0; i < m; ++i){
		vector[i] = ((double) std::rand())/RAND_MAX;
		
	}
	norming(m, vector);
	double prevvector[n];
	for(int i = 0; i < n; ++i){
		prevvector[i] = 1e10;
	}
	double diff = 1e10;
	while(diff > eps){
		matrix_dot_vector(n, m, out, matrix, vector);
		norming(m, out);
		
		std::copy(vector, vector + n, prevvector);
		std::copy(out, out + n, vector);
		
		diff = mean_diff(m, vector, prevvector);
		std::cout << diff << std::endl;
	}	
}
