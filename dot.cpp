#include <fstream>
#include <vector>
#include <iostream>

template<class T>
void matrix_dot_vector(int n, int m, T out, T M, T V){

	for(int i = 0; i < n; ++i){
		out[i] = 0;
		for(int j = 0; j < m; ++j){
			out[i] += M[i*n+j]*V[j]; 
		}
	}
}

int main(){
	std::ifstream in_matrix;
	std::ifstream in_vector;

	in_matrix.open("matrix");
	in_vector.open("vector");

        int n, m;
	in_matrix >> n >> m;

	double M[m*n];
	double V[m];
	double out[n];

	for(int i = 0; i < n; ++i){
		for(int j = 0; j < m; ++j){
			in_matrix >> M[i*n+j];
		}
	}

	for(int i = 0; i < m; ++i){
		in_vector >> V[i];
	}

	matrix_dot_vector(n, m, out, M, V);

	for(int i = 0; i < n; ++i){
		std::cout << out[i] << " ";
	}


}
