void print(int n, double* v);
void print(int n, int m, double* matrix);
void make_sth_matrix(double* matrix, int n, int m);
void matrix_dot_vector(int n, int m, double* out, double* matrix, double* vector);
double norm(int n, double* vector);
void norming(int n, double* vector);
double mean_diff(int n, double* vector1, double* vector2);
void power_method(int n, int m, double* out, double* matrix, double eps = 1e-16);
