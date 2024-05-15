#include <vector>


void vec_sub_vec(double* f_vec, double* s_vec, int n, double val);
double norm_vec(const double* vec, int n);
void scal_vec(double* vec, int n, double val);
double rotatation(double a, double b);
void scal_prod(double* res, double* vec, double* mat, int m, int n, std::thread* threads, int num_threads); 
void create_block_mas(int* block_count_size_mas, int num_threads, int count_block_size);