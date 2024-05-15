#include <vector>
#include <algorithm>
#include <thread>


void 
eye_mat_place(std::vector<double>& mat, int n, int m, int place) 
{
    for(int i = place; i < m; i++) {
        mat[i * n + i] = 1;
    }
}

void 
both_position_mat_mul(std::vector<double>& res, std::vector<double>& vec_res, int n, int m, int position, double c, double s)
{
    double val_1;

    for(int j = 0; j < m; j++) {
        val_1 = c * res[j * n + position] + s * res[j * n + position + 1];
        res[j * n + position + 1] = - s * res[j * n + position] + c * res[j * n + position + 1];
        res[j * n + position] = val_1;

        val_1 = c * vec_res[j * n + position] + s * vec_res[j * n + position + 1];
        vec_res[j * n + position + 1] = - s * vec_res[j * n + position] + c * vec_res[j * n + position + 1];
        vec_res[j * n + position] = val_1;
    }

    val_1 = c * vec_res[m * n + position] + s * vec_res[m * n + position + 1];
    vec_res[m * n + position + 1] = - s * vec_res[m * n + position] + c * vec_res[m * n + position + 1];
    vec_res[m * n + position] = val_1;
}

void
place_mat_vec_foo(int down, int up, int m, double* res, double* mat, double* vec, int block_size, int size_of_matrix)
{
    for (int j = down; j < up; j += block_size) {
        int top = std::min(j + block_size, up);
        for (int i = 0; i < m; ++i) {
            for (int jj = j; jj < top; ++jj) {
                res[jj] += mat[i * size_of_matrix + jj] * vec[i];
            }
        }
    }
}

void 
place_mat_vec(std::vector<double> mat, double* vec, int m, int size_of_matrix, int block_size, std::thread* threads, int num_threads, int* block_mas)
{
    double* res = new double[m];

    for(int i = 0; i < m; i++) {
        res[i] = 0;
    }

    if(num_threads * block_size <= size_of_matrix) {
        for(int j = 0; j < num_threads; j++) {
            int top = std::min(block_size * block_mas[j + 1], m);

            threads[j] = std::thread(place_mat_vec_foo, block_size * block_mas[j], top, m, res, &mat[0], vec, block_size, size_of_matrix);
        }

        for(int i = 0; i < num_threads; i++) {
            threads[i].join();
        }
    } else {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; j++) {
                res[j] += mat[i * size_of_matrix + j] * vec[i];
            }
        }
    }

    for(int i = 0; i < m; i++) {
        vec[i] = res[i];
    }

    delete[] res;
}

void
vec_sub_mat_foo(int down, int up, int n, int m, double* res, double* mat, double* vec, int block_size)
{
    for (int j = down; j < up; j += block_size) {
        int top = std::min(j + block_size, up);
        for (int i = 0; i < n; ++i) {
            for (int jj = j; jj < top; ++jj) {
                res[jj] -= mat[i * m + jj] * vec[i];
            }
        }
    }
}

void
vec_sub_mat(double* res, double* mat, double* vec, int m, int n, int block_size, std::thread* threads, int num_threads, int* block_mas)
{   
    n += 1;

    if(num_threads * block_size <= m) {
        for (int j = 0; j < num_threads; j++) {
            int top = std::min(block_size * block_mas[j + 1], m);
            
            threads[j] = std::thread(vec_sub_mat_foo, block_size * block_mas[j], top, n, m, res, mat, vec, block_size);
        }

        for(int i = 0; i < num_threads; i++) {
            threads[i].join();
        }
    } else {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; j++) {
                res[j] -= mat[i * m + j] * vec[i];
            }
        }
    }
}