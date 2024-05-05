#include <algorithm>
#include <thread>
#include "world_class.h"


void
Matvec::mat_vec_foo(int down, int up, double* res, double* vec, int block_size, int size_of_matrix) const
{
    for (int j = down; j < up; j += block_size) {
        int top = std::min(j + block_size, up);
        for(int i = 0; i < size_of_matrix; i++) {
            for (int jj = j; jj < top; ++jj) {
                res[jj] += matrix[i * size_of_matrix + jj] * vec[i];
            }
        }
    }
}

void 
Matvec::mat_vec(double* vec, double* res, std::thread* threads, int num_threads, int block_size, int* block_mas) const
{
    if(num_threads * block_size <= size_of_matrix) {
        for(int j = 0; j < num_threads; j++) {
            int top = std::min(block_size * block_mas[j + 1], width_of_matrix);

            threads[j] = std::thread([=](){ this->mat_vec_foo(block_size * block_mas[j], top, res, vec, block_size, size_of_matrix); });
        }

        for(int i = 0; i < num_threads; i++) {
            threads[i].join();
        }
    } else {
        for (int i = 0; i < size_of_matrix; ++i) {
            for (int j = 0; j < width_of_matrix; j++) {
                res[j] += matrix[i * size_of_matrix + j] * vec[i];
            }
        }
    }

}

double*
Matvec::get_matrix() const
{
    return matrix;
}

int
Matvec::get_size() const
{
    return size_of_matrix;
}

int 
Matvec::get_width() const
{
    return width_of_matrix;
}

Matvec::Matvec(double* mtx, const int size, const int width)
{
    matrix = mtx;
    size_of_matrix = size;
    width_of_matrix = width;
}