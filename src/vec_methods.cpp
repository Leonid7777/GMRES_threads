#include <cmath>
#include <thread>
#include <algorithm>


void
scal_prod_foo(int down, int top, double* res, double* vec, double* mat, int m)
{
    for(int i = down; i < top; i++) {
        for(int j = 0; j < m; j++) {
            res[i] += vec[j] * mat[i * m + j];
        }
    }
}

void
scal_prod(double* res, double* vec, double* mat, int m, int n, std::thread* threads, int num_threads)
{
    n += 1;

    if(n >= 100) {
        int count = (n + num_threads - 1) / num_threads;

        for(int i = 0; i < num_threads; i++) {
            int top = std::min((i + 1) * count, n);

            threads[i] = std::thread(scal_prod_foo, i * count, top, res, vec, mat, m);
        }

        for(int i = 0; i < num_threads; i++) {
            threads[i].join();
        }
    } else {
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < m; j++) {
                res[i] += vec[j] * mat[i * m + j];
            }
        }
    }
}

double
norm_vec(const double* vec, int n)
{
    double val = 0;
    for(int i = 0; i < n; i++) {
        val += vec[i] * vec[i]; 
    }

    return std::sqrt(val);
}

void
scal_vec(double* vec, int n, double val)
{
    for(int i = 0; i < n; i++) {
        vec[i] /= val;
    }
}

double
rotatation(double a, double b)
{
    return std::sqrt(a * a + b * b);
}

void
create_block_mas(int* block_count_size_mas, int num_threads, int count_block_size)
{
    for(int i = 0; i < num_threads + 1; i++) {
        block_count_size_mas[i] = 0;
    }
    int bl_c = 0;
    while(bl_c < count_block_size) {
        block_count_size_mas[bl_c % num_threads + 1]++;
        bl_c++;
    }
    for(int i = 2; i < num_threads + 1; i++) {
        block_count_size_mas[i] += block_count_size_mas[i - 1];
    }
}