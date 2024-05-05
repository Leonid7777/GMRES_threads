#include <iostream>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <thread>
#include "create_mat_r.h"
#include "vec_methods.h"
#include "gmres.hpp"


int 
main()
{
    int num_threads = 5;
    int size_of_matrix = 4000;
    int block_size = 128;
    double norm = 0;

    double* matrix = new double[size_of_matrix * size_of_matrix];
    double* right_part = new double[size_of_matrix];
    double* res = new double[size_of_matrix];

    matrix_make(matrix, size_of_matrix);

    Matvec A(matrix, size_of_matrix, size_of_matrix);

    right_part_make(right_part, size_of_matrix, norm);

    auto start = std::chrono::high_resolution_clock::now();

    GMRES<Matvec>(A, right_part, res, num_threads);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Время работы программы: " << duration.count() << " миллисекунд" << std::endl;

    std::thread threads[num_threads];

    int count_block_size = (size_of_matrix + block_size - 1) / block_size;
    int* block_count_size_mas = new int[num_threads + 1];
    create_block_mas(block_count_size_mas, num_threads, count_block_size);

    double* h_i_g = new double[size_of_matrix];
    A.mat_vec(res, h_i_g, threads, num_threads, block_size, block_count_size_mas);

    double normis = 0;

    for(int j = 0; j < size_of_matrix; j++) {
        normis += (h_i_g[j] - right_part[j]) * (h_i_g[j] - right_part[j]);
    }

    std::cout << std::sqrt(normis) / norm << std::endl;

    delete[] right_part;
    delete[] res;
    delete[] matrix;
    delete[] h_i_g;
    delete[] block_count_size_mas;

    return 0;
}