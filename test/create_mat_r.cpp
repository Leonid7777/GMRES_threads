#include <chrono>
#include <cstdlib>
#include <cmath>

void
matrix_make(double* matrix, int size_of_matrix)
{
    std::srand(std::time(nullptr));
    for(int j = 0; j < size_of_matrix; j++) {
        for(int i = 0; i < size_of_matrix; i++) {
            matrix[i * size_of_matrix + j] = -1000.0 + 2000.0 * double(std::rand()) / RAND_MAX;
            // matrix[i * size_of_matrix + j] = i + j * size_of_matrix + 1;
        }
    }
}

void
right_part_make(double* right_part, int size_of_matrix, double& norm)
{
    std::srand(std::time(nullptr));
    for(int i = 0; i < size_of_matrix; i++) {
        right_part[i] = -1000.0 + 2000.0 * double(std::rand()) / RAND_MAX;
        // right_part[i] = i + 1;
        norm += right_part[i] * right_part[i];
    }
    norm = std::sqrt(norm);
}