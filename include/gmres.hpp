#include <vector>
#include <cmath>
#include <thread>
#include "vec_methods.h"
#include "mat_methods.h"
#include "../class/world_class.h"


template <class T>
void
GMRES(const T& A, const double* right_part, double* res, int num_threads)
{
    static constexpr int block_size = 128;

    int size_of_matrix = A.get_size();
    double norm = norm_vec(right_part, size_of_matrix);
    int krylov_count = 0;

    int kr_v = 2;
    std::vector<double> Q_vec ((size_of_matrix + 1) * kr_v);
    std::vector<double> H_vec ((size_of_matrix + 1) * (kr_v - 1));
    std::vector<double> krylov_subspaces_vec (size_of_matrix * kr_v);

    std::thread threads[num_threads];

    int count_block_size = (size_of_matrix + block_size - 1) / block_size;
    int* block_count_size_mas = new int[num_threads + 1];
    create_block_mas(block_count_size_mas, num_threads, count_block_size);


    int count_block_one_size = (size_of_matrix + block_size) / block_size;
    int* block_count_one_size_mas = new int[num_threads + 1];
    create_block_mas(block_count_one_size_mas, num_threads, count_block_one_size);

    eye_mat_place(Q_vec, size_of_matrix + 1, kr_v, 0);

    for(int i = 0; i < size_of_matrix; i++) {
        res[i] = 0;
        krylov_subspaces_vec[i] = right_part[i] / norm;
    }

    double denom, c, s;
    double err = 10;

    while(err >= 0.001 && krylov_count < size_of_matrix) {

        if(krylov_count + 1 == kr_v) {
            int place = kr_v;
            kr_v = std::min(2 * kr_v, size_of_matrix + 1);
            Q_vec.resize(kr_v * (size_of_matrix + 1));
            H_vec.resize((kr_v - 1) * (size_of_matrix + 1));
            krylov_subspaces_vec.resize(kr_v * size_of_matrix);
            eye_mat_place(Q_vec, size_of_matrix + 1, kr_v, place);
        }

        A.mat_vec(&krylov_subspaces_vec[krylov_count * size_of_matrix], &krylov_subspaces_vec[(krylov_count + 1) * size_of_matrix], threads, num_threads, block_size, block_count_size_mas);

        scal_prod(&H_vec[krylov_count * (size_of_matrix + 1)], &krylov_subspaces_vec[(krylov_count + 1) * size_of_matrix], &krylov_subspaces_vec[0], size_of_matrix, krylov_count, threads, num_threads);

        vec_sub_mat(&krylov_subspaces_vec[(krylov_count + 1) * size_of_matrix], &krylov_subspaces_vec[0], &H_vec[krylov_count * (size_of_matrix + 1)], size_of_matrix, krylov_count, block_size, threads, num_threads, block_count_size_mas);

        H_vec[krylov_count * (size_of_matrix + 1) + krylov_count + 1] = norm_vec(&krylov_subspaces_vec[(krylov_count + 1) * size_of_matrix], size_of_matrix);

        if(H_vec[krylov_count * (size_of_matrix + 1) + krylov_count + 1] == 0) {
            krylov_count++;
            break;
        } else {
            scal_vec(&krylov_subspaces_vec[(krylov_count + 1) * size_of_matrix], size_of_matrix, H_vec[krylov_count * (size_of_matrix + 1) + krylov_count + 1]);
        }

        place_mat_vec(Q_vec, &H_vec[krylov_count * (size_of_matrix + 1)], krylov_count + 1, size_of_matrix + 1, block_size, threads, num_threads, block_count_one_size_mas);

        denom = rotatation(H_vec[krylov_count * (size_of_matrix + 1) + krylov_count], H_vec[krylov_count * (size_of_matrix + 1) + krylov_count + 1]);

        c =  H_vec[krylov_count * (size_of_matrix + 1) + krylov_count] / denom;
        s =  H_vec[krylov_count * (size_of_matrix + 1) + krylov_count + 1] / denom;

        both_position_mat_mul(H_vec, Q_vec, size_of_matrix + 1, krylov_count + 1, krylov_count, c, s);

        err = std::abs(Q_vec[krylov_count + 1]);
        krylov_count++;
    }

    double* y = new double[size_of_matrix];

    for(int i = krylov_count - 1; i >= 0; i--) {
        y[i] = Q_vec[i] * norm;
        for(int j = i + 1; j < krylov_count; j++) {
            y[i] -= H_vec[j * (size_of_matrix + 1) + i] * y[j];
        }
        y[i] /=  H_vec[i * (size_of_matrix + 1) + i];
    }

    for(int i = 0; i < krylov_count; i++) {
        for(int j = 0; j < size_of_matrix; j++) {
            res[j] += y[i] * krylov_subspaces_vec[i * size_of_matrix + j];
        }
    }

    delete[] y;
    delete[] block_count_size_mas;
    delete[] block_count_one_size_mas;
}