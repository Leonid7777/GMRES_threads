#include <algorithm>


class Matvec
{
private:
    int size_of_matrix;
    int width_of_matrix;
    double* matrix;

public:
    Matvec(double* mtx, const int size, const int width);

    double* get_matrix() const;

    int get_size() const;

    int get_width() const;

    void mat_vec_foo(int down, int up, double* res, double* vec, int block_size, int size_of_matrix) const;

    void mat_vec(double* vec, double* res, std::thread* threads, int num_threads, int block_size, int* block_count_size_mas) const;

};