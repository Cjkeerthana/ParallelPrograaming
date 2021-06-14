extern "C"{
  double seconds( void );
  void print_matrix(double *M, size_t rows, size_t cols);
  void cuda_initialize(double* A, size_t m, size_t n, size_t p, int rank, double* t_cpy);
  void cuda_vector_mul(double* B, double** C, size_t m, size_t n, size_t p, double* t_cpy, double* t_calc);
  void cuda_stop();
}
