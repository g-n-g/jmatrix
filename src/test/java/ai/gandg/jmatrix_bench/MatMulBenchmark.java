package ai.gandg.jmatrix_bench;

import ai.gandg.jmatrix.Matrix;


/** Superclass of matrix multiplication benchmarks. */
public abstract class MatMulBenchmark extends Benchmark
{
  protected abstract Matrix result();

  @Override
  protected double check(Matrix A, Matrix B) throws BenchmarkException {
    final int rows = A.rows(), cols = A.cols();
    // compute C = A'*bB
    Matrix C = Matrix.zeros(cols, cols);
    for (int i = 0; i < cols; ++i) {
      for (int j = 0; j < cols; ++j) {
        double s = 0.0;
        for (int k = 0; k < rows; ++k) {
          s += A.get(k,i) * B.get(k,j);
        }
        C.set(i, j, s);
      }
    }
    return checkMatrixEquals(C, result());
  }
}

