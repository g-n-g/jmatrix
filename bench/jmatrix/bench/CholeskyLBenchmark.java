package jmatrix.bench;

import jmatrix.Matrix;

/** LL Cholesky decomposition benchmark. */
public final class CholeskyLBenchmark extends Benchmark
{
  private Matrix L;

  @Override
  protected BenchmarkType type() {
    return BenchmarkType.A_PD;
  }

  @Override
  protected void compute(Matrix A, Matrix bB) throws BenchmarkException {
    L = A.choleskyL();
  }

  @Override
  protected double check(Matrix A, Matrix bB) throws BenchmarkException {
    checkMatrixSquareSize(A);
    checkMatrixSquareSize(L);
    double delta = checkMatrixEquals(A, L.mul(L.T()));
    delta = Math.max(delta, checkMatrixLT(L));
    return delta;
  }
}
