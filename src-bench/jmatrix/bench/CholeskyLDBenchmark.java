package jmatrix.bench;

import jmatrix.Matrix;

/** LDL Cholesky decomposition benchmark. */
public final class CholeskyLDBenchmark extends Benchmark
{
  private Matrix L, D;

  @Override
  protected BenchmarkType type() {
    return BenchmarkType.A_PD;
  }

  @Override
  protected void compute(Matrix A, Matrix bB) throws BenchmarkException {
    Matrix[] LD = A.choleskyLD();
    L = LD[0];
    D = LD[1];
  }

  @Override
  protected double check(Matrix A, Matrix bB) throws BenchmarkException {
    checkMatrixSquareSize(A);
    double delta = checkMatrixEquals(A, L.mul(D).mul(L.T()));
    delta = Math.max(delta, checkMatrixUnitLT(L));
    delta = Math.max(delta, checkMatrixDiag(D));
    return delta;
  }
}
