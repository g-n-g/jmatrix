package jmatrix;

import static jmatrix.BasicBinaryOperation.MUL;

/** Reduced singular value decomposition benchmark. */
public final class ReducedSVDBenchmark extends Benchmark
{
  private Matrix U, S, V;

  @Override
  protected BenchmarkType type() {
    return BenchmarkType.A_RG;
  }

  @Override
  protected void compute(Matrix A, Matrix bB) throws BenchmarkException {
    Matrix[] USV = A.reducedSVD();
    U = USV[0];
    S = USV[1];
    V = USV[2];
  }

  @Override
  protected double check(Matrix A, Matrix bB) throws BenchmarkException {
    double delta = checkMatrixEquals(A, U.ewb(MUL, S.T()).mul(V.T()));
    delta = Math.max(delta, checkMatrixOrthoOrZeroCols(U));
    delta = Math.max(delta, checkMatrixOrthoOrZeroCols(V));
    return delta;
  }
}
