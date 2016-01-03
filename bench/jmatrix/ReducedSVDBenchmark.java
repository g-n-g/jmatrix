package jmatrix;

import static jmatrix.BasicBinaryOperation.MUL;

/** Reduced singular value decomposition benchmark. */
public final class ReducedSVDBenchmark extends Benchmark
{
  private Matrix U, S, V;

  @Override
  public String name() {
    return "SVD (reduced)";
  }

  @Override
  protected void compute(Matrix A) {
    Matrix[] USV = A.reducedSVD();
    U = USV[0];
    S = USV[1];
    V = USV[2];
  }

  @Override
  protected double check(Matrix A) throws BenchmarkException {
    double delta = checkMatrixEquals(A, U.ewb(MUL, S.T()).mul(V.T()));
    delta = Math.max(delta, checkMatrixOrthoOrZeroCols(U));
    delta = Math.max(delta, checkMatrixOrthoOrZeroCols(V));
    return delta;
  }

  public static void main(String[] args) {
    new SvdReducedBenchmark().run(args);
  }
}
