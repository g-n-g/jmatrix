package jmatrix;

/** Matrix multiplication benchmark for computing A'*A. */
public final class MatMulAtABenchmark extends MatMulBenchmark
{
  private Matrix AtA;

  @Override
  protected BenchmarkType type() {
    return BenchmarkType.A_RG;
  }

  @Override
  protected void compute(Matrix A, Matrix B) throws BenchmarkException {
    AtA = A.T().mul(A);
  }

  @Override
  protected Matrix result() {
    return AtA;
  }

  @Override
  protected double check(Matrix A, Matrix B) throws BenchmarkException {
    double delta = super.check(A, A);
    delta = Math.max(delta, checkMatrixSymmetric(AtA));
    return delta;
  }
}
