package jmatrix;

/** Benchmark of inverting a positive definite matrix by a specialized method. */
public final class MatInvPsdBenchmark extends Benchmark
{
  private Matrix Ainv;

  @Override
  protected BenchmarkType type() {
    return BenchmarkType.A_PD;
  }

  @Override
  protected void compute(Matrix A, Matrix bB) throws BenchmarkException {
    Ainv = A.invPsd();
  }

  @Override
  protected double check(Matrix A, Matrix bB) throws BenchmarkException {
    double delta = checkMatrixEye(A.mul(Ainv));
    delta = Math.max(delta, checkMatrixEye(Ainv.mul(A)));
    return delta;
  }
}
