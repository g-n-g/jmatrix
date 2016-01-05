package jmatrix;

/** Benchmark of inverting a positive definite matrix by the general method. */
public final class MatInvBenchmark extends Benchmark
{
  private Matrix Ainv;

  @Override
  public String name() {
    return "MatInv";
  }

  @Override
  protected BenchmarkType type() {
    return BenchmarkType.A_PD;
  }

  @Override
  protected void compute(Matrix A, Matrix bB) throws BenchmarkException {
    Ainv = A.inv();
  }

  @Override
  protected double check(Matrix A, Matrix bB) throws BenchmarkException {
    double delta = checkMatrixEye(A.mul(Ainv));
    delta = Math.max(delta, checkMatrixEye(Ainv.mul(A)));
    return delta;
  }
}
