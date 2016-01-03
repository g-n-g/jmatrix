package jmatrix;

/** Benchmark of inverting a positive definite matrix by a specialized method. */
public final class MatInvPsdBenchmark extends Benchmark
{
  private Matrix Ainv;

  @Override
  public String name() {
    return "MatInvPsd";
  }

  @Override
  public boolean isPd() {
    return true;
  }

  @Override
  protected void compute(Matrix A) {
    Ainv = A.invPsd();
  }

  @Override
  protected double check(Matrix A) throws BenchmarkException {
    double delta = checkMatrixEye(A.mul(Ainv));
    delta = Math.max(delta, checkMatrixEye(Ainv.mul(A)));
    return delta;
  }

  public static void main(String[] args) {
    new MatInvPsdBenchmark().run(args);
  }
}
