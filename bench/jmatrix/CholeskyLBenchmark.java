package jmatrix;

/** LL Cholesky decomposition benchmark. */
public final class CholeskyLBenchmark extends Benchmark
{
  private Matrix L;

  @Override
  public String name() {
    return "Cholesky LL";
  }

  @Override
  public boolean isPd() {
    return true;
  }

  @Override
  protected void compute(Matrix A) {
    L = A.choleskyL();
  }

  @Override
  protected double check(Matrix A) throws BenchmarkException {
    checkMatrixSquareSize(A);
    checkMatrixSquareSize(L);
    double delta = checkMatrixEquals(A, L.mul(L.T()));
    delta = Math.max(delta, checkMatrixLT(L));
    return delta;
  }

  public static void main(String[] args) {
    new CholeskyLBenchmark().run(args);
  }
}
