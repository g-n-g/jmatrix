package jmatrix;

/** Benchmark of equation solving by Cholesky LL decomposition. */
public final class SolveEqnLLBenchmark extends Benchmark
{
  private Matrix x, b;

  @Override
  public String name() {
    return "SolveEqn (LL)";
  }

  @Override
  public boolean isPd() {
    return true;
  }

  @Override
  protected void compute(Matrix A) {
    b = Matrix.ones(A.rows(), 1);
    Matrix L = A.choleskyL();
    Matrix y = L.backsL(b, false);
    x = L.T().backsU(y, false, y);
  }

  @Override
  protected double check(Matrix A) throws BenchmarkException {
    return checkMatrixEquals(b, A.mul(x));
  }

  public static void main(String[] args) {
    new SolveEqnLLBenchmark().run(args);
  }
}
