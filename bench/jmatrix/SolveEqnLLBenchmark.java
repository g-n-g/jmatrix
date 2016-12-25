package jmatrix;

/** Benchmark of equation solving by Cholesky LL decomposition. */
public final class SolveEqnLLBenchmark extends Benchmark
{
  private Matrix x;

  @Override
  protected BenchmarkType type() {
    return BenchmarkType.Ab_PD;
  }

  @Override
  protected void compute(Matrix A, Matrix b) throws BenchmarkException {
    Matrix L = A.choleskyL();
    Matrix y = L.backsL(b, false);
    x = L.T().backsU(y, false, y);
  }

  @Override
  protected double check(Matrix A, Matrix b) throws BenchmarkException {
    return checkMatrixEquals(b, A.mul(x));
  }
}
