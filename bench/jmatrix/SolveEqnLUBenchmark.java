package jmatrix;

/** Benchmark of equation solving by LU decomposition. */
public final class SolveEqnLUBenchmark extends Benchmark
{
  private Matrix x, b;

  @Override
  public String name() {
    return "SolveEqn (LU)";
  }

  @Override
  public boolean isPd() {
    return true;
  }

  @Override
  protected void compute(Matrix A) {
    b = Matrix.ones(A.rows(), 1);
    Matrix[] LU = A.LU();
    Matrix y = LU[0].backsL(b, true);
    x = LU[1].backsU(y, false, y);
  }

  @Override
  protected double check(Matrix A) throws BenchmarkException {
    return checkMatrixEquals(b, A.mul(x));
  }

  public static void main(String[] args) {
    new SolveEqnLUBenchmark().run(args);
  }
}
