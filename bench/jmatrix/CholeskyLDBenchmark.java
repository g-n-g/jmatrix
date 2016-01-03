package jmatrix;

/** LDL Cholesky decomposition benchmark. */
public final class CholeskyLDBenchmark extends Benchmark
{
  private Matrix L, D;

  @Override
  public String name() {
    return "Cholesky LDL";
  }

  @Override
  public boolean isPd() {
    return true;
  }

  @Override
  protected void compute(Matrix A) {
    Matrix[] LD = A.choleskyLD();
    L = LD[0];
    D = LD[1];
  }

  @Override
  protected double check(Matrix A) throws BenchmarkException {
    checkMatrixSquareSize(A);
    double delta = checkMatrixEquals(A, L.mul(D).mul(L.T()));
    delta = Math.max(delta, checkMatrixUnitLT(L));
    delta = Math.max(delta, checkMatrixDiag(D));
    return delta;
  }

  public static void main(String[] args) {
    new CholeskyLDBenchmark().run(args);
  }
}
