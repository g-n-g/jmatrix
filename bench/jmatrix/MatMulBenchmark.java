package jmatrix;

/** Matrix multiplication benchmark. */
public final class MatMulBenchmark extends Benchmark
{
  private Matrix AtA;

  @Override
  public String name() {
    return "MatMul";
  }

  @Override
  protected void compute(Matrix A) {
    AtA = A.T().mul(A);
  }

  @Override
  protected double check(Matrix A) throws BenchmarkException {
    return checkMatrixSymmetric(AtA);
  }

  public static void main(String[] args) {
    new MatMulBenchmark().run(args);
  }
}
