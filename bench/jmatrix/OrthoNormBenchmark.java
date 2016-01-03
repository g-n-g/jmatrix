package jmatrix;

/** Orthonormalization benchmark. */
public final class OrthoNormBenchmark extends Benchmark
{
  private Matrix M;

  @Override
  public String name() {
    return "OrthoNorm";
  }

  @Override
  protected void compute(Matrix A) {
    M = A.orthonormalize();
  }

  @Override
  protected double check(Matrix A) throws BenchmarkException {
    return checkMatrixOrthoOrZeroCols(M);
  }

  public static void main(String[] args) {
    new OrthoNormBenchmark().run(args);
  }
}
