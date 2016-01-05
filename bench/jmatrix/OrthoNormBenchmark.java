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
  protected void compute(Matrix A, Matrix bB) throws BenchmarkException {
    M = A.orthonormalize();
  }

  @Override
  protected double check(Matrix A, Matrix bB) throws BenchmarkException {
    return checkMatrixOrthoOrZeroCols(M);
  }
}
