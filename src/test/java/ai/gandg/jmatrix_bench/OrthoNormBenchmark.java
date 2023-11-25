package ai.gandg.jmatrix_bench;

import ai.gandg.jmatrix.Matrix;


/** Orthonormalization benchmark. */
public final class OrthoNormBenchmark extends Benchmark
{
  private Matrix M;

  @Override
  protected BenchmarkType type() {
    return BenchmarkType.A_RG;
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
