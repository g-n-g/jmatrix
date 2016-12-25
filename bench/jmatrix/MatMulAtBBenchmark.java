package jmatrix;

/** Matrix multiplication benchmark for computing A'*B. */
public final class MatMulAtBBenchmark extends MatMulBenchmark
{
  private Matrix AtB;

  @Override
  protected BenchmarkType type() {
    return BenchmarkType.AB_RG;
  }

  @Override
  protected Matrix result() {
    return AtB;
  }

  @Override
  protected void compute(Matrix A, Matrix B) throws BenchmarkException {
    AtB = A.T().mul(B);
  }
}
