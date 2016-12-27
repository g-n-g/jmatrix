package jmatrix;

/** Benchmark of computing the least squares solution
 *  for the full column rank case by reduced QR decomposition.
 */
public final class SolveLSFCREqnReducedQRBenchmark extends Benchmark
{
  private Matrix x;

  @Override
  protected BenchmarkType type() {
    return BenchmarkType.Ab_FCR;
  }

  @Override
  protected void compute(Matrix A, Matrix b) throws BenchmarkException {
    Matrix[] QR = A.reducedQR();
    Matrix y = QR[0].T().mul(b);
    x = QR[1].backsU(y);
  }

  @Override
  protected double check(Matrix A, Matrix b) throws BenchmarkException {
    return checkMatrixEquals(A.T().mul(b), A.T().mul(A).mul(x));
  }
}
