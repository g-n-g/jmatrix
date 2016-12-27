package jmatrix;

/** Benchmark of computing the least squares solution
 *  for the full column rank case by QR decomposition without computing Q.
 */
public final class SolveLSFCREqnQRnoQBenchmark extends Benchmark
{
  private Matrix x;

  @Override
  protected BenchmarkType type() {
    return BenchmarkType.Ab_FCR;
  }

  @Override
  protected void compute(Matrix A, Matrix b) throws BenchmarkException {
    Matrix[] QR = A.QR(false); // without computing Q
    Matrix y = QR[1].T().backsL(A.T().mul(b));
    x = QR[1].backsU(y);
  }

  @Override
  protected double check(Matrix A, Matrix b) throws BenchmarkException {
    return checkMatrixEquals(A.T().mul(b), A.T().mul(A).mul(x));
  }
}
