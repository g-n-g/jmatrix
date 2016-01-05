package jmatrix;

/** Benchmark of QR decomposition when only R is computed. */
public final class QRnoQBenchmark extends Benchmark
{
  private Matrix Q, R;

  @Override
  public String name() {
    return "QR (noQ)";
  }

  @Override
  protected void compute(Matrix A, Matrix bB) throws BenchmarkException {
    R = Matrix.create(A.rows(), A.cols());
    final int t = Math.min(A.rows()-1, A.cols());
    A.QR(null, R, Matrix.create(t,1));
  }

  @Override
  protected double check(Matrix A, Matrix bB) throws BenchmarkException {
    double delta = checkMatrixUT(R);
    Matrix[] QR = A.QR();
    Matrix Q = QR[0];
    delta = Math.max(delta, checkMatrixEquals(A, Q.mul(R)));
    return delta;
  }
}
