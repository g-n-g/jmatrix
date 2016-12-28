package jmatrix.bench;

import jmatrix.Matrix;

/** QR decomposition benchmark. */
public final class QRBenchmark extends Benchmark
{
  private Matrix Q, R;

  @Override
  protected BenchmarkType type() {
    return BenchmarkType.A_RG;
  }

  @Override
  protected void compute(Matrix A, Matrix bB) throws BenchmarkException {
    Matrix[] QR = A.QR();
    Q = QR[0];
    R = QR[1];
  }

  @Override
  protected double check(Matrix A, Matrix bB) throws BenchmarkException {
    double delta = checkMatrixEquals(A, Q.mul(R));
    delta = Math.max(delta, checkMatrixOrtho(Q));
    delta = Math.max(delta, checkMatrixUT(R));
    return delta;
  }
}
