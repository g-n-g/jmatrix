package ai.gandg.jmatrix_bench;

import ai.gandg.jmatrix.Matrix;


/** Reduced QR decomposition benchmark. */
public final class ReducedQRBenchmark extends Benchmark
{
  private Matrix Q, R;

  @Override
  protected BenchmarkType type() {
    return BenchmarkType.A_RG;
  }

  @Override
  protected void compute(Matrix A, Matrix bB) throws BenchmarkException {
    Matrix[] QR = A.reducedQR();
    Q = QR[0];
    R = QR[1];
  }

  @Override
  protected double check(Matrix A, Matrix bB) throws BenchmarkException {
    double delta = checkMatrixEquals(A, Q.mul(R));
    delta = Math.max(delta, checkMatrixOrthoCols(Q));
    delta = Math.max(delta, checkMatrixUT(R));
    return delta;
  }
}
