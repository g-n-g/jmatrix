package jmatrix;

/** Reduced QR decomposition benchmark. */
public final class ReducedQRBenchmark extends Benchmark
{
  private Matrix Q, R;

  @Override
  public String name() {
    return "QR (reduced)";
  }

  @Override
  protected void compute(Matrix A) {
    Matrix[] QR = A.reducedQR();
    Q = QR[0];
    R = QR[1];
  }

  @Override
  protected double check(Matrix A) throws BenchmarkException {
    double delta = checkMatrixEquals(A, Q.mul(R));
    delta = Math.max(delta, checkMatrixOrthoCols(Q));
    delta = Math.max(delta, checkMatrixUT(R));
    return delta;
  }

  public static void main(String[] args) {
    new ReducedQRBenchmark().run(args);
  }
}
