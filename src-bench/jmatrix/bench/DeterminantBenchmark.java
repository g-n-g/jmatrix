package jmatrix.bench;

import jmatrix.Matrix;

/** Determinant benchmark. */
public final class DeterminantBenchmark extends Benchmark
{
  private double det;

  @Override
  protected BenchmarkType type() {
    return BenchmarkType.A_PD;
  }

  @Override
  protected double customScaling(Matrix A, Matrix bB) {
    double cs = 0.0;
    for (int i = 0; i < Math.min(A.rows(), A.cols()); ++i) {
      cs += Math.abs(A.get(i,i));
    }
    return 1.0 / cs;
  }

  @Override
  protected void compute(Matrix A, Matrix bB) throws BenchmarkException {
    det = A.det();
  }

  @Override
  protected double check(Matrix A, Matrix bB) throws BenchmarkException {
    Matrix L = A.choleskyL();
    double detchk = L.prodDiag();
    detchk = detchk*detchk;
    double delta = Math.abs(det - detchk);
    checkDelta(delta);
    return delta;
  }
}
