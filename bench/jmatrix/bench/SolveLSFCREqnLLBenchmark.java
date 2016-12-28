package jmatrix.bench;

import jmatrix.Matrix;

/** Benchmark of computing the least squares solution
 *  for the full column rank case by Cholesky decomposition.
 */
public final class SolveLSFCREqnLLBenchmark extends Benchmark
{
  private Matrix x;

  @Override
  protected BenchmarkType type() {
    return BenchmarkType.Ab_FCR;
  }

  @Override
  protected void compute(Matrix A, Matrix b) throws BenchmarkException {
    Matrix L = A.T().mul(A).choleskyL();
    Matrix y = L.backsL(A.T().mul(b));
    x = L.T().backsU(y);
  }

  @Override
  protected double check(Matrix A, Matrix b) throws BenchmarkException {
    return checkMatrixEquals(A.T().mul(b), A.T().mul(A).mul(x));
  }
}
