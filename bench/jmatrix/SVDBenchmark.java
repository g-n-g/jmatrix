package jmatrix;

/** Singular value decomposition benchmark. */
public final class SVDBenchmark extends Benchmark
{
  private Matrix U, S, V;

  @Override
  public String name() {
    return "SVD";
  }

  @Override
  protected void compute(Matrix A) {
    Matrix[] USV = A.SVD();
    U = USV[0];
    S = USV[1];
    V = USV[2];
  }

  @Override
  protected double check(Matrix A) throws BenchmarkException {
    final int rows = A.rows(), cols = A.cols();

    S = Matrix.diag(S);
    if (rows > cols) {
      S = Matrix.vertcat(S, Matrix.zeros(rows-cols,cols));
    }
    else if (rows < cols) {
      S = Matrix.horzcat(S, Matrix.zeros(rows,cols-rows));
    }

    double delta = checkMatrixEquals(A, U.mul(S).mul(V.T()));
    delta = Math.max(delta, checkMatrixOrtho(U));
    delta = Math.max(delta, checkMatrixOrtho(V));
    return delta;
  }

  public static void main(String[] args) {
    new SvdBenchmark().run(args);
  }
}
