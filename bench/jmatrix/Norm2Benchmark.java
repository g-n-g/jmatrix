package jmatrix;

/** Spectral norm benchmark. */
public final class Norm2Benchmark extends Benchmark
{
  private double norm2;

  @Override
  public String name() {
    return "Norm2";
  }

  @Override
  protected void compute(Matrix A, Matrix bB) throws BenchmarkException {
    norm2 = A.norm2();
  }

  @Override
  protected double check(Matrix A, Matrix bB) throws BenchmarkException {
    final int rows = A.rows(), cols = A.cols();
    // check distance from Gershgorin discs
    double delta = Double.MAX_VALUE;
    for (int i = 0; i < Math.min(rows,cols); ++i) {
      double rRow = 0.0;
      for (int j = 0; j < cols; ++j) {
        if (i == j) { continue; }
        rRow += Math.abs(A.get(i,j));
      }
      double rCol = 0.0;
      for (int j = 0; j < rows; ++j) {
        if (i == j) { continue; }
        rCol += Math.abs(A.get(j,i));
      }
      double r = Math.min(rRow, rCol);
      double Aii = A.get(i,i);
      delta = Math.min(delta, distance(norm2,  Aii, r));
      delta = Math.min(delta, distance(norm2, -Aii, r));
    }
    return delta;
  }

  private double distance(double norm2, double c, double r) {
    return Math.max(0.0, Math.abs(norm2 - c) - r);
  }
}
