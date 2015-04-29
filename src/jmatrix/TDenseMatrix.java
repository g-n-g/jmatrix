package jmatrix;

final class TDenseMatrix extends Matrix {

  TDenseMatrix(DenseMatrix mat) {
    this.mat = mat;
  }

  //----------------------------------------------------------------------------

  @Override
  public int rows() {
    return mat.cols();
  }

  @Override
  public int cols() {
    return mat.rows();
  }

  @Override
  public double get(int i, int j) {
    assert (0 <= i && i < rows());
    assert (0 <= j && j < cols());
    return mat.get(j, i);
  }

  @Override
  public void set(int i, int j, double value) {
    assert (0 <= i && i < rows());
    assert (0 <= j && j < cols());
    mat.set(j, i, value);
  }

  //----------------------------------------------------------------------------

  @Override
  public Matrix T() {
    return mat;
  }

  //----------------------------------------------------------------------------

  private final DenseMatrix mat;
}
