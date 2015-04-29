package jvecmat;

final class DenseMatrix extends Matrix {

  DenseMatrix(double[][] data) {
    this.data = (data == null) ? new double[0][0] : data;
    trMat = new TDenseMatrix(this);
  }

  //---------------------------------------------------------------------------

  @Override
  public int rows() {
    return data.length;
  }

  @Override
  public int cols() {
    return (data.length != 0) ? data[0].length : 0;
  }

  @Override
  public double get(int i, int j) {
    assert (0 <= i && i < rows());
    assert (0 <= j && j < cols());
    return data[i][j];
  }

  @Override
  public void set(int i, int j, double value) {
    assert (0 <= i && i < rows());
    assert (0 <= j && j < cols());
    data[i][j] = value;
  }

  //---------------------------------------------------------------------------

  @Override
  public Matrix T() {
    return trMat;
  }

  //---------------------------------------------------------------------------

  private final double[][] data;
  private final TDenseMatrix trMat;
}
