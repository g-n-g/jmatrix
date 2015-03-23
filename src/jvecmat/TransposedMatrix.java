package jvecmat;

public final class TransposedMatrix extends Matrix {

  TransposedMatrix(NonTransposedMatrix mat) {
    super(mat.array(), mat.cols(), mat.rows());
    this.mat = mat;
  }

  //----------------------------------------------------------------------------
    
  @Override
  public double get(int i, int j) {
    assert (0 <= i && i < rows());
    assert (0 <= j && j < cols());
    return array()[j][i];
  }
    
  @Override
  public void set(int i, int j, double value) {
    assert (0 <= i && i < rows());
    assert (0 <= j && j < cols());
    array()[j][i] = value;
  }

  //----------------------------------------------------------------------------

  @Override
  public Matrix T() {
    return mat;
  }
    
  //----------------------------------------------------------------------------
    
  private final NonTransposedMatrix mat;
}
