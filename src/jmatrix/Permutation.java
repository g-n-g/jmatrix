package jmatrix;

import java.util.Random;
import java.lang.IllegalArgumentException;

/**
 * Representation for permuting vectors and matrix rows/columns.
 */
public class Permutation {

  /**
   * Creates a new permutation encapsulating <code>data</code>.
   * The consistency of the permutation is verified, i.e. <code>data</code>
   * has to hold a permutation of <code>0, 1, ..., data.length-1</code>.
   *
   * @param data permutation data (not <code>null</code>)
   * @return new permutation encapsulating data
   * @throws IllegalArgumentException if <code>data</code> does not hold
   *         a permutation of <code>[0..data.length-1]</code>
   */
  public static Permutation create(int... data) {
    assert (data != null);
    final int n = data.length;
    int s = 0;
    for (int i = 0; i < n; ++i) {
      int v = data[i];
      if (v < 0) {
        throw new IllegalArgumentException("Negative element in permutation: " + v + "!");
      }
      if (v >= n) {
        throw new IllegalArgumentException("Out of range element in permutation: " + v + "!");
      }
      s += data[i];
    }
    if (2*s != (n-1)*n) {
      throw new IllegalArgumentException("Inconsistent permutation!");
    }
    return new Permutation(data);
  }

  /**
   * Creates a new identity permutation of the given <code>length</code>.
   *
   * @param length length of the permutation (non-negative)
   * @return identity permutation
   */
  public static Permutation eye(int length) {
    assert (length >= 0);
    int[] data = new int[length];
    return new Permutation(data).setToEye();
  }

  /**
   * Creates a random permutation by the Knuth shuffle algorithm.
   *
   * @param length length of the permutation (non-negative)
   * @param rng random number generator
   * @return random permutation
   */
  public static Permutation rand(int length, Random rng) {
    return eye(length).setToRand(rng);
  }

  //----------------------------------------------------------------------------

  /**
   * Returns a copy of the permutation placed into <code>result</code>.
   * If <code>result</code> is <code>null</code>, a new object is created.
   * The operation is skipped if <code>result</code> is equal to <code>this</code>.
   *
   * @param result appropriately sized storage for the copy
   *               or <code>null</code>
   * @return copy of the permutation
   */
  public Permutation copy(Permutation result) {
    if (result == null) { result = eye(length()); }
    assert(result.length() == length());
    if (result != this) {
      for (int i = 0; i < length(); ++i) { result.data[i] = data[i]; }
    }
    return result;
  }

  /**
   * Returns a new copy of the permutation.
   *
   * @return copy of the permutation
   */
  public Permutation copy() {
    return copy(null);
  }

  //----------------------------------------------------------------------------

  /**
   * Returns the element number of the permutation.
   *
   * @return length of permutation
   */
  public int length() {
    return data.length;
  }

  /**
   * Returns the permutation element at the specified <code>index</code>.
   *
   * @param index the index of the permutation element
   * @return the permutation element at <code>index</code>
   */
  public int get(int index) {
    assert (0 <= index && index < length());
    return data[index];
  }

  /**
   * Swaps two permutation elements at <code>index1</code> and <code>index2</code>.
   * This operation changes <code>this</code> permutation object,
   * use <code>copy().swap(...)</code> for a non-changing alternative.
   *
   * @param index1 the index of a permutation element
   * @param index2 the index of a permutation element
   * @return changed <code>this</code> object
   */
  public Permutation swap(int index1, int index2) {
    assert (0 <= index1 && index1 < length());
    assert (0 <= index2 && index2 < length());
    int i = data[index1];
    data[index1] = data[index2];
    data[index2] = i;
    return this;
  }

  //----------------------------------------------------------------------------

  /**
   * Sets <code>this</code> permutation to the identity permutation.
   *
   * @return <code>this</code> permutation set to identity
   */
  public Permutation setToEye() {
    for (int i = length()-1; i >= 0; --i) { data[i] = i; }
    return this;
  }

  /**
   * Sets <code>this</code> permutation to a random permutation
   * by the Knuth shuffle algorithm.
   *
   * @param rng random number generator
   * @return <code>this</code> permutation set randomly
   */
  public Permutation setToRand(Random rng) {
    for (int i = length()-1; i >= 0; --i) {
      int j = rng.nextInt(i+1);
      swap(i, j);
    }
    return this;
  }

  //----------------------------------------------------------------------------

  /**
   * Permutes a vector or the rows of a matrix <code>m</code>,
   * i.e. multiplying the matrix from the left with a permutation matrix
   * represented by <code>this</code> (in <code>result</code>).
   *
   * The length of <code>this</code> permutation has to be equal to the row
   * number of matrix <code>m</code>.
   * Furthermore, <code>result</code> has to have the same size
   * as <code>this</code> matrix.
   *
   * @param m matrix to permute (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>
   *                                      and not equal to <code>m</code>)
   * @return <code>result</code> having matrix <code>m</code> multiplied from
   *         the left with a permutation matrix represented by <code>this</code>
   * @see Matrix#mul(Permutation, Matrix)
   */
  public Matrix mul(Matrix m, Matrix result) {
    assert (m != null && result != null && m != result);
    Matrix r = result;
    if (m.rows() == 1 && m.cols() > 1) {
      m = m.T();
      r = result.T();
    }
    final int rows = m.rows(), cols = m.cols();
    assert (length() == rows);
    assert (r.rows() == rows && r.cols() == cols);
    for (int i = 0; i < rows; ++i) {
      int row = get(i);
      for (int j = 0; j < cols; ++j) {
        r.set(i, j, m.get(row,j));
      }
    }
    return result;
  }

  /**
   * Permutes the rows of matrix <code>m</code>, i.e. multiplying the matrix
   * from the left with a permutation matrix represented by <code>this</code>
   * (in new matrix).
   *
   * The length of <code>this</code>  permutation has to be equal to the row
   * number of matrix <code>m</code>.
   *
   * @param m matrix to permute (not <code>null</code>)
   * @return new matrix being equal to matrix <code>m</code> multiplied from
   *         the left with a permutation matrix represented by <code>this</code>
   * @see Matrix#mul(Permutation)
   */
  public Matrix mul(Matrix m) {
    return mul(m, Matrix.create(m.rows(), m.cols()));
  }

  //----------------------------------------------------------------------------

  /**
   * Returns the inverse of the permutation (in <code>result</code>).
   * The new permutation represents the transpose of the permutation matrix
   * corresponding to <code>this</code> permutation.
   *
   * The length of permutation <code>result</code> has to be equal to the length
   * of <code>this</code> permutation.
   *
   * @param result storage of the result (not <code>null</code>
   *                                      and not equal to <code>this</code>)
   * @return <code>result</code> holding the inverse permutation
   */
  public Permutation inv(Permutation result) {
    assert (result != null && result.length() == length());
    for (int i = length()-1; i >= 0; --i) { result.data[get(i)] = i; }
    return result;
  }

  /**
   * Returns the inverse of the permutation (in new permutation).
   * The new permutation represents the transpose of the permutation matrix
   * corresponding to <code>this</code> permutation.
   *
   * @return inverse permutation in new object
   */
  public Permutation inv() {
    return inv(new Permutation(new int[length()]));
  }

  //----------------------------------------------------------------------------

  @Override
  public String toString() {
    String str = "[";
    for (int i = 0; i < length(); ++i) {
      if (i != 0) { str += " "; }
      str += get(i);
    }
    return str + "]";
  }

  /**
   * Converts <code>this</code> permutation to a permutation matrix
   * having the ones at the (i,perm(i)) positions (in <code>result</code>).
   *
   * Parameter <code>result</code> has to be a square matrix with side length
   * equal to the length of this permutation.
   *
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> holding the permutation matrix
   */
  public Matrix toMatrix(Matrix result) {
    assert(result != null && result.rows() == length() && result.cols() == length());
    result.setToZeros();
    for (int i = length()-1; i >= 0; --i) { result.set(i, get(i), 1.0); }
    return result;
  }

  /**
   * Converts <code>this</code> permutation to a permutation matrix
   * having the ones at the (i,perm(i)) positions (in new matrix).
   *
   * @return the permutation matrix
   */
  public Matrix toMatrix() {
    return toMatrix(Matrix.create(length(), length()));
  }

  //----------------------------------------------------------------------------

  /**
   * Creates a permutation object which encapsulates <code>data</code>.
   *
   * @param data permutation data
   */
  Permutation(int[] data) {
    this.data = data;
  }

  /**
   * Returns the array representation of the permutation.
   * The data is not copied, so any changes to the returned array
   * will change (and might violate) the permutation.
   *
   * @return array representation of the permutation
   */
  final int[] array() {
    return data;
  }

  private final int[] data;
}
