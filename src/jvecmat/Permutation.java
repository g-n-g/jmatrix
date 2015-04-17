package jvecmat;

/**
 * Representation for permuting vectors and matrix rows/columns.
 */
public class Permutation {

  /**
   * Creates a new identity permutation of the given <code>length</code>.
   */
  public static Permutation eye(int length) {
    data = new int[length];
    for (int i = 0; i < length; ++i) { data[i] = i; }
  }

  // TODO: rand

  //----------------------------------------------------------------------------

  /**
   * Creates a permutation object which encapsulates <code>data</code>.
   *
   * @param data permutation data
   */
  private Permutation(int[] data) {
    this.data = data;
  }

  /**
   * Returns the array representation of the permutation.
   * The data is not copied, so any changes to the returned array
   * will change (and might violate) the permutation.
   *
   * @return array representation of the permutation
   */
  public final int[] array() {
    return data;
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
    if (result == null) { result = new Permutation(length()); }
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
   *
   * @param index1 the index of a permutation element
   * @param index2 the index of a permutation element
   */
  public void swap(int index1, int index2) {
    assert (0 <= index1 && index1 < length());
    assert (0 <= index2 && index2 < length());
    int i = data[index1];
    data[index1] = data[index2];
    data[index2] = i;
  }

  //----------------------------------------------------------------------------

  // TODO: mul(Vector)
  // TODO: mul(Matrix)

  //----------------------------------------------------------------------------

  private final int[] data;
}
