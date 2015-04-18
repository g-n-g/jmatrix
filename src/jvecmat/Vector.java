package jvecmat;

import java.util.Random;

/**
 * Representation of column vectors.
 */
public class Vector implements VecMat {
  /**
   * Shared object for all empty vectors.
   */
  public static final Vector EMPTY = new Vector(new double[0]);

  /**
   * Creates a vector using the provided data.
   * An empty vector is created if <code>data</code> is <code>null</code>
   * or has zero length.
   *
   * @param data elements of the vector
   * @return vector which encapsulates <code>data</code> (not copied)
   */
  public static Vector create(double[] data) {
    if (data == null || data.length == 0) { return EMPTY; }
    return new Vector(data);
  }

  /**
   * Creates a vector of given length having uninitialized elements.
   * An empty vector is created if <code>length</code> is not positive.
   *
   * @param length number of elements
   * @return vector of size <code>length</code> having uninitialized elements
   */
  public static Vector create(int length) {
    if (0 >= length) return EMPTY;
    return create(new double[length]);
  }

  /**
   * Creates a vector of the given <code>length</code>
   * and initializes its elements to <code>value</code>.
   * An empty vector is created if <code>length</code> is not positive.
   *
   * @param length number of elements
   * @param value value of elements
   * @return constant (column) vector of size <code>length</code>
   *         having its elements set to <code>value</code>
   */
  public static Vector constant(int length, double value) {
    Vector v = create(length);
    v.setToConstant(value);
    return v;
  }

  /**
   * Creates a vector of the given <code>length</code>
   * and initializes its elements to zero.
   * An empty vector is created if <code>length</code> is not positive.
   *
   * @param length number of elements
   * @return zero (column) vector of size <code>length</code>
   */
  public static Vector zero(int length) {
    return constant(length, 0.0);
  }

  /**
   * Creates a vector of the given <code>length</code>
   * and initializes its elements to one.
   * An empty vector is created if <code>length</code> is not positive.
   *
   * @param length number of elements
   * @return all one (column) vector of size <code>length</code>
   */
  public static Vector one(int length) {
    return constant(length, 1.0);
  }

  /**
   * Creates a standard unit vector of the given <code>length</code>
   * having value one on position <code>onePosition</code> and zeros elsewhere.
   * An empty vector is created if <code>length</code> is not positive.
   *
   * @param length number of elements
   * @param onePosition position of the one value
   * @return standard unit (column) vector of size <code>length</code>
   */
  public static Vector unit(int length, int onePosition) {
    Vector v = zero(length);
    v.data[onePosition] = 1.0;
    return v;
  }

  /**
   * Creates a random vector of the given <code>length</code>
   * drawing its elements from the uniform distribution on [0,1].
   * An empty vector is created if <code>length</code> is not positive.
   *
   * @param length number of elements
   * @param rng random number generator
   * @return uniform random vector of size <code>length</code>
   */
  public static Vector rand(int length, Random rng) {
    Vector v = create(length);
    v.setToRand(rng);
    return v;
  }

  /**
   * Creates a random vector of the given <code>length</code>
   * drawing its elements from the standard normal distribution.
   * An empty vector is created if <code>length</code> is not positive.
   *
   * @param length number of elements
   * @param rng random number generator
   * @return standard normal random vector of size <code>length</code>
   */
  public static Vector randN(int length, Random rng) {
    Vector v = create(length);
    v.setToRandN(rng);
    return v;
  }

  //----------------------------------------------------------------------------

  /**
   * Creates a (column) vector object
   * which encapsulates <code>data</code>.
   *
   * @param data vector data (not copied)
   */
  protected Vector(double[] data) {
    this.data = data;
  }

  /**
   * Returns the array representation of the vector. The data is not copied,
   * so any changes to the returned array will change the vector too.
   *
   * @return array representation of the vector
   */
  public final double[] array() {
    return data;
  }

  //----------------------------------------------------------------------------

  @Override
  public final boolean hasNaN() {
    for (int i = 0; i < length(); ++i) {
      if (Double.isNaN(get(i))) { return true; }
    }
    return false;
  }

  @Override
  public final boolean hasInf() {
    for (int i = 0; i < length(); ++i) {
      if (Double.isInfinite(get(i))) { return true; }
    }
    return false;
  }

  @Override
  public final void replaceNaNandInf(double nan, double negInf, double posInf) {
    double e;
    for (int i = 0; i < length(); ++i) {
      e = get(i);
      if (Double.isNaN(e)) { set(i, nan); }
      else if (Double.POSITIVE_INFINITY == e) { set(i, posInf); }
      else if (Double.NEGATIVE_INFINITY == e) { set(i, negInf); }
    }
  }

  //----------------------------------------------------------------------------

  /**
   * Returns a copy of the vector placed into <code>result</code>.
   * If <code>result</code> is <code>null</code>, a new object is created.
   * The operation is skipped if <code>result</code> is equal to <code>this</code>.
   *
   * @param result appropriately sized storage for the copy
   *               or <code>null</code>
   * @return copy of the vector
   */
  public Vector copy(Vector result) {
    if (result == null) { result = create(length()); }
    assert(result.length() == length());
    if (result != this) {
      for (int i = 0; i < length(); ++i) { result.set(i, get(i)); }
    }
    return result;
  }

  @Override
  public Vector copy() {
    return copy(null);
  }

  //----------------------------------------------------------------------------

  /**
   * Returns the number of elements of the vector.
   *
   * @return length of the vector
   */
  public final int length() {
    return data.length;
  }

  /**
   * Returns the vector element at the specified <code>index</code>.
   *
   * @param index the index of the vector element
   * @return the vector element at <code>index</code>
   */
  public final double get(int index) {
    assert (0 <= index && index < length());
    return data[index];
  }

  /**
   * Set the vector element at the specified <code>index</code>
   * to <code>value</code>.
   *
   * @param index the index of the vector element
   * @param value the new value
   */
  public final void set(int index, double value) {
    assert (0 <= index && index < length());
    data[index] = value;
  }

  //----------------------------------------------------------------------------

  @Override
  public Vector setToConstant(double c) {
    for (int i = 0; i < length(); ++i) { set(i, c); }
    return this;
  }

  @Override
  public Vector setToZero() {
    return setToConstant(0.0);
  }

  @Override
  public Vector setToOne() {
    return setToConstant(1.0);
  }

  @Override
  public Vector setToRand(Random rng) {
    for (int i = 0; i < length(); ++i) { set(i, rng.nextDouble()); }
    return this;
  }

  @Override
  public Vector setToRandN(Random rng) {
    for (int i = 0; i < length(); ++i) { set(i, rng.nextGaussian()); }
    return this;
  }

  //----------------------------------------------------------------------------

  /**
   * Returns the elements between <code>from</code> and <code>to</code>
   * (both inclusive) in <code>result</code>.
   * An empty vector is returned if <code>from</code> is larger than
   * <code>to</code> (regardless how <code>result</code> is set).
   * If <code>result</code> is <code>null</code> a new vector is created.
   * Otherwise, the length of <code>result</code> has to be equal to
   * <code>to-from+1</code>.
   * So <code>getVec(0,length()-1,result)</code> is the same as
   * <code>copy(result)</code>.
   *
   * @param from index of the first element
   * @param to index of the last element
   * @param result storage of the result or <code>null</code> (new vector)
   * @return elements between <code>from</code> and <code>to</code>
   *         (both inclusive)
   */
  public final Vector getVec(int from, int to, Vector result) {
    assert (0 <= from && from < length());
    assert (0 <= to && to < length());
    if (from > to) { return EMPTY; }
    if (result == null) { result = create(to-from+1); }
    assert(result.length() == to-from+1);
    for (int k1 = from, k2 = 0; k1 <= to; ++k1, ++k2) {
      result.set(k2, get(k1));
    }
    return result;
  }

  /**
   * Returns the elements between <code>from</code> and <code>to</code>
   * (both inclusive).
   * So <code>getVec(0,length()-1)</code> is the same as <code>copy()</code>.
   * An empty vector is returned if <code>from</code> is larger than
   * <code>to</code>.
   *
   * @param from index of the first element
   * @param to index of the last element
   * @return elements between <code>from</code> and <code>to</code>
   *         (both inclusive)
   */
  public final Vector getVec(int from, int to) {
    return getVec(from, to, null);
  }

  /**
   * Returns the elements starting by <code>from</code> (inclusive).
   * If <code>result</code> is <code>null</code>, all elements until the end are
   * placed into a new vector. Otherwise, the vector elements are placed into
   * <code>result</code>, which has to have the same length as <code>this</code>
   * vector.
   *
   * @param from index of the first element
   * @param result storage of the result or <code>null</code> (new vector)
   * @return elements starting by <code>from</code>
   */
  public final Vector getVec(int from, Vector result) {
    return getVec(from, length()-1, result);
  }

  /**
   * Returns the elements starting by <code>from</code> (inclusive).
   *
   * @param from index of the first element
   * @return elements starting by <code>from</code>
   */
  public final Vector getVec(int from) {
    return getVec(from, null);
  }

  //----------------------------------------------------------------------------

  /**
   * Entrywise absolute value operation placed into <code>result</code>.
   *
   * Vector <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the absolute values
   */
  public Vector abs(Vector result) {
    assert (result != null && result.length() == length());
    for (int i = 0; i < length(); ++i) {
      result.set(i, Math.abs(get(i)));
    }
    return result;
  }

  @Override
  public Vector abs() {
    return abs(create(length()));
  }

  //----------------------------------------------------------------------------

  /**
   * Entrywise sign operation with zero replacement
   * placed into <code>result</code>.
   * 
   * Vector <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param zeroReplacement value replacing 0.0 values
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the sign and zero replacement values
   */
  public Vector sign(double zeroReplacement, Vector result) {
    assert (result != null && result.length() == length());
    double value;
    for (int i = 0; i < length(); ++i) {
      value = get(i);
      result.set(i, (0.0 == value) ? zeroReplacement : Math.signum(value));
    }
    return result;
  }

  /**
   * Entrywise sign operation placed into <code>result</code>.
   * 
   * Vector <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the sign values
   */
  public Vector sign(Vector result) {
    assert (result != null && result.length() == length());
    for (int i = 0; i < length(); ++i) {
      result.set(i, Math.signum(get(i)));
    }
    return result;
  }

  @Override
  public Vector sign(double zeroReplacement) {
    return sign(zeroReplacement, create(length()));
  }

  @Override
  public Vector sign() {
    return sign(create(length()));
  }

  //----------------------------------------------------------------------------

  /**
   * Entrywise negation placed into <code>result</code>.
   *
   * Vector <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the negated vector
   */
  public Vector neg(Vector result) {
    assert (result != null && result.length() == length());
    for (int i = 0; i < length(); ++i) { result.set(i, -get(i)); }
    return result;
  }

  @Override
  public Vector neg() {
    return neg(create(length()));
  }

  //----------------------------------------------------------------------------

  /**
   * Vector-vector addition (in <code>result</code>).
   * Adds vector <code>v</code> to <code>this</code> vector
   * in <code>result</code>.
   *
   * Vectors <code>v</code> and <code>result</code> have to have the same size
   * as <code>this</code>.
   * Vector <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * or <code>v</code> providing in-place operation.
   *
   * @param v vector to add (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the sum of
   *         <code>this</code> and <code>v</code>
   */
  public Vector add(Vector v, Vector result) {
    assert (v != null && v.length() == length());
    assert (result != null && result.length() == length());
    for (int i = 0; i < length(); ++i) {
      result.set(i, get(i) + v.get(i));
    }
    return result;
  }

  /**
   * Vector-vector addition (in new vector).
   * Adds vector <code>v</code> to <code>this</code> vector in a new vector.
   *
   * Vector <code>v</code> has to have the same size as <code>this</code>.
   *
   * @param v vector to add (not <code>null</code>)
   * @return the sum of <code>this</code> and <code>v</code> in a new vector
   */
  public Vector add(Vector v) {
    return add(v, create(length()));
  }

  /**
   * Vector-constant addition (in <code>result</code>).
   * Adds constant <code>c</code> to all elements of <code>this</code> vector
   * in <code>result</code>.
   *
   * Vector <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param c constant to add
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the shifted values of
   *         <code>this</code> by <code>c</code>
   */
  public Vector add(double c, Vector result) {
    assert (result != null && result.length() == length());
    for (int i = 0; i < length(); ++i) { result.set(i, get(i) + c); }
    return result;
  }

  @Override
  public Vector add(double c) {
    return add(c, create(length()));
  }

  //----------------------------------------------------------------------------

  /**
   * Vector-vector subtraction (in <code>result</code>).
   * Subtracts vector <code>v</code> from <code>this</code> vector
   * in <code>result</code>.
   *
   * Vectors <code>v</code> and <code>result</code> have to have the same size
   * as <code>this</code>.
   * Vector <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * or <code>v</code> providing in-place operation.
   *
   * @param v vector to subtract (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the difference of
   *         <code>this</code> and <code>v</code>
   */
  public Vector sub(Vector v, Vector result) {
    assert (v != null && v.length() == length());
    assert (result != null && result.length() == length());
    for (int i = 0; i < length(); ++i) {
      result.set(i, get(i) - v.get(i));
    }
    return result;
  }

  /**
   * Vector-vector subtraction (in new vector).
   * Subtracts vector <code>v</code> from <code>this</code> vector
   * in a new vector.
   *
   * Vector <code>v</code> has to have the same size as <code>this</code>.
   *
   * @param v vector to subtract (not <code>null</code>)
   * @return the difference of <code>this</code> and <code>v</code>
   *         in a new vector
   */
  public Vector sub(Vector v) {
    return sub(v, create(length()));
  }

  /**
   * Vector-constant subtraction (in <code>result</code>).
   * Subtracts constant <code>c</code> from all elements of <code>this</code>
   * vector in <code>result</code>.
   *
   * Vector <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param c constant to subtract
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the shifted values of
   *         <code>this</code> by <code>-c</code>
   */
  public Vector sub(double c, Vector result) {
    return add(-c, result);
  }

  @Override
  public Vector sub(double c) {
    return add(-c);
  }

  //----------------------------------------------------------------------------

  /**
   * Vector-constant division (in <code>result</code>). Divides all elements of
   * <code>this</code> vector by constant <code>c</code>.
   *
   * Vector <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param c constant to divide with
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the values of
   *         <code>this</code> divided by <code>c</code>
   */
  public Vector div(double c, Vector result) {
    return mul(1.0 / c, result);
  }

  @Override
  public Vector div(double c) {
    return div(c, create(length()));
  }

  //----------------------------------------------------------------------------

  /**
   * Vector-constant multiplication (in <code>result</code>). Multiplies all
   * elements of <code>this</code> vector by constant <code>c</code>.
   *
   * Vector <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param c constant to multiply with
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the values of
   *         <code>this</code> multiplied by <code>c</code>
   */
  public Vector mul(double c, Vector result) {
    assert (result != null && result.length() == length());
    for (int i = 0; i < length(); ++i) { result.set(i, c * get(i)); }
    return result;
  }

  @Override
  public Vector mul(double c) {
    return mul(c, create(length()));
  }

  /**
   * Vector-matrix multiplication (in <code>result</code>).
   * Multiplying the (transpose of the column) vector <code>this</code>
   * by matrix <code>m</code> from the right.
   *
   * Matrix <code>m</code> has to have the same number of rows
   * as the length of <code>this</code> vector. Furthermore, the length of
   * <code>result</code> has to be equal to the column number of <code>m</code>.
   *
   * @param m matrix multiplier from the right (not <code>null</code>)
   * @param result storage of the result
   *               (not <code>null</code> and not equal to <code>this</code>)
   * @return <code>result</code> being vector <code>v</code> multiplied by
   *         matrix <code>m</code> from the right
   * @see Matrix#mul(Vector, Vector)
   */    
  public Vector mul(Matrix m, Vector result) {
    assert (m != null && length() == m.rows());
    assert (result != null && result != this && result.length() == m.cols());
    double vi = get(0); // i = 0
    for (int j = 0; j < m.cols(); ++j) { // initialize "result"
      result.set(j, vi * m.get(0,j));
    }
    for (int i = 1; i < length(); ++i) {
      vi = get(i);
      for (int j = 0; j < m.cols(); ++j) {
        result.set(j, result.get(j) + vi * m.get(i,j));
      }
    }
    return result;
  }

  /**
   * Vector-matrix multiplication (in new vector).
   * Multiplying the (transpose of the column) vector <code>this</code>
   * by matrix <code>m</code> from the right.
   *
   * Matrix <code>m</code> has to have the same number of rows
   * as the length of <code>this</code> vector.
   *
   * @param m matrix multiplier from the right (not <code>null</code>)
   * @return new vector being the result of vector <code>v</code> multiplied by
   *         matrix <code>m</code> from the right
   * @see Matrix#mul(Vector)
   */    
  public Vector mul(Matrix m) {
    return mul(m, create(m.cols()));
  }

  /**
   * Quadratic vector-matrix-vector product.
   *
   * Matrix <code>m</code> has to be square with its row/column number
   * equal to the length of <code>this</code> vector.
   *
   * @param m matrix multiplier (not <code>null</code>)
   * @return matrix <code>m</code> multiplied by <code>this</code> from both sides
   */
  public final double mulQ(Matrix m) {
    assert (m != null && m.rows() == length() && m.cols() == length());
    double r = 0.0;
    for (int i = 0; i < length(); ++i) {
      for (int j = 0; j < length(); ++j) {
        r += get(i) * get(j) * m.get(i,j);
      }
    }
    return r;
  }

  /**
   * Multiplies matrix <code>m</code> from the left
   * with a diagonal matrix represented by <code>this</code> vector
   * (in <code>result</code>).
   *
   * Matrix <code>m</code> has to have the same number of rows as the length of
   * <code>this</code> vector.
   * Furthermore, <code>result</code> has to have the same size
   * as matrix <code>m</code>. The <code>result</code> parameter can also be set
   * to <code>m</code> providing in-place operation.
   *
   * @param m matrix multiplier (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> being matrix <code>m</code> multiplied
   *         with a diagonal matrix represented by <code>this</code> vector
   * @see Matrix#mulD(Vector, Matrix)
   */
  public Matrix mulD(Matrix m, Matrix result) {
    assert (m != null && length() == m.rows());
    assert (result != null && result.rows() == m.rows() && result.cols() == m.cols());
    for (int i = 0; i < length(); ++i) {
      double d = get(i);
      for (int j = 0; j < m.cols(); ++j) {
        result.set(i, j, m.get(i,j) * d);
      }
    }
    return result;
  }

  /**
   * Multiplies matrix <code>m</code> from the left
   * with a diagonal matrix represented by <code>this</code> vector
   * (in new matrix).
   *
   * Matrix <code>m</code> has to have the same number of rows as the length of
   * <code>this</code> vector.
   *
   * @param m matrix multiplier (not <code>null</code>)
   * @return new matrix being the result of matrix <code>m</code> multiplied
   *         with a diagonal matrix represented by <code>this</code> vector
   * @see Matrix#mulD(Vector)
   */
  public Matrix mulD(Matrix m) {
    return mulD(m, Matrix.create(m.rows(), m.cols()));
  }

  /**
   * Permutes the elements of the vector, i.e. multiplying the vector
   * (left or right) with a permutation matrix represented by <code>p</code>
   * (in <code>result</code>).
   *
   * The length of permutation <code>p</code> has to be equal to the length
   * of <code>this</code> vector.
   * Furthermore, <code>result</code> has to have the same length
   * as <code>this</code> vector.
   *
   * @param p permutation (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>
   *                                      and not equal to <code>this</code>)
   * @return <code>result</code>
   *         having <code>this</code> vector permuted by <code>p</code>
   * @see Permutation#mul(Vector, Vector)
   */
  public Vector mul(Permutation p, Vector result) {
    assert (p != null && p.length() == length());
    assert (result != null && result != this && result.length() == length());
    for (int i = 0; i < length(); ++i) {
      result.set(i, get(p.get(i)));
    }
    return result;
  }

  /**
   * Permutes the elements of the vector, i.e. multiplying the vector
   * (left or right) with a permutation matrix represented by <code>p</code>
   * (in new vector).
   *
   * The length of permutation <code>p</code> has to be equal to the length
   * of <code>this</code> vector.
   *
   * @param p permutation (not <code>null</code>)
   * @return new vector being equal to
   *         <code>this</code> vector permuted by <code>p</code>
   * @see Permutation#mul(Vector)
   */
  public Vector mul(Permutation p) {
    return mul(p, create(length()));
  }

  //----------------------------------------------------------------------------

  /**
   * Takes the reciproc of all elements (in <code>result</code>).
   *
   * Vector <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> holding the elementwise reciproc vector
   */
  public Vector reciproc(Vector result) {
    assert (result != null && result.length() == length());
    for (int i = 0; i < length(); ++i) {
      result.set(i, 1.0 / get(i));
    }
    return result;
  }

  @Override
  public Vector reciproc() {
    return reciproc(create(length()));
  }

  //----------------------------------------------------------------------------

  /**
   * Entrywise (signed) remainder with respect to modulus <code>m</code>
   * (in <code>result</code>).
   *
   * Vector <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param m modulus
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the remainders of the values of
   *         <code>this</code> with respect to modulus <code>m</code>
   */
  public Vector mod(double m, Vector result) {
    assert (result != null && result.length() == length());
    for (int i = 0; i < length(); ++i) { result.set(i, get(i) % m); }
    return result;
  }

  @Override
  public Vector mod(double m) {
    return mod(m, create(length()));
  }

  //----------------------------------------------------------------------------

  /**
   * Inner product (dot product, scalar product) with vector <code>v</code>.
   *
   * The length of <code>this</code> and <code>v</code> has to be the same.
   *
   * @param v vector to multiply with
   * @return inner product of <code>this</code> and vector <code>v</code>
   */
  public final double inner(Vector v) {
    assert (v != null && v.length() == length());
    double r = 0.0;
    for (int i = 0; i < length(); ++i) { r += get(i) * v.get(i); }
    return r;
  }

  /**
   * Outer product with vector <code>v</code> (in <code>result</code>).
   *
   * The row and column number of matrix <code>result</code> has to be equal to
   * the length <code>this</code> and <code>v</code>, respectively.
   *
   * @param v vector to multiply with (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> being the
   *         outer product of <code>this</code> and vector <code>v</code>
   */
  public Matrix outer(Vector v, Matrix result) {
    assert (result != null && result.rows() == length());
    assert (v != null && result.cols() == v.length());
    for (int i = 0; i < length(); ++i) {
      for (int j = 0; j < v.length(); ++j) {
        result.set(i, j, get(i) * v.get(j));
      }
    }
    return result;
  }

  /**
   * Outer product with vector <code>v</code> (in new matrix).
   *
   * @param v vector to multiply with (not <code>null</code>)
   * @return new square matrix being the
   *         outer product of <code>this</code> and vector <code>v</code>
   */
  public Matrix outer(Vector v) {
    return outer(v, Matrix.create(length(), v.length()));
  }

  //----------------------------------------------------------------------------

  /**
   * Entrywise multiplication (Hadamard product) by vector <code>v</code>
   * (in <code>result</code>).
   *
   * Vectors <code>v</code> and <code>result</code> have to have the same size
   * as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param v vector to multiply with entrywise (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the elements of <code>this</code>
   *         and <code>v</code> multiplied entrywise
   */
  public Vector emul(Vector v, Vector result) {
    assert (v != null && v.length() == length());
    assert (result != null && result.length() == length());
    for (int i = 0; i < length(); ++i) {
      result.set(i, get(i) * v.get(i));
    }
    return result;
  }

  /**
   * Entrywise multiplication (Hadamard product) by vector <code>v</code>
   * (in new vector).
   *
   * Vector <code>v</code> has to have the same size as <code>this</code>.
   *
   * @param v vector to multiply with entrywise (not <code>null</code>)
   * @return new vector having the elements of <code>this</code>
   *         and <code>v</code> multiplied entrywise
   */
  public Vector emul(Vector v) {
    return emul(v, create(length()));
  }

  //----------------------------------------------------------------------------

  /**
   * Normalize the vector with respect to the 1-norm.
   * The zero vector is left unchanged.
   *
   * @return (1-norm) normalized <code>this</code>
   */
  public Vector normalize1() {
    double norm1 = norm1();
    if (0.0 < norm1) { div(norm1, this); }
    return this;
  }

  /**
   * Normalize the vector with respect to the 2-norm.
   * The zero vector is left unchanged.
   *
   * @return (2-norm) normalized <code>this</code>
   */
  public Vector normalize2() {
    double norm2 = norm2();
    if (0.0 < norm2) { div(norm2, this); }
    return this;
  }

  /**
   * Normalize the vector with respect to the inf-norm.
   * The zero vector is left unchanged.
   *
   * @return (inf-norm) normalized <code>this</code>
   */
  public Vector normalizeI() {
    double normI = normI();
    if (0.0 < normI) { div(normI, this); }
    return this;
  }

  /**
   * Computes the 1-norm of <code>this</code> vector.
   *
   * @return 1-norm of <code>this</code> vector
   */
  public double norm1() {
    double s = 0.0;
    for (int i = 0; i < length(); ++i) { s += Math.abs(get(i)); }
    return s;
  }

  /**
   * Computes the 2-norm of <code>this</code> vector.
   *
   * @return 2-norm of <code>this</code> vector
   */
  public double norm2() {
    double s = 0.0, v;
    for (int i = 0; i < length(); ++i) {
      v = get(i);
      s += v * v;
    }
    return Math.sqrt(s);
  }

  /**
   * Computes the infinity-norm of <code>this</code> vector.
   *
   * @return infinity-norm of <code>this</code> vector
   */
  public double normI() {
    double s = 0.0, v;
    for (int i = 0; i < length(); ++i) {
      v = Math.abs(get(i));
      if (v > s) { s = v; }
    }
    return s;
  }

  //----------------------------------------------------------------------------

  /**
   * Returns the sum of the vector elements.
   *
   * @return sum of vector elements
   */
  public double sum() {
    double s = 0.0;
    for (int i = 0; i < length(); ++i) { s += get(i); }
    return s;
  }

  /**
   * Returns the product of the vector elements.
   *
   * @return product of vector elements
   */
  public double prod() {
    double p = 1.0;
    for (int i = 0; i < length(); ++i) { p *= get(i); }
    return p;
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

  //----------------------------------------------------------------------------

  private final double[] data;
}
