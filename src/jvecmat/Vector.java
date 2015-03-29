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
    Vector v = Vector.create(length);
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
    Vector v = Vector.create(length);
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
    Vector v = Vector.create(length);
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
   * Returns the array representation of the vector.
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
   *
   * @param result appropriately sized storage for the copy
   *               or <code>null</code>
   * @return copy of the vector
   */
  public Vector copy(Vector result) {
    if (result == null) { result = Vector.create(length()); }
    assert(result.length() == length());
    for (int i = 0; i < length(); ++i) { result.set(i, get(i)); }
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
    if (result == null) { result = Vector.create(to-from+1); }
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
   * Normalize (1-norm) the vector if it is not the zero vector.
   * @return (1-norm) normalized "this"
   */
  public Vector normalize1() {
    double norm1 = norm1();
    if (0.0 < norm1) { divL(norm1); }
    return this;
  }

  /**
   * Normalize (2-norm) the vector if it is not the zero vector.
   * @return (2-norm) normalized "this"
   */
  public Vector normalize2() {
    double norm2 = norm2();
    if (0.0 < norm2) { divL(norm2); }
    return this;
  }

  /**
   * Normalize (inf-norm) the vector if it is not the zero vector.
   * @return (inf-norm) normalized "this"
   */
  public Vector normalizeI() {
    double normI = normI();
    if (0.0 < normI) { divL(normI); }
    return this;
  }

  //----------------------------------------------------------------------------

  /**
   * Entrywise absolute value.
   * @return |this| (placed into "result")
   */
  public <T extends Vector> T abs(T result) {
    assert (result.length() == length());
    for (int i = 0; i < length(); ++i) {
      result.set(i, Math.abs(get(i)));
    }
    return result;
  }

  /**
   * Entrywise absolute value.
   * @return |this| (placed into a new vector)
   */
  public Vector abs() {
    return abs(create(length()));
  }

  /**
   * Entry absolute value.
   * @return |this| (placed into "this")
   */
  public Vector absL() {
    for (int i = 0; i < length(); ++i) {
      if (0.0 > get(i)) { set(i, -get(i)); }
    }
    return this;
  }

  //----------------------------------------------------------------------------

  protected final double sign(double value, double zeroReplacement) {
    return (0.0 == value) ? zeroReplacement : Math.signum(value);
  }

  /**
   * Entrywise sign operation with zero replacement.
   * @return sign(this) (placed into "result")
   */
  public <T extends Vector> T sign(double zeroReplacement, T result) {
    assert (result.length() == length());
    for (int i = 0; i < length(); ++i) {
      result.set(i, sign(get(i), zeroReplacement));
    }
    return result;
  }

  /**
   * Entrywise sign operation.
   * @return sign(this) (placed into "result")
   */
  public <T extends Vector> T sign(T result) {
    assert (result.length() == length());
    for (int i = 0; i < length(); ++i) {
      result.set(i, Math.signum(get(i)));
    }
    return result;
  }

  /**
   * Entrywise sign operation with zero replacement.
   * @return sign(this) (placed into a new vector)
   */
  public Vector sign(double zeroReplacement) {
    return sign(zeroReplacement, create(length()));
  }

  /**
   * Entrywise sign operation.
   * @return sign(this) (placed into a new vector)
   */
  public Vector sign() {
    return sign(create(length()));
  }

  /**
   * Entrywise sign operation with zero replacement.
   * @return sign(this) (placed into "this")
   */
  public Vector signL(double zeroReplacement) {
    return sign(zeroReplacement, this);
  }

  /**
   * Entrywise sign operation.
   * @return sign(this) (placed into "this")
   */
  public Vector signL() {
    return sign(this);
  }

  //----------------------------------------------------------------------------

  /**
   * @return -this (placed into "result")
   */
  public <T extends Vector> T neg(T result) {
    assert (result.length() == length());
    for (int i = 0; i < length(); ++i) { result.set(i, -get(i)); }
    return result;
  }

  /**
   * @return -this (placed into a new vector)
   */
  public Vector neg() {
    return neg(create(length()));
  }

  /**
   * @return -this (placed into "this")
   */
  public Vector negL() {
    return neg(this);
  }

  //----------------------------------------------------------------------------

  /**
   * @return this + v (placed into "result")
   */
  public <T extends Vector> T add(Vector v, T result) {
    assert (length() == v.length());
    assert (result.length() == length());
    for (int i = 0; i < length(); ++i) {
      result.set(i, get(i) + v.get(i));
    }
    return result;
  }

  /**
   * @return this + v (placed into a new vector)
   */
  public Vector add(Vector v) {
    return add(v, Vector.create(length()));
  }

  /**
   * @return this + v (placed into "this")
   */
  public Vector addL(Vector v) {
    return add(v, this);
  }

  /**
   * @return this + v (placed into "v")
   */
  public <T extends Vector> T addR(T v) {
    return add(v, v);
  }

  //----------------------------------------------------------------------------

  /**
   * @return this + c * Vector.one (placed into "result")
   */
  public <T extends Vector> T add(double c, T result) {
    assert (result.length() == length());
    for (int i = 0; i < length(); ++i) { result.set(i, get(i) + c); }
    return result;
  }

  /**
   * @return this + c * Vector.one (placed into a new vector)
   */
  public Vector add(double c) {
    return add(c, Vector.create(length()));
  }

  /**
   * @return this + c * Vector.one (placed into "this")
   */
  public Vector addL(double c) {
    return add(c, this);
  }

  //----------------------------------------------------------------------------

  /**
   * @return this - v (placed into "result")
   */
  public <T extends Vector> T sub(Vector v, T result) {
    assert (length() == v.length());
    assert (result.length() == length());
    for (int i = 0; i < length(); ++i) {
      result.set(i, get(i) - v.get(i));
    }
    return result;
  }

  /**
   * @return this - v (placed into a new vector)
   */
  public Vector sub(Vector v) {
    return sub(v, Vector.create(length()));
  }

  /**
   * @return this - v (placed into "this")
   */
  public Vector subL(Vector v) {
    return sub(v, this);
  }

  /**
   * @return this - v (placed into "v")
   */
  public <T extends Vector> T subR(T v) {
    return sub(v, v);
  }

  //----------------------------------------------------------------------------

  /**
   * @return this - c * Vector.one (placed into "result")
   */
  public <T extends Vector> T sub(double c, T result) {
    return add(-c, result);
  }

  /**
   * @return this - c * Vector.one (placed into a new vector)
   */
  public Vector sub(double c) {
    return add(-c);
  }

  /**
   * @return this - c * Vector.one (placed into "this")
   */
  public Vector subL(double c) {
    return addL(-c);
  }

  //----------------------------------------------------------------------------

  /**
   * @return this * c (placed into "result")
   */
  public <T extends Vector> T mul(double c, T result) {
    assert (result.length() == length());
    for (int i = 0; i < length(); ++i) { result.set(i, c * get(i)); }
    return result;
  }

  /**
   * @return this * c (placed into a new vector)
   */
  public Vector mul(double c) {
    return mul(c, Vector.create(length()));
  }

  /**
   * @return this * c (placed into "this")
   */
  public Vector mulL(double c) {
    return mul(c, this);
  }

  //----------------------------------------------------------------------------

  /**
   * @return this / c (placed into "result")
   */
  public <T extends Vector> T div(double c, T result) {
    return mul(1.0 / c, result);
  }

  /**
   * @return this / c (placed into a new vector)
   */
  public Vector div(double c) {
    return div(c, Vector.create(length()));
  }

  /**
   * @return this / c (placed into "this")
   */
  public Vector divL(double c) {
    return div(c, this);
  }

  //----------------------------------------------------------------------------

  /**
   * @return this .% m (placed into "result")
   */
  public <T extends Vector> T mod(double m, T result) {
    assert (result.length() == length());
    for (int i = 0; i < length(); ++i) { result.set(i, get(i) % m); }
    return result;
  }

  /**
   * @return this .% m (placed into a new vector)
   */
  public Vector mod(double m) {
    return mod(m, create(length()));
  }

  /**
   * @return this .% m (placed into "this")
   */
  public Vector modL(double m) {
    return mod(m, this);
  }
    
  //----------------------------------------------------------------------------

  /**
   * @return inner product (dot product), <this,v>
   */
  public final double iprod(Vector v) {
    assert (length() == v.length());
    double r = 0.0;
    for (int i = 0; i < length(); ++i) { r += get(i) * v.get(i); }
    return r;
  }

  //----------------------------------------------------------------------------

  /**
   * Entrywise multiplication.
   * @return this .* v (placed into "result")
   */
  public <T extends Vector> T emul(Vector v, T result) {
    assert (length() == v.length());
    assert (result.length() == length());
    for (int i = 0; i < length(); ++i) {
      result.set(i, get(i) * v.get(i));
    }
    return result;
  }

  /**
   * Entrywise multiplication.
   * @return this .* v (placed into a new vector)
   */
  public Vector emul(Vector v) {
    return emul(v, Vector.create(length()));
  }

  /**
   * Entrywise multiplication.
   * @return this .* v (placed into "this") 
   */
  public Vector emulL(Vector v) {
    return emul(v, this);
  }

  /**
   * Entrywise multiplication.
   * @return this .* v (placed into "v")
   */
  public <T extends Vector> T emulR(T v) {
    return emul(v, v);
  }

  //----------------------------------------------------------------------------

  /**
   * Vector-matrix multiplication.
   * The "result" parameter has to be different from "this".
   * @return (this^T * m)^T (placed into "result")
   */    
  public <T extends Vector> T mul(Matrix m, T result) {
    assert (length() == m.rows());
    assert (result != this);
    assert (result.length() == m.cols());
    final int cols = m.cols();
    int i, j;
    double vi = get(0); // i = 0
    for (j = 0; j < cols; ++j) { // initialize "result"
      result.set(j, vi * m.get(0,j));
    }
    for (i = 1; i < length(); ++i) {
      vi = get(i);
      for (j = 0; j < cols; ++j) {
        result.set(j, result.get(j) + vi * m.get(i,j));
      }
    }
    return result;
  }

  /**
   * Vector-matrix multiplication.
   * @return (this^T * m)^T (placed into a new vector)
   */
  public Vector mul(Matrix m) {
    return mul(m, Vector.create(m.cols()));
  }

  //----------------------------------------------------------------------------

  /**
   * Quadratic vector-matrix-vector product.
   * @return this^T * m * this
   */
  public final double mulQ(Matrix m) {
    assert (m.rows() == length());
    assert (m.cols() == length());
    int i, j;
    double r = 0.0;
    for (i = 0; i < length(); ++i) {
      for (j = 0; j < length(); ++j) {
        r += get(i) * get(j) * m.get(i,j);
      }
    }
    return r;
  }

  //----------------------------------------------------------------------------

  /**
   * Matrix product with a diagonal matrix represented by "this".
   * @return diag(this) * m (placed into "result")
   */
  public Matrix mulD(Matrix m, Matrix result) {
    assert (length() == m.rows());
    assert (result.rows() == m.rows());
    assert (result.cols() == m.cols());
    final int cols = m.cols();
    int i, j;
    double d;
    for (i = 0; i < length(); ++i) {
      d = get(i);
      for (j = 0; j < cols; ++j) {
        result.set(i, j, m.get(i,j) * d);
      }
    }
    return result;
  }

  /**
   * Matrix product with a diagonal matrix represented by "this".
   * @return diag(this) * m (placed into a new matrix)
   */
  public Matrix mulD(Matrix m) {
    return mulD(m, Matrix.create(m.rows(),m.cols()));
  }

  /**
   * @return determinant of the diag("this") matrix
   */
  public double detD() {
    double det = 1.0;
    for (int i = 0; i < length(); ++i) { det *= get(i); }
    return det;
  }

  //----------------------------------------------------------------------------

  /**
   * Outer product of two vectors.
   * @return this * v^T (placed into "result")
   */
  public Matrix outp(Vector v, Matrix result) {
    assert (result.rows() == length());
    assert (result.cols() == v.length());
    int i, j;
    for (i = 0; i < length(); ++i) {
      for (j = 0; j < v.length(); ++j) {
        result.set(i, j, get(i) * v.get(j));
      }
    }
    return result;
  }

  /**
   * Outer product of two vectors.
   * @return this * v^T (placed into a new matrix)
   */
  public Matrix outp(Vector v) {
    return outp(v, Matrix.create(length(), v.length()));
  }

  //----------------------------------------------------------------------------

  /**
   * Take the reciproc of all elements.
   * It is assumed that the elements are non-zero.
   * @return 1 ./ this (placed into "result")
   */
  public <T extends Vector> T reciproc(T result) {
    assert (result.length() == length());
    for (int i = 0; i < length(); ++i) {
      result.set(i, 1.0 / get(i));
    }
    return result;
  }

  /**
   * Take the reciproc of all elements.
   * It is assumed that the elements are non-zero.
   * @return 1 ./ this (placed into a new vector)
   */
  public Vector reciproc() {
    return reciproc(Vector.create(length()));
  }

  /**
   * Take the reciproc of all elements.
   * It is assumed that the elements are non-zero.
   * @return 1 ./ this (placed into "this")
   */
  public Vector reciprocL() {
    return reciproc(this);
  }

  //----------------------------------------------------------------------------

  /**
   * @return 1-norm of the vector
   */
  public double norm1() {
    double s = 0.0;
    for (int i = 0; i < length(); ++i) { s += Math.abs(get(i)); }
    return s;
  }

  /**
   * @return 2-norm of the vector
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
   * @return inf-norm of the vector
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

  public String logString() {
    String str = "";
    for (int i = 0; i < length(); ++i) {
      if (i != 0) { str += " "; }
      str += get(i);
    }
    return str;
  }

  @Override
  public String toString() {
    return "[" + logString() + "]";
  }

  //----------------------------------------------------------------------------

  private final double[] data;
}
