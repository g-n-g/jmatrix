package jvecmat;

import java.util.Random;

/**
 * Abstract matrix representation.
 */
public abstract class Matrix {

  /**
   * Creates a matrix using the provided data.
   *
   * @param data elements of the vector
   * @return matrix which encapsulates <code>data</code> (not copied)
   */
  public static Matrix create(double[][] data) {
    return new NonTransposedMatrix(data);
  }

  /**
   * Creates a matrix of given size having uninitialized elements.
   *
   * @param rows number of rows
   * @param cols number of columns
   * @return matrix with size <code>rows</code> x <code>cols</code>
   *         having uninitialized elements
   */
  public static Matrix create(int rows, int cols) {
    return create(new double[rows][cols]);
  }

  /**
   * Creates a matrix by its <code>rows</code>.
   *
   * @param rows vectors forming the rows of the matrix
   * @return matrix formed by <code>rows</code> (not copied)
   */
  public static Matrix createByRows(Vector[] rows) {
    double [][]d = new double[rows.length][];
    for (int i = 0; i < rows.length; ++i) { d[i] = rows[i].array(); }
    return create(d);
  }

  /**
   * Creates a matrix by its <code>columns</code>.
   *
   * @param columns vectors forming the columns of the matrix
   * @return matrix formed by <code>columns</code> (not copied)
   */
  public static Matrix createByCols(Vector[] columns) {
    return createByRows(columns).T();
  }

  /**
   * Creates a matrix of size <code>rows</code> x <code>cols</code>
   * and initializes its elements to <code>values</code>.
   *
   * @param rows number of rows
   * @param cols number of columns
   * @param value value of elements
   * @return constant matrix of size "rows" x "cols" having elements "c"
   */
  public static Matrix constant(int rows, int cols, double value) {
    Matrix m = create(rows, cols);
    m.setToConstant(value);
    return m;
  }

  /**
   * Creates a matrix of size <code>rows</code> x <code>cols</code>
   * and initializes its elements to zero.
   *
   * @param rows number of rows
   * @param cols number of columns
   * @return zero matrix of size <code>rows</code> x <code>cols</code>
   */
  public static Matrix zero(int rows, int cols) {
    return constant(rows, cols, 0.0);
  }

  /**
   * Creates a matrix of size <code>rows</code> x <code>cols</code>
   * and initializes its elements to zero.
   *
   * @param rows number of rows
   * @param cols number of columns
   * @return all one matrix of size <code>rows</code> x <code>cols</code>
   */
  public static Matrix one(int rows, int cols) {
    return constant(rows, cols, 1.0);
  }

  /**
   * Creates a square diagonal matrix
   * using the elements of vector <code>v</code>.
   *
   * @param v vector forming the diagonal
   * @return square diagonal matrix defined by vector <code>v</code>
   */
  public static Matrix diag(Vector v) {
    int n = v.length();
    double[][] mat = new double[n][n];
    for (int i = 0; i < n; ++i) {
      mat[i][i] = v.get(i);
      for (int j = 0; j < i; ++j)
        mat[i][j] = mat[j][i] = 0.0;
    }
    return Matrix.create(mat);
  }

  /**
   * Creates an identity matrix of size <code>dim</code> x <code>dim</code>.
   *
   * @param dim dimension
   * @return identity matrix of size <code>dim</code> x <code>dim</code>
   */
  public static Matrix eye(int dim) {
    double[][] mat = new double[dim][dim];
    for (int i = 0; i < dim; ++i) {
      mat[i][i] = 1.0;
      for (int j = 0; j < i; ++j) { mat[i][j] = mat[j][i] = 0.0; }
    }
    return Matrix.create(mat);
  }

  /**
   * Creates a random matrix of size <code>rows</code> x <code>cols</code>
   * drawing its elements from the uniform distribution on [0,1].
   *
   * @param rows number of rows
   * @param cols number of columns
   * @param rng random number generator
   * @return uniform random matrix
   *         of size <code>rows</code> x <code>cols</code>
   */
  public static Matrix rand(int rows, int cols, Random rng) {
    Matrix m = Matrix.create(rows,cols);
    m.setToRand(rng);
    return m;
  }

  /**
   * Creates a random matrix of size <code>rows</code> x <code>cols</code>
   * drawing its elements from the standard normal distribution.
   *
   * @param rows number of rows
   * @param cols number of columns
   * @param rng random number generator
   * @return standard normal random matrix
   *         of size <code>rows</code> x <code>cols</code>
   */
  public static Matrix randN(int rows, int cols, Random rng) {
    Matrix m = Matrix.create(rows,cols);
    m.setToRandN(rng);
    return m;
  }

  //----------------------------------------------------------------------------

  /**
   * Creates a matrix object of size <code>rows</code> x <code>cols</code>
   * which encapsulates <code>data</code>.
   *
   * @param rows number of rows
   * @param cols number of columns
   * @param data matrix data (not copied)
   */
  protected Matrix(double[][] data, int rows, int cols) {
    this.data = data;
    this.rows = rows;
    this.cols = cols;
  }

  //----------------------------------------------------------------------------

  /**
   * @return true if there are NaN elements in the matrix
   */
  public boolean hasNaN() {
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        if (Double.isNaN(get(i,j))) { return true; }
      }
    }
    return false;
  }

  /**
   * @return true if there are Infinite elements in the matrix
   */
  public boolean hasInf() {
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        if (Double.isInfinite(get(i,j))) { return true; }
      }
    }
    return false;
  }

  /**
   * Replace NaN and infinite elements in the matrix by the given ones.
   */
  public void replaceNaNandInf(double nan, double negInf, double posInf) {
    int i, j;
    double e;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        e = get(i,j);
        if (Double.isNaN(e)) { set(i, j, nan); }
        else if (Double.POSITIVE_INFINITY == e) { set(i, j, posInf); }
        else if (Double.NEGATIVE_INFINITY == e) { set(i, j, negInf); }
      }
    }
  }

  //----------------------------------------------------------------------------

  /**
   * Copy "this" matrix into "result".
   * @return copied matrix (placed into "result")
   */
  public Matrix copy(Matrix result) {
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        result.set(i, j, get(i,j));
      }
    }
    return result;
  }

  /**
   * @return new matrix with copied data
   */
  public Matrix copy() {
    return copy(Matrix.create(rows(),cols()));
  }

  //----------------------------------------------------------------------------

  /**
   * @return number of rows of the matrix
   */
  public final int rows() {
    return rows;
  }

  /**
   * @return number of columns of the matrix
   */
  public final int cols() {
    return cols;
  }

  /**
   * @return representation array of the matrix
   */
  public final double[][] array() {
    return data;
  }

  /**
   * @return return the (i,j) element of the matrix
   */
  abstract public double get(int i, int j);

  /**
   * Set matrix element at ("i","j") to "value".
   */
  abstract public void set(int i, int j, double value);

  //----------------------------------------------------------------------------

  /**
   * Set all matrix elements to "c".
   * @return "this"
   */
  public Matrix setToConstant(double c) {
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        set(i, j, c);
      }
    }
    return this;
  }

  /**
   * Set all matrix elements to zero.
   * @return "this"
   */
  public Matrix setToZero() {
    return setToConstant(0.0);
  }

  /**
   * Set all matrix elements to one.
   * @return "this"
   */
  public Matrix setToOne() {
    return setToConstant(1.0);
  }

  /**
   * Set to the identity matrix.
   * @return "this"
   */
  public Matrix setToEye() {
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        set(i, j, i==j ? 1.0 : 0.0);
      }
    }
    return this;
  }

  /**
   * Set all matrix elements U(0,1) randomly.
   * @return "this"
   */
  public Matrix setToRand(Random rng) {
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        set(i, j, rng.nextDouble());
      }
    }
    return this;
  }

  /**
   * Set all matrix elements N(0,1) randomly.
   * @return "this"
   */
  public Matrix setToRandN(Random rng) {
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        set(i, j, rng.nextGaussian());
      }
    }
    return this;
  }

  //----------------------------------------------------------------------------

  /**
   * @return diagonal elements of a squared matrix (placed into "result")
   */
  public Vector getDiag(Vector result) {
    assert (rows() == cols());
    for (int i = 0; i < rows(); ++i) { result.set(i, get(i,i)); }
    return result;
  }

  /**
   * @return diagonal elements of a squared matrix (placed into a new vector)
   */
  public Vector getDiag() {
    return getDiag(Vector.create(rows()));
  }

  //----------------------------------------------------------------------------

  /**
   * @return i-th row of the matrix (placed into "result")
   */
  public Vector getRow(int i, Vector result) {
    assert (0 <= i && i <= rows());
    assert (result.length() == cols());
    for (int j = 0; j < cols(); ++j) { result.set(j, get(i,j)); }
    return result;
  }

  /**
   * @return i-th row of the matrix (placed into a new vector)
   */
  public Vector getRow(int i) {
    return getRow(i, Vector.create(cols()));
  }

  /**
   * Set "this"[i:0][i:v.length-1] from "v".
   * @return "this"
   */
  public Matrix setRow(int i, Vector v) {
    assert (0 <= i && i <= rows());
    assert (cols() == v.length());
    for (int j = 0; j < v.length(); ++j) { set(i, j, v.get(j)); }
    return this;
  }

  /**
   * Set "this"[i:0][i:m.cols-1] from "m".
   * @return "this"
   */
  public Matrix setRow(int i, Matrix m) {
    assert (0 <= i && i <= rows());
    assert (cols() == m.cols());
    for (int j = 0; j < m.cols(); ++j) { set(i, j, m.get(i, j)); }
    return this;
  }

  //----------------------------------------------------------------------------

  /**
   * @return j-th column of the matrix (placed into "result")
   */
  public Vector getCol(int j, Vector result) {
    assert (0 <= j && j <= cols());
    assert (result.length() == rows());
    for (int i = 0; i < rows(); ++i) { result.set(i, get(i,j)); }
    return result;
  }

  /**
   * @return j-th column of the matrix (placed into a new vector)
   */
  public Vector getCol(int j) {
    return getCol(j, Vector.create(rows()));
  }

  /**
   * Set "this"[0:j][v.length-1:j] from "v".
   * @return "this"
   */
  public Matrix setCol(int j, Vector v) {
    assert (0 <= j && j <= cols());
    assert (rows() == v.length());
    for (int i = 0; i < v.length(); ++i) { set(i, j, v.get(i)); }
    return this;
  }

  /**
   * Set "this"[0:j][m.rows-1:j] from "m".
   * @return "this"
   */
  public Matrix setCol(int j, Matrix m) {
    assert (0 <= j && j <= cols());
    assert (rows() == m.rows());
    for (int i = 0; i < m.rows(); ++i) { set(i, j, m.get(i, j)); }
    return this;
  }

  //----------------------------------------------------------------------------

  /**
   * Get the sub-matrix having rows iF-iT and columns jF-jT.
   * Both interval end-points are inclusive.
   * @return [iF:iT][jF:jT] sub-matrix (placed into "result")
   */
  public Matrix getMat(int iF, int iT, int jF, int jT, Matrix result) {
    assert (0 <= iF && iF <= iT && iT < rows());
    assert (0 <= jF && jF <= jT && jT < cols());
    assert (result.rows() == iT-iF+1);
    assert (result.cols() == jT-jF+1);
    int i, j, ir, jr;
    for (i = iF, ir = 0; i <= iT; ++i, ++ir) {
      for (j = jF, jr = 0; j <= jT; ++j, ++jr) {
        result.set(ir, jr, get(i,j));
      }
    }
    return result;
  }

  /**
   * Get the sub-matrix having rows iF-iT and columns jF-jT.
   * Both interval end-points are inclusive.
   * @return [iF:iT][jF:jT] sub-matrix (placed into a new matrix)
   */
  public Matrix getMat(int iF, int iT, int jF, int jT) {
    return getMat(iF, iT, jF, jT, create(iT-iF+1,jT-jF+1));
  }

  /**
   * Get the sub-matrix having rows iF-(iF+result.rows-1)
   * and columns jF-(jF+result.cols-1).
   * Both interval end-points are inclusive.
   * @return [iF:iF+result.rows-1][jF:jF+result.cols-1] sub-matrix
   *         (placed into "result")
   */
  public Matrix getMat(int iF, int jF, Matrix result) {
    assert (0 <= iF && iF + result.rows() <= rows());
    assert (0 <= jF && jF + result.cols() <= cols());
    for (int i = 0; i < result.rows(); ++i) {
      for (int j = 0; j < result.cols(); ++j) {
        result.set(i, j, get(iF+i, jF+j));
      }
    }
    return result;
  }

  /**
   * Set "this"[iFrom:iFrom+m.rows-1][jFrom:jFrom+m.cols-1] from "m".
   * @return "this"
   */
  public Matrix setMat(int iF, int jF, Matrix m) {
    assert (0 <= iF && iF+m.rows() <= rows());
    assert (0 <= jF && jF+m.cols() <= cols());
    for (int i = 0; i < m.rows(); ++i) {
      for (int j = 0; j < m.cols(); ++j) {
        set(iF+i, jF+j, m.get(i, j));
      }
    }
    return this;
  }

  //----------------------------------------------------------------------------

  /**
   * Entrywise absolute value.
   * @return |this| (placed into "result")
   */
  public Matrix abs(Matrix result) {
    assert (result.rows() == rows());
    assert (result.cols() == cols());
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        result.set(i, j, Math.abs(get(i,j)));
      }
    }
    return result;
  }

  /**
   * Entrywise absolute value.
   * @return |this| (placed into a new matrix)
   */
  public Matrix abs() {
    return abs(create(rows(),cols()));
  }

  /**
   * Entrywise absolute value.
   * @return |this| (placed into "this")
   */
  public Matrix absL() {
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        if (0.0 > get(i,j)) { set(i, j, -get(i,j)); }
      }
    }
    return this;
  }

  //----------------------------------------------------------------------------

  /**
   * Entrywise sign operation with zero replacement.
   * @return sign(this) (placed into "result")
   */
  public Matrix sign(double zeroReplacement, Matrix result) {
    assert (result.rows() == rows());
    assert (result.cols() == cols());
    int i, j;
    double e;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        e = get(i,j);
        result.set(i, j, (0.0 == e) ? zeroReplacement : Math.signum(e));
      }
    }
    return result;
  }

  /**
   * Entrywise sign operation.
   * @return sign(this) (placed into "result")
   */
  public Matrix sign(Matrix result) {
    assert (result.rows() == rows());
    assert (result.cols() == cols());
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        result.set(i, j, Math.signum(get(i,j)));
      }
    }
    return result;
  }

  /**
   * Entrywise sign operation with zero replacement.
   * @return sign(this) (placed into a new matrix)
   */
  public Matrix sign(double zeroReplacement) {
    return sign(zeroReplacement, create(rows(),cols()));
  }

  /**
   * Entrywise sign operation.
   * @return sign(this) (placed into a new matrix)
   */
  public Matrix sign() {
    return sign(create(rows(),cols()));
  }

  /**
   * Entrywise sign operation with zero replacement.
   * @return sign(this) (placed into "this")
   */
  public Matrix signL(double zeroReplacement) {
    return sign(zeroReplacement, this);
  }

  /**
   * Entrywise sign operation.
   * @return sign(this) (placed into "this")
   */
  public Matrix signL() {
    return sign(this);
  }

  //----------------------------------------------------------------------------

  /**
   * @return -this (placed into "result")
   */
  public Matrix neg(Matrix result) {
    assert (result.rows() == rows());
    assert (result.cols() == cols());
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        result.set(i, j, -get(i,j));
      }
    }
    return result;
  }

  /**
   * @return -this (placed into a new matrix)
   */
  public Matrix neg() {
    return neg(create(rows(),cols()));
  }

  /**
   * @return -this (placed into "this")
   */
  public Matrix negL() {
    return neg(this);
  }

  //----------------------------------------------------------------------------

  /**
   * @return this + m (placed into "result")
   */
  public Matrix add(Matrix m, Matrix result) {
    assert (rows() == m.rows());
    assert (cols() == m.cols());
    assert (result.rows() == rows());
    assert (result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        result.set(i, j, get(i,j) + m.get(i,j));
      }
    }
    return result;
  }

  /**
   * @return this + m (placed into a new matrix)
   */
  public Matrix add(Matrix m) {
    return add(m, create(rows(),cols()));
  }

  /**
   * @return this + m (placed into "this")
   */
  public Matrix addL(Matrix m) {
    return add(m, this);
  }

  /**
   * @return this + m (placed into "m")
   */
  public Matrix addR(Matrix m) {
    return add(m, m);
  }

  //----------------------------------------------------------------------------

  /**
   * @return this - m (placed into "result")
   */
  public Matrix sub(Matrix m, Matrix result) {
    assert (rows() == m.rows());
    assert (cols() == m.cols());
    assert (result.rows() == rows());
    assert (result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        result.set(i, j, get(i,j) - m.get(i,j));
      }
    }
    return result;
  }

  /**
   * @return this - m (placed into a new matrix)
   */
  public Matrix sub(Matrix m) {
    return sub(m, create(rows(),cols()));
  }

  /**
   * @return this - m (placed into "this")
   */
  public Matrix subL(Matrix m) {
    return sub(m, this);
  }

  /**
   * @return this - m (placed into "m")
   */
  public Matrix subR(Matrix m) {
    return sub(m, m);
  }

  //----------------------------------------------------------------------------

  /**
   * @return this * c (placed into "result")
   */
  public Matrix mul(double c, Matrix result) {
    assert (result.rows() == rows());
    assert (result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        result.set(i, j, c * get(i,j));
      }
    }
    return result;
  }

  /**
   * @return this * c (placed into a new matrix)
   */
  public Matrix mul(double c) {
    return mul(c, create(rows(),cols()));
  }

  /**
   * @return this * c (placed into "this")
   */
  public Matrix mulL(double c) {
    return mul(c, this);
  }

  //----------------------------------------------------------------------------

  /**
   * @return this / c (placed into "result")
   */
  public Matrix div(double c, Matrix result) {
    return mul(1.0 / c, result);
  }

  /**
   * @return this / c (placed into a new matrix)
   */
  public Matrix div(double c) {
    return div(c, create(rows(),cols()));
  }

  /**
   * @return this / c (placed into "this")
   */
  public Matrix divL(double c) {
    return div(c, this);
  }

  //----------------------------------------------------------------------------

  /**
   * @return this .% m (placed into "result")
   */
  public Matrix mod(double m, Matrix result) {
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        result.set(i, j, get(i,j) % m);
      }
    }
    return result;
  }

  /**
   * @return this .% m (placed into a new matrix)
   */
  public Matrix mod(double m) {
    return mod(m, create(rows(),cols()));
  }

  /**
   * @return this .% m (placed into "this")
   */
  public Matrix modL(double m) {
    return mod(m, this);
  }

  //----------------------------------------------------------------------------

  /**
   * Matrix-vector product.
   * The "result" parameter has to be different from "v".
   * @return this * v (placed into "result")
   */
  public Vector mul(Vector v, Vector result) {
    assert (cols() == v.length());
    assert (result != v);
    assert (result.length() == rows());
    int i, j;
    double vj = v.get(0); // j = 0
    for (i = 0; i < rows(); ++i) { // initialize "result"
      result.set(i, vj * get(i,0));
    }
    for (j = 1; j < cols(); ++j) {
      vj = v.get(j);
      for (i = 0; i < rows(); ++i) {
        result.set(i, result.get(i) + vj * get(i,j));
      }
    }
    return result;
  }

  /**
   * Matrix-vector product.
   * @return this * v (placed into a new vector)
   */
  public Vector mul(Vector v) {
    return mul(v, Vector.create(rows()));
  }

  /**
   * Matrix product with a diagonal matrix represented by "v".
   * @return this * diag(v) (placed into "result")
   */
  public Matrix mulD(Vector v, Matrix result) {
    assert (cols() == v.length());
    assert (result.rows() == rows());
    assert (result.cols() == cols());
    int i, j;
    double d;
    for (j = 0; j < cols(); ++j) {
      d = v.get(j);
      for (i = 0; i < rows(); ++i)
        result.set(i, j, get(i,j) * d);
    }
    return result;
  }

  /**
   * Matrix product with a diagonal matrix represented by "v".
   * @return this * diag(v) (placed into a new matrix)
   */
  public Matrix mulD(Vector v) {
    return mulD(v, create(rows(),cols()));
  }

  /**
   * Matrix product.
   * The "result" parameter has to be different from "this" and "m".
   * @return this * m (placed into "result")
   */
  public Matrix mul(Matrix m, Matrix result) {
    assert (cols() == m.rows());
    assert (result != this);
    assert (result != m);
    assert (result.rows() == rows());
    assert (result.cols() == m.cols());
    int i, j, k;
    final int cols = m.cols();
    if (cols <= rows()) {
      double tik;
      for (i = 0; i < rows(); ++i) {
        tik = get(i,0); // k = 0
        for (j = 0; j < cols; ++j) { // initialize "result[i:*]"
          result.set(i, j, tik * m.get(0,j));
        }
        for (k = 1; k < cols(); ++k) {
          tik = get(i,k);
          for (j = 0; j < cols; ++j) {
            result.set(i, j, result.get(i,j) + tik * m.get(k,j));
          }
        }
      }
    }
    else { m.T().mul(T(), result.T()); }
    return result;
  }

  /**
   * Matrix product.
   * @return this * m (placed into a new matrix)
   */
  public Matrix mul(Matrix m) {
    return mul(m, create(rows(),m.cols()));
  }

  //----------------------------------------------------------------------------

  /**
   * Entrywise multiplication.
   * @return this .* m (placed into "result")
   */
  public Matrix emul(Matrix m, Matrix result) {
    assert (rows() == m.rows());
    assert (cols() == m.cols());
    assert (result.rows() == rows());
    assert (result.cols() == cols());
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        result.set(i, j, get(i,j) * m.get(i,j));
      }
    }
    return result;
  }

  /**
   * Entrywise multiplication.
   * @return this .* m (placed into a new matrix)
   */
  public Matrix emul(Matrix m) {
    return emul(m, create(rows(),cols()));
  }

  /**
   * Entrywise multiplication.
   * @return this .* m (placed into "this")
   */
  public Matrix emulL(Matrix m) {
    return emul(m, this);
  }

  /**
   * Entrywise multiplication.
   * @return this .* m (placed into "m")
   */
  public Matrix emulR(Matrix m) {
    return emul(m, m);
  }

  //----------------------------------------------------------------------------

  /**
   * @return 1-norm, maximum absolute column sum
   */
  public double norm1() {
    int i, j;
    double s, maxs = 0.0;
    for (j = 0; j < cols(); ++j) {
      s = 0.0;
      for (i = 0; i < rows(); ++i) { s += Math.abs(get(i,j)); }
      if (s > maxs) { maxs = s; }
    }
    return maxs;
  }

  /**
   * @return inf-norm, maximum absolute row sum
   */
  public double normi() {
    int i, j;
    double s, maxs = 0.0;
    for (i = 0; i < rows(); ++i) {
      s = 0.0;
      for (j = 0; j < cols(); ++j) { s += Math.abs(get(i,j)); }
      if (s > maxs) { maxs = s; }
    }
    return maxs;
  }

  /**
   * @return Frobenius norm
   */
  public double normf() {
    int i, j;
    double s = 0.0, v;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        v = get(i,j);
        s += v * v;
      }
    }
    return Math.sqrt(s);
  }

  //----------------------------------------------------------------------------

  /**
   * @return trace of a square matrix
   */
  public double tr() {
    assert (rows() == cols());
        
    double trace = 0.0;
    for (int i = 0; i < rows(); ++i) { trace += get(i,i); }
    return trace;
  }

  /**
   * @return trace of the square ("this" * A) matrix
   */
  public double trMul(Matrix A) {
    assert (cols() == A.rows());
    assert (rows() == A.cols());
        
    int i, j;
    double trace = 0.0;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        trace += get(i,j) * A.get(j,i);
      }
    }
    return trace;
  }

  //----------------------------------------------------------------------------

  /**
   * Cholesky decomposition of a (symmetric) positive-definite matrix.
   * The provided "result" matrix is assumed to be zero on its
   * upper-triangular half.
   * The "result" parameter also has to be different from "this".
   * @return Cholesky lower-triangular matrix (placed into "result")
   */
  public Matrix choleskyL(Matrix result) {
    assert (rows() == cols());
    assert (result.rows() == rows());
    assert (result.rows() == result.cols());
    final int n = rows();
    int i, j, k;
    double Ljj, Lij, v;
    for (j = 0; j < n; ++j) {
      Ljj = get(j,j);
      for (k = 0; k < j; ++k) {
        v = result.get(j,k);
        Ljj -= v * v;
      }
      Ljj = Math.sqrt(Ljj);
      result.set(j, j, Ljj);
            
      for (i = j+1; i < n; ++i) {
        Lij = get(i,j);
        for (k = 0; k < j; ++k) {
          Lij -= result.get(i,k) * result.get(j,k);
        }
        result.set(i, j, Lij / Ljj);
      }
    }
    return result;
  }

  /**
   * Cholesky decomposition of a (symmetric) positive-definite matrix.
   * @return Cholesky lower-triangular matrix (placed into a new matrix)
   */
  public Matrix choleskyL() {
    return choleskyL(zero(rows(), cols()));
  }

  /**
   * Cholesky decomposition of a (symmetric) positive-definite matrix.
   * Provide the result in matrix "L" and vector "D" for which 
   * this = L * diag(D) * L.T holds. Here the diagonal elements of L are all 
   * one. The provided L matrix has to be zero above the diagonal.
   */
  public void choleskyLD(Matrix L, Vector D) {
    assert (rows() == cols());
    assert (L.rows() == rows());
    assert (L.rows() == L.cols());
    assert (D.length() == rows());
    final int n = rows();
    int i, j, k;
    double Dj, Lij, v;
    for (j = 0; j < n; ++j) {
      Dj = get(j,j);
      for (k = 0; k < j; ++k) {
        v = L.get(j,k);
        Dj -= v * v * D.get(k);
      }
      D.set(j, Dj);
      L.set(j, j, 1.0);
            
      for (i = j+1; i < n; ++i) {
        Lij = get(i,j);
        for (k = 0; k < j; ++k) {
          Lij -= L.get(i,k) * L.get(j,k) * D.get(k);
        }
        L.set(i, j, Lij / Dj);
      }
    }
  }

  //----------------------------------------------------------------------------

  /**
   * QR decomposition of an arbitrary matrix using Hauseholder transformations.
   * @return {Q R} (placed into new matrices)
   */
  public Matrix[] QR() {
    Matrix Q = Matrix.create(rows(), rows());
    Matrix R = Matrix.create(rows(), cols());
    return QR(new Matrix[]{Q, R}, Vector.create(rows()));
  }

  /**
   * QR decomposition of an arbitrary matrix using Hauseholder transformations.
   * The "tmpV" vector is a temporary storage for the computation
   * for which tmpV.length >= min(rows-1,cols) has to hold.
   * @return {Q R} (placed into "result")
   */
  public Matrix[] QR(Matrix result[], Vector tmpV) {
    QR(result[0], result[1], true, tmpV);
    return result;
  }

  /**
   * QR decomposition of an arbitrary matrix using Hauseholder transformations.
   * Only the R matrix is computed.
   * The "tmpV" vector is a temporary storage for the computation
   * for which tmpV.length >= min(rows-1,cols) has to hold.
   * @return R (placed into "R")
   */
  public Matrix QR(Matrix R, Vector tmpV) {
    QR(null, R, false, tmpV);
    return R;
  }

  /**
   * QR decomposition of an arbitrary matrix using Hauseholder transformations.
   * Both the Q and R matrices are computed into "Q" and "R", respectively.
   * The "tmpV" vector is a temporary storage for the computation
   * for which tmpV.length >= min(rows-1,cols) has to hold.
   */
  public void QR(Matrix Q, Matrix R, Vector tmpV) {
    QR(Q, R, true, tmpV);
  }

  private void QR(Matrix Q, Matrix R, boolean isComputeQ, Vector tmpV) {
    assert (Q.rows() == rows() && Q.cols() == rows());
    assert (R.rows() == rows() && R.cols() == cols());

    final int rows = rows(), cols = cols();
    final int t = Math.min(rows-1, cols);
    assert (t <= tmpV.length());

    int i, j, k;
    double norm, s;

    copy(R);
    for (k = 0; k < t; ++k) {
      norm = 0.0;
      for (i = k; i < rows; ++i) { norm = hypot(norm, R.get(i, k)); }

      if (norm != 0.0) {
        if (R.get(k, k) < 0) norm = -norm;

        for (i = k; i < rows; ++i) { R.set(i, k, R.get(i, k) / norm); }
        R.set(k, k, R.get(k, k) + 1.0); // 1 <= R[k,k] <= 2

        for (j = k+1; j < cols; ++j) {
          s = 0.0;
          for (i = k; i < rows; ++i) {
            s += R.get(i, k) * R.get(i, j);
          }
          s = -s / R.get(k, k);
          for (i = k; i < rows; ++i) {
            R.set(i, j, R.get(i, j) + s*R.get(i, k));
          }
        }
      }
      tmpV.set(k, -norm);
    }

    if (isComputeQ) {
      Q.setToEye();
      for (k = t-1; k >= 0; --k) {
        for (j = k; j < rows; ++j) {
          if (R.get(k, k) != 0.0) {
            s = 0.0;
            for (i = k; i < rows; ++i) {
              s += R.get(i, k) * Q.get(i, j);
            }
            s = -s / R.get(k, k);
            for (i = k; i < rows; ++i) {
              Q.set(i, j, Q.get(i, j) + s*R.get(i, k));
            }
          }
        }
      }
    }

    for (k = 0; k < t; ++k) {
      R.set(k, k, tmpV.get(k));
      for (i = k+1; i < rows; ++i) { R.set(i, k, 0.0); }
    }
  }

  private double hypot(double x, double y) {
    final double absX = Math.abs(x), absY = Math.abs(y);

    double r;
    if (absX > absY) {
      r = y/x;
      r = absX * Math.sqrt(1.0 + r*r);
    }
    else if (y != 0.0) {
      r = x/y;
      r = absY * Math.sqrt(1.0 + r*r);
    }
    else { r = 0.0; }
    return r;
  }

  //----------------------------------------------------------------------------

  /**
   * @return determinant of a lower/upper triangular matrix
   */
  public double detT() {
    if (rows() != cols()) return 0.0;
        
    double det = 1.0;
    for (int i = 0; i < rows(); ++i) { det *= get(i,i); }
    return det;
  }

  /**
   * @return determinant of a 2x2 matrix
   */
  public double det2x2() {
    assert (rows() == 2 && cols() == 2);
    return get(0,0)*get(1,1) - get(0,1)*get(1,0);
  }

  /**
   * @return determinant of a 3x3 matrix
   */
  public double det3x3() {
    assert (rows() == 3 && cols() == 3);
    return get(0,0) * (get(1,1)*get(2,2) - get(1,2)*get(2,1))
         + get(0,1) * (get(1,2)*get(2,0) - get(1,0)*get(2,2))
         + get(0,2) * (get(1,0)*get(2,1) - get(1,1)*get(2,0));
  }

  //----------------------------------------------------------------------------

  /**
   * @return inverse of a 2x2 matrix (placed into "result")
   *         or null if the determinant is zero ("result" remains unchanged)
   */
  public Matrix inv2x2(Matrix result) {
    assert (rows() == 2 && cols() == 2);
    double a = get(0,0), b = get(0,1),
           c = get(1,0), d = get(1,1);
    double det = a*d - b*c;
    if (0.0 == det) return null;
    if (null == result) result = create(2,2);
    result.set(0, 0,  d); result.set(0, 1, -b);
    result.set(1, 0, -c); result.set(1, 1,  a);
    return result.divL(det);
  }

  /**
   * @return inverse of a 2x2 matrix (placed into a new matrix)
   *         or null if the determinant is zero
   */
  public Matrix inv2x2() {
    return inv2x2(null);
  }

  /**
   * @return inverse of a 2x2 matrix (placed into "this")
   *         or null if the determinant is zero ("this" remains unchanged)
   */
  public Matrix inv2x2L() {
    return inv2x2(this);
  }

  /**
   * @return inverse of a 3x3 matrix (placed into "result")
   *         or null if the determinant is zero ("result" remains unchanged)
   */
  public Matrix inv3x3(Matrix result) {
    assert (rows() == 3 && cols() == 3);
    double a = get(0,0), b = get(0,1), c = get(0,2),
           d = get(1,0), e = get(1,1), f = get(1,2),
           g = get(2,0), h = get(2,1), k = get(2,2);
    double A = e*k-f*h, D = c*h-b*k, G = b*f-c*e,
           B = f*g-d*k, E = a*k-c*g, H = c*d-a*f,
           C = d*h-e*g, F = g*b-a*h, K = a*e-b*d;
    double det = a*A + b*B + c*C;
    if (0.0 == det) return null;
    if (null == result) result = create(3,3);
    result.set(0, 0, A); result.set(0, 1, D); result.set(0, 2, G);
    result.set(1, 0, B); result.set(1, 1, E); result.set(1, 2, H);
    result.set(2, 0, C); result.set(2, 1, F); result.set(2, 2, K);
    return result.divL(det);
  }

  /**
   * @return inverse of a 3x3 matrix (placed into a new matrix)
   *         or null if the determinant is zero
   */    
  public Matrix inv3x3() {
    return inv3x3(null);
  }

  /**
   * @return inverse of a 3x3 matrix (placed into "this")
   *         or null if the determinant is zero ("this" remains unchanged)
   */
  public Matrix inv3x3L() {
    return inv3x3(this);
  }

  /**
   * Invert a diagonal matrix having non-zero diagonal elements.
   * The provided "result" matrix is assumed to be zero on its
   * non-diagonal elements.
   * @return diagonal matrix inverse (placed into "result")
   */
  public Matrix invD(Matrix result) {
    for (int i = 0; i < rows(); ++i) {
      result.set(i, i, 1.0 / get(i,i));
    }
    return result;
  }

  /**
   * Invert a diagonal matrix having non-zero diagonal elements.
   * @return diagonal matrix inverse (placed into a new matrix)
   */
  public Matrix invD() {
    return invD(zero(rows(), cols()));
  }

  /**
   * Invert a lower triangular matrix having non-zero diagonal elements.
   * The provided "result" matrix is assumed to be zero on its
   * upper-triangular half.
   * The "result" parameter also has to be different from "this".
   * @return lower-triangular matrix inverse (placed into "result")
   */
  public Matrix invLT(Matrix result) {
    final int n = rows();
    int i, j, k;
    double Rii, Rij;
    for (i = 0; i < n; ++i) { result.set(i, i, 1.0 / get(i,i)); }
    for (i = 0; i < n; ++i) {
      Rii = result.get(i,i);
      for (j = 0; j < i; ++j) {
        Rij = 0.0;
        for (k = 0; k < i; ++k) { Rij -= get(i,k) * result.get(k,j); }
        result.set(i, j, Rii * Rij);
      }
    }
    return result;
  }

  /**
   * Invert a lower triangular matrix.
   * @return lower-triangular matrix inverse (placed into a new matrix)
   */
  public Matrix invLT() {
    return invLT(zero(rows(), cols()));
  }

  /**
   * Invert a (symmetric) positive-definite matrix.
   * This version is slightly faster, but not as stable numerically as the 
   * invPD(Matrix,Matrix,Vector,Matrix) one.
   * The provided "invL" matrix is assumed to be zero on its
   * upper-triangular half. They also have to be different from "this" and
   * each other. Also, "result" = invL^T * invL will hold.
   * @return positive-definite matrix inverse (placed into "result")
   */
  public Matrix invPD(Matrix result, Matrix invL) {
    assert (result != this && result != invL && invL != this);
    choleskyL(result);
    result.invLT(invL);
    return invL.T().mul(invL, result);
  }

  /**
   * Invert a (symmetric) positive-definite matrix.
   * The provided "invL" matrix is assumed to be zero on its upper-triangular
   * half. All matrix parameters have to be different from "this" and each 
   * other. Also, "result" = invL^T * diag(invD) * invL and
   * "tmp" = invL^T * diag(invD) will hold.
   * @return positive-definite matrix inverse (placed into "result")
   */
  public Matrix invPD(Matrix result, Matrix invL, Vector invD, Matrix tmp) {
    assert (result != this && result != invL && invL != this);
    choleskyLD(result, invD);
    result.invLT(invL);
    invD.reciprocL();
    return invL.T().mulD(invD, tmp).mul(invL, result);
  }

  /**
   * Invert a (symmetric) positive-definite matrix.
   * @return positive-definite matrix inverse (placed into a new matrix)
   */
  public Matrix invPD() {
    return invPD(zero(rows(), cols()),
                 zero(rows(), cols()),
                 Vector.zero(rows()),
                 zero(rows(), cols()));
  }

  //----------------------------------------------------------------------------

  /**
   * @return transpose of the matrix
   */
  abstract public Matrix T();

  //----------------------------------------------------------------------------

  @Override
  public String toString() {
    String str = "[";
    int i, j;
    for (i = 0; i < rows(); ++i) {
      if (0 != i) str += "; ";
      for (j = 0; j < cols(); ++j) {
        if (0 != j) str += " ";
        str += get(i,j);
      }
    }
    str += "]";
    return str;
  }

  //----------------------------------------------------------------------------

  protected final int rows, cols;
  private final double[][] data;
}
