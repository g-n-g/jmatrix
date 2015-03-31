package jvecmat;

import java.util.Random;

/**
 * Abstract matrix representation.
 */
public abstract class Matrix implements VecMat {

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

  @Override
  public boolean hasNaN() {
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        if (Double.isNaN(get(i,j))) { return true; }
      }
    }
    return false;
  }

  @Override
  public boolean hasInf() {
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        if (Double.isInfinite(get(i,j))) { return true; }
      }
    }
    return false;
  }

  @Override
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
   * Returns a copy of the vector placed into <code>result</code>.
   * If <code>result</code> is <code>null</code>, a new object is created.
   *
   * @param result appropriately sized storage for the copy
   *               or <code>null</code>
   * @return copy of the matrix
   */
  public Matrix copy(Matrix result) {
    if (result == null) { result = Matrix.create(rows(), cols()); }
    assert(result.rows() == rows());
    assert(result.cols() == cols());
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        result.set(i, j, get(i,j));
      }
    }
    return result;
  }

  @Override
  public Matrix copy() {
    return copy(null);
  }

  //----------------------------------------------------------------------------

  /**
   * Returns the number of rows of the matrix.
   *
   * @return number of rows of the matrix
   */
  public final int rows() {
    return rows;
  }

  /**
   * Returns the number of columns of the matrix.
   *
   * @return number of columns of the matrix
   */
  public final int cols() {
    return cols;
  }

  /**
   * Returns the array representation of the matrix.
   * The data is not copied, so any changes to the returned array
   * will change the matrix too.
   *
   * @return the array representation of the matrix
   */
  public final double[][] array() {
    return data;
  }

  /**
   * Returns the element of the matrix at row <code>i</code>
   * and column <code>j</code>.
   *
   * @param i the row index of the matrix element
   * @param j the column index of the matrix element
   * @return the (i,j) element of the matrix
   */
  abstract public double get(int i, int j);

  /**
   * Set the matrix element at row <code>i</code> and column <code>j</code>
   * to <code>value</code>.
   *
   * @param i the row index of the matrix element
   * @param j the column index of the matrix element
   * @param value the new value
   */
  abstract public void set(int i, int j, double value);

  //----------------------------------------------------------------------------

  @Override
  public Matrix setToConstant(double c) {
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        set(i, j, c);
      }
    }
    return this;
  }

  @Override
  public Matrix setToZero() {
    return setToConstant(0.0);
  }

  @Override
  public Matrix setToOne() {
    return setToConstant(1.0);
  }

  /**
   * Set the diagonal elements to one and all other elements to zero.
   * If the matrix is square, this sets the identity matrix.
   *
   * @return <code>this</code> matrix
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

  @Override
  public Matrix setToRand(Random rng) {
    int i, j;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        set(i, j, rng.nextDouble());
      }
    }
    return this;
  }

  @Override
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
   * Returns the diagonal elements as a vector.
   * A new vector is created if <code>result</code> is <code>null</code>.
   * Otherwise, the length of <code>result</code> has to match the number
   * of diagonal elements (the minimum of the number of rows and columns).
   *
   * @param result storage of the result or <code>null</code>
   * @return vector of the diagonal elements
   */
  public Vector getDiag(Vector result) {
    int n = Math.min(rows(), cols());
    if (result == null) { result = Vector.create(n); }
    assert(result.length() == n);
    for (int i = 0; i < n; ++i) { result.set(i, get(i,i)); }
    return result;
  }

  /**
   * Returns the diagonal elements in a new vector.
   *
   * @return vector of the diagonal elements
   */
  public Vector getDiag() {
    return getDiag(null);
  }

  /**
   * Returns the <code>i</code>-th row as a vector.
   * A new vector is created if <code>result</code> is <code>null</code>.
   * Otherwise, the length of <code>result</code> has to match the number
   * of columns.
   *
   * @param i row index
   * @param result storage of the result or <code>null</code>
   * @return vector of the <code>i</code>-th row
   */
  public Vector getRow(int i, Vector result) {
    assert (0 <= i && i <= rows());
    if (result == null) { result = Vector.create(cols()); }
    assert (result.length() == cols());
    for (int j = 0; j < cols(); ++j) { result.set(j, get(i,j)); }
    return result;
  }

  /**
   * Returns the <code>i</code>-th row in a new vector.
   *
   * @param i row index
   * @return vector of the <code>i</code>-th row
   */
  public Vector getRow(int i) {
    return getRow(i, null);
  }

  /**
   * Sets the <code>i</code>-th row to <code>v</code>.
   *
   * @param i row index
   * @param v new row elements
   * @return <code>this</code> matrix
   */
  public Matrix setRow(int i, Vector v) {
    assert (0 <= i && i <= rows());
    assert (v != null && cols() == v.length());
    for (int j = 0; j < v.length(); ++j) { set(i, j, v.get(j)); }
    return this;
  }

  /**
   * Sets the <code>i</code>-th row to the <code>i</code>-th row
   * of matrix <code>m</code>.
   *
   * @param i row index
   * @param m matrix containing the new row elements
   * @return <code>this</code> matrix
   */
  public Matrix setRow(int i, Matrix m) {
    assert (0 <= i && i <= rows());
    assert (m != null && cols() == m.cols());
    for (int j = 0; j < m.cols(); ++j) { set(i, j, m.get(i, j)); }
    return this;
  }

  /**
   * Returns the <code>j</code>-th column as a vector.
   * A new vector is created if <code>result</code> is <code>null</code>.
   * Otherwise, the length of <code>result</code> has to match the number
   * of rows.
   *
   * @param j column index
   * @param result storage of the result or <code>null</code>
   * @return vector of the <code>j</code>-th column
   */
  public Vector getCol(int j, Vector result) {
    assert (0 <= j && j <= cols());
    if (result == null) { result = Vector.create(rows()); }
    assert (result.length() == rows());
    for (int i = 0; i < rows(); ++i) { result.set(i, get(i,j)); }
    return result;
  }

  /**
   * Returns the <code>j</code>-th column in a new vector.
   *
   * @param j column index
   * @return vector of the <code>j</code>-th column
   */
  public Vector getCol(int j) {
    return getCol(j, null);
  }

  /**
   * Sets the <code>j</code>-th column to <code>v</code>.
   *
   * @param j column index
   * @param v new column elements
   * @return <code>this</code> matrix
   */
  public Matrix setCol(int j, Vector v) {
    assert (0 <= j && j <= cols());
    assert (v != null && rows() == v.length());
    for (int i = 0; i < v.length(); ++i) { set(i, j, v.get(i)); }
    return this;
  }

  /**
   * Sets the <code>j</code>-th column to the <code>j</code>-th row
   * of matrix <code>m</code>.
   *
   * @param j column index
   * @param m matrix containing the new column elements
   * @return <code>this</code> matrix
   */
  public Matrix setCol(int j, Matrix m) {
    assert (0 <= j && j <= cols());
    assert (m != null && rows() == m.rows());
    for (int i = 0; i < m.rows(); ++i) { set(i, j, m.get(i, j)); }
    return this;
  }

  /**
   * Returns the submatrix having rows from <code>iF</code> to <code>iT</code>
   * and columns from <code>jF</code> to <code>jT</code> (all inclusive).
   * A new matrix is created if <code>result</code> is <code>null</code>.
   * Otherwise, the size of <code>result</code> has to match the submatrix.
   *
   * @param iF index of the first row
   * @param iT index of the last row
   * @param jF index of the first column
   * @param jT index of the last column
   * @param result storage of the result or <code>null</code>
   * @return submatrix formed by rows from <code>iF</code> to <code>iT</code>
   *                   and columns from <code>jF</code> to <code>jT</code>
   */
  public Matrix getMat(int iF, int iT, int jF, int jT, Matrix result) {
    assert (0 <= iF && iF <= iT && iT < rows());
    assert (0 <= jF && jF <= jT && jT < cols());
    if (result == null) { result = create(iT-iF+1, jT-jF+1); }
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
   * Returns the submatrix having rows from <code>iF</code> to <code>iT</code>
   * and columns from <code>jF</code> to <code>jT</code> (all inclusive).
   *
   * @param iF index of the first row
   * @param iT index of the last row
   * @param jF index of the first column
   * @param jT index of the last column
   * @return submatrix formed by rows from <code>iF</code> to <code>iT</code>
   *                   and columns from <code>jF</code> to <code>jT</code>
   */
  public Matrix getMat(int iF, int iT, int jF, int jT) {
    return getMat(iF, iT, jF, jT, null);
  }

  /**
   * Sets the region specified by rows from <code>iF</code> to <code>iT</code>
   * and columns from <code>jF</code> to <code>jT</code> (all inclusive)
   * using matrix <code>m</code>.
   * Matrix <code>m</code> has to be the same size as the specified region.
   *
   * @param iF index of the first row
   * @param iT index of the last row
   * @param jF index of the first column
   * @param jT index of the last column
   * @param m matrix containing the new elements
   * @return <code>this</code> matrix
   */
  public Matrix setMat(int iF, int iT, int jF, int jT, Matrix m) {
    assert (0 <= iF && iF <= iT && iT < rows());
    assert (0 <= jF && jF <= jT && jT < cols());
    assert (m != null && m.rows() == iT-iF+1 && m.cols() == jT-jF+1);
    for (int i = 0, ii = iF; ii <= iT; ++i, ++ii) {
      for (int j = 0, jj = jF; jj <= jT; ++j, ++jj) {
        set(ii, jj, m.get(i, j));
      }
    }
    return this;
  }

  //----------------------------------------------------------------------------

  /**
   * Entrywise absolute value operation placed into <code>result</code>.
   *
   * @param result storage of the result (should not be <code>null</code>)
   * @return entrywise absolute value
   */
  public Matrix abs(Matrix result) {
    assert (result != null);
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

  @Override
  public Matrix abs() {
    return abs(create(rows(), cols()));
  }

  @Override
  public Matrix absL() {
    return abs(this);
  }

  //----------------------------------------------------------------------------

  /**
   * Entrywise sign operation with zero replacement
   * placed into <code>result</code>.
   *
   * @param zeroReplacement value replacing 0.0 values
   * @param result storage of the result (should not be <code>null</code>)
   * @return entrywise sign value
   */
  public Matrix sign(double zeroReplacement, Matrix result) {
    assert (result != null);
    assert (result.rows() == rows());
    assert (result.cols() == cols());
    int i, j;
    double value;
    for (i = 0; i < rows(); ++i) {
      for (j = 0; j < cols(); ++j) {
        value = get(i,j);
        result.set(i, j, (0.0 == value) ? zeroReplacement : Math.signum(value));
      }
    }
    return result;
  }

  /**
   * Entrywise sign operation placed into <code>result</code>.
   *
   * @param result storage of the result (should not be <code>null</code>)
   * @return entrywise sign value
   */
  public Matrix sign(Matrix result) {
    assert (result != null);
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

  @Override
  public Matrix sign(double zeroReplacement) {
    return sign(zeroReplacement, create(rows(), cols()));
  }

  @Override
  public Matrix sign() {
    return sign(create(rows(), cols()));
  }

  @Override
  public Matrix signL(double zeroReplacement) {
    return sign(zeroReplacement, this);
  }

  @Override
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
    Matrix Q = create(rows(), rows());
    Matrix R = create(rows(), cols());
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
