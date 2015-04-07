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
    return new DenseMatrix(data);
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
   * Entrywise absolute value (in <code>result</code>).
   *
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param result storage of the result (not <code>null</code>)
   * @return entrywise absolute value
   */
  public Matrix abs(Matrix result) {
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        result.set(i, j, Math.abs(get(i,j)));
      }
    }
    return result;
  }

  @Override
  public Matrix abs() {
    return abs(create(rows(), cols()));
  }

  //----------------------------------------------------------------------------

  /**
   * Entrywise sign operation with zero replacement (in <code>result</code>).
   *
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param zeroReplacement value replacing 0.0 values
   * @param result storage of the result (not <code>null</code>)
   * @return entrywise sign value
   */
  public Matrix sign(double zeroReplacement, Matrix result) {
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        double value = get(i,j);
        result.set(i, j, (0.0 == value) ? zeroReplacement : Math.signum(value));
      }
    }
    return result;
  }

  /**
   * Entrywise sign operation (in <code>result</code>).
   *
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param result storage of the result (not <code>null</code>)
   * @return entrywise sign value
   */
  public Matrix sign(Matrix result) {
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
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

  //----------------------------------------------------------------------------

  /**
   * Entrywise negation (in <code>result</code>).
   *
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the negated matrix
   */
  public Matrix neg(Matrix result) {
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        result.set(i, j, -get(i,j));
      }
    }
    return result;
  }

  @Override
  public Matrix neg() {
    return neg(create(rows(), cols()));
  }

  //----------------------------------------------------------------------------

  /**
   * Matrix-matrix addition (in <code>result</code>).
   * Adds matrix <code>m</code> to <code>this</code> matrix.
   *
   * Matrices <code>m</code> and <code>result</code> have to have the same size
   * as <code>this</code>.
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * or <code>m</code> providing in-place operation.
   *
   * @param m matrix to add (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the sum of
   *         <code>this</code> and <code>m</code>
   */
  public Matrix add(Matrix m, Matrix result) {
    assert (m != null && rows() == m.rows() && cols() == m.cols());
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        result.set(i, j, get(i,j) + m.get(i,j));
      }
    }
    return result;
  }

  /**
   * Matrix-matrix addition (in new matrix).
   * Adds matrix <code>m</code> to <code>this</code> matrix in a new matrix.
   *
   * Matrix <code>m</code> has to have the same size as <code>this</code>.
   *
   * @param m matrix to add (not <code>null</code>)
   * @return the sum of <code>this</code> and <code>m</code> in a new matrix
   */
  public Matrix add(Matrix m) {
    return add(m, create(rows(), cols()));
  }

  // TODO: add constant
  // TODO: addRow(Vector)
  // TODO: addCol(Vector)

  //----------------------------------------------------------------------------

  /**
   * Matrix-matrix subtraction (in <code>result</code>).
   * Subtracts matrix <code>m</code> from <code>this</code> matrix.
   *
   * Matrices <code>m</code> and <code>result</code> have to have the same size
   * as <code>this</code>.
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * or <code>m</code> providing in-place operation.
   *
   * @param m matrix to subtract (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the difference of
   *         <code>this</code> and <code>m</code>
   */
  public Matrix sub(Matrix m, Matrix result) {
    assert (m != null && rows() == m.rows() && cols() == m.cols());
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        result.set(i, j, get(i,j) - m.get(i,j));
      }
    }
    return result;
  }

  /**
   * Matrix-matrix subtraction (in new matrix).
   * Subtracts matrix <code>m</code> from <code>this</code> matrix
   * in a new matrix.
   *
   * Matrix <code>m</code> has to have the same size as <code>this</code>.
   *
   * @param m matrix to subtract (not <code>null</code>)
   * @return the difference of <code>this</code> and <code>m</code>
   *         in a new matrix
   */
  public Matrix sub(Matrix m) {
    return sub(m, create(rows(),cols()));
  }

  // TODO: sub constant
  // TODO: subRow(Vector)
  // TODO: subCol(Vector)

  //----------------------------------------------------------------------------

  /**
   * Matrix-constant division (in <code>result</code>). Divides all elements of
   * <code>this</code> matrix by constant <code>c</code>.
   *
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param c constant to divide with
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the values of
   *         <code>this</code> divided by <code>c</code>
   */
  public Matrix div(double c, Matrix result) {
    return mul(1.0 / c, result);
  }

  @Override
  public Matrix div(double c) {
    return div(c, create(rows(),cols()));
  }

  // TODO: divRow(Vector)
  // TODO: divCol(Vector)

  //----------------------------------------------------------------------------

  /**
   * Matrix-constant multiplication (in <code>result</code>). Multiplies all
   * elements of <code>this</code> matrix by constant <code>c</code>.
   *
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param c constant to multiply with
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the values of
   *         <code>this</code> multiplied by <code>c</code>
   */
  public Matrix mul(double c, Matrix result) {
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        result.set(i, j, c * get(i,j));
      }
    }
    return result;
  }

  @Override
  public Matrix mul(double c) {
    return mul(c, create(rows(), cols()));
  }

  // TODO: mulRow(Vector)
  // TODO: mulCol(Vector)

  /**
   * Matrix-vector multiplication (in <code>result</code>).
   *
   * The length of vector <code>v</code> has to be equal to the column number
   * of <code>this</code> matrix.
   * The length of <code>result</code> has to be equal to the row number
   * of <code>this</code> matrix.
   *
   * @param v vector to multiply with on the right
   *               (not <code>null</code> and not equal to <code>result</code>)
   * @param result storage of the result
   *               (not <code>null</code> and not equal to <code>v</code>)
   * @return <code>this</code> matrix multiplied
   *         by vector <code>v</code> from the right side
   */
  public Vector mul(Vector v, Vector result) {
    assert (v != null && v.length() == cols());
    assert (result != null && result != v && result.length() == rows());
    double vj = v.get(0); // j = 0
    for (int i = 0; i < rows(); ++i) { // initialize "result"
      result.set(i, vj * get(i,0));
    }
    for (int j = 1; j < cols(); ++j) {
      vj = v.get(j);
      for (int i = 0; i < rows(); ++i) {
        result.set(i, result.get(i) + vj * get(i,j));
      }
    }
    return result;
  }

  /**
   * Matrix-vector multiplication (in new vector).
   *
   * The length of vector <code>v</code> has to be equal to the column number
   * of <code>this</code> matrix.
   *
   * @param v vector to multiply with on the right
   * @return <code>this</code> matrix multiplied
   *         by vector <code>v</code> from the right side
   */
  public Vector mul(Vector v) {
    return mul(v, Vector.create(rows()));
  }

  /**
   * Matrix multiplication from the right with a diagonal matrix
   * represented by vector <code>v</code> (in <code>result</code>).
   *
   * The length of vector <code>v</code> has to be equal to the column number
   * of <code>this</code> matrix.
   * Furthermore, <code>result</code> has to have the same size
   * as <code>this</code> matrix. The <code>result</code> parameter can also
   * be set to <code>this</code> providing in-place operation.
   *
   * @param v vector representing the multiplier diagonal matrix
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having <code>this</code> matrix multiplied by
   *         <code>diag(v)</code> from the right
   */
  public Matrix mulD(Vector v, Matrix result) {
    assert (v != null && v.length() == cols());
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int j = 0; j < cols(); ++j) {
      double d = v.get(j);
      for (int i = 0; i < rows(); ++i)
        result.set(i, j, get(i,j) * d);
    }
    return result;
  }

  /**
   * Matrix multiplication from the right with a diagonal matrix
   * represented by vector <code>v</code> (in new matrix).
   *
   * The length of vector <code>v</code> has to be equal to the column number
   * of <code>this</code> matrix.
   *
   * @param v vector representing the multiplier diagonal matrix
   * @return <code>result</code> having <code>this</code> matrix multiplied by
   *         <code>diag(v)</code> from the right
   */
  public Matrix mulD(Vector v) {
    return mulD(v, create(rows(), cols()));
  }

  /**
   * Matrix-matrix multiplication (in <code>result</code>).
   *
   * The row number of matrix <code>m</code> has to match the column number
   * of <code>this</code> matrix.
   * The size of <code>result</code> has to match the row number of
   * <code>this</code> matrix and the column number of <code>m</code>.
   *
   * @param m matrix multiplier from the right
   *        (not <code>null</code> and not equal to <code>result</code>)
   * @param result storage of the result
   *        (not <code>null</code>
   *         and not equal to <code>this</code> or <code>m</code>)
   * @return <code>result</code> having <code>this</code> matrix multiplied
   *         by matrix <code>m</code> from the right
   */
  public Matrix mul(Matrix m, Matrix result) {
    assert (m != null && m.rows() == cols());
    assert (result != null && result != this && result != m);
    assert (result.rows() == rows() && result.cols() == m.cols());
    final int cols = m.cols();
    if (cols <= rows()) {
      double tik;
      for (int i = 0; i < rows(); ++i) {
        tik = get(i,0); // k = 0
        for (int j = 0; j < cols; ++j) { // initialize "result[i:*]"
          result.set(i, j, tik * m.get(0,j));
        }
        for (int k = 1; k < cols(); ++k) {
          tik = get(i,k);
          for (int j = 0; j < cols; ++j) {
            result.set(i, j, result.get(i,j) + tik * m.get(k,j));
          }
        }
      }
    }
    else {
      m.T().mul(T(), result.T());
    }
    return result;
  }

  /**
   * Matrix-matrix multiplication (in new matrix).
   *
   * The row number of matrix <code>m</code> has to match the column number
   * of <code>this</code> matrix.
   *
   * @param m matrix multiplier from the right
   * @return <code>result</code> having <code>this</code> matrix multiplied
   *         by matrix <code>m</code> from the right
   */
  public Matrix mul(Matrix m) {
    return mul(m, create(rows(), m.cols()));
  }

  //----------------------------------------------------------------------------

  /**
   * Takes the reciproc of all elements (in <code>result</code>).
   *
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> holding the elementwise reciproc matrix
   */
  public Matrix reciproc(Matrix result) {
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        result.set(i, j, 1.0 / get(i,j));
      }
    }
    return result;
  }

  @Override
  public Matrix reciproc() {
    return reciproc(create(rows(), cols()));
  }

  //----------------------------------------------------------------------------

  /**
   * Entrywise (signed) remainder with respect to modulus <code>m</code>
   * (in <code>result</code>).
   *
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param m modulus
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the remainders of the values of
   *         <code>this</code> with respect to modulus <code>m</code>
   */
  public Matrix mod(double m, Matrix result) {
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        result.set(i, j, get(i,j) % m);
      }
    }
    return result;
  }

  @Override
  public Matrix mod(double m) {
    return mod(m, create(rows(), cols()));
  }

  //----------------------------------------------------------------------------

  /**
   * Entrywise multiplication (Hadamard product) by matrix <code>m</code>
   * (in <code>result</code>).
   *
   * Matrices <code>v</code> and <code>result</code> have to have the same size
   * as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param m matrix to multiply with entrywise (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the elements of <code>this</code>
   *         and <code>m</code> multiplied entrywise
   */
  public Matrix emul(Matrix m, Matrix result) {
    assert (m != null && rows() == m.rows() && cols() == m.cols());
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        result.set(i, j, get(i,j) * m.get(i,j));
      }
    }
    return result;
  }

  /**
   * Entrywise multiplication (Hadamard product) by matrix <code>m</code>
   * (in new matrix).
   *
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   *
   * @param m matrix to multiply with entrywise (not <code>null</code>)
   * @return new matrix having the elements of <code>this</code>
   *         and <code>m</code> multiplied entrywise
   */
  public Matrix emul(Matrix m) {
    return emul(m, create(rows(), cols()));
  }

  //----------------------------------------------------------------------------

  /**
   * Computes the 1-norm (maximum absolute column sum).
   *
   * @return 1-norm of <code>this</code> matrix
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
   * Computes the infinity-norm (maximum absolute row sum).
   *
   * @return infinity-norm of <code>this</code> matrix
   */
  public double normI() {
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
   * Computes the Frobenius norm (maximum absolute row sum).
   *
   * @return Frobenius norm of <code>this</code> matrix
   */
  public double normF() {
    double s = 0.0;
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        double v = get(i,j);
        s += v * v;
      }
    }
    return Math.sqrt(s);
  }

  //----------------------------------------------------------------------------

  /**
   * Returns the trace (sum of diagonal elements).
   *
   * @return trace of <code>this</code> matrix
   */
  public double trace() {
    double trace = 0.0;
    final int limit = Math.min(rows(), cols());
    for (int i = 0; i < limit; ++i) { trace += get(i,i); }
    return trace;
  }

  /**
   * Returns the trace of the result given by <code>this</code> matrix
   * multiplied with matrix <code>m</code> on the right.
   *
   * The row number of <code>m</code> has to be equal to the column number
   * of <code>this</code> matrix.
   *
   * @param m matrix multiplying on the right (not <code>null</code>)
   * @return trace of <code>this</code> times <code>m</code>
   */
  public double traceMul(Matrix m) {
    assert (m != null && cols() == m.rows());
    double trace = 0.0;
    final int limit = Math.min(rows(), m.cols());
    for (int i = 0; i < limit; ++i) {
      for (int j = 0; j < cols(); ++j) {
        trace += get(i,j) * m.get(j,i);
      }
    }
    return trace;
  }

  /**
   * Returns the product of the diagonal elements.
   *
   * @return product of diagonal elements
   */
  public double prodDiag() {
    final int limit = Math.min(rows(), cols());
    double prod = 1.0;
    for (int i = 0; i < limit; ++i) { prod *= get(i,i); }
    return prod;
  }

  /**
   * Returns the sum of the rows (in <code>result</code>).
   *
   * The length of <code>result</code> has to be the column number of
   * <code>this</code> matrix.
   *
   * @param result storage of the result (not <code>null</code>)
   * @return sum of rows
   */
  public Vector rowSum(Vector result) {
    assert (result != null && result.length() == cols());
    result.setToZero();
    for (int j = 0; j < cols(); ++j) {
      double s = 0.0;
      for (int i = 0; i < rows(); ++i) { s += get(i,j); }
      result.set(j, s);
    }
    return result;
  }

  /**
   * Returns the sum of the rows (in new vector).
   *
   * @return sum of rows
   */
  public Vector rowSum() {
    return rowSum(Vector.create(cols()));
  }

  /**
   * Returns the sum of the columns (in <code>result</code>).
   *
   * The length of <code>result</code> has to be the row number of
   * <code>this</code> matrix.
   *
   * @param result storage of the result (not <code>null</code>)
   * @return sum of columns
   */
  public Vector colSum(Vector result) {
    assert (result != null && result.length() == rows());
    result.setToZero();
    for (int i = 0; i < rows(); ++i) {
      double s = 0.0;
      for (int j = 0; j < cols(); ++j) { s += get(i,j); }
      result.set(i, s);
    }
    return result;
  }

  /**
   * Returns the sum of the columns (in new vector).
   *
   * @return sum of columns
   */
  public Vector colSum() {
    return colSum(Vector.create(rows()));
  }

  //----------------------------------------------------------------------------

  /**
   * Cholesky decomposition of a (symmetric) positive-definite matrix.
   *
   * The <code>result</code> matrix is assumed to be zero above the diagonal
   * and has the same size as <code>this</code>.
   *
   * @param result storage of the result
   *        (not <code>null</code> and not equal to <code>this</code>)
   * @return <code>result</code> holding the lower-triangular Cholesky factor
   */
  public Matrix choleskyL(Matrix result) {
    assert (rows() == cols());
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    final int n = rows();
    double Ljj, Lij, v;
    for (int j = 0; j < n; ++j) {
      Ljj = get(j,j);
      for (int k = 0; k < j; ++k) {
        v = result.get(j,k);
        Ljj -= v * v;
      }
      Ljj = Math.sqrt(Ljj);
      result.set(j, j, Ljj);
            
      for (int i = j+1; i < n; ++i) {
        Lij = get(i,j);
        for (int k = 0; k < j; ++k) {
          Lij -= result.get(i,k) * result.get(j,k);
        }
        result.set(i, j, Lij / Ljj);
      }
    }
    return result;
  }

  /**
   * Cholesky decomposition of a (symmetric) positive-definite matrix.
   *
   * @return the lower-triangular Cholesky factor in a new matrix
   */
  public Matrix choleskyL() {
    return choleskyL(zero(rows(), cols()));
  }

  /**
   * LDL (Cholesky) decomposition of a (symmetric) positive-definite matrix.
   * Returns the result in matrix <code>L</code> and vector <code>D</code>.
   *
   * The <code>L</code> matrix is assumed to be zero above the diagonal
   * and has the same size as <code>this</code>.
   * The diagonal elements of <code>L</code> will be set to one.
   * The length of vector D has to be equal to the row/column number of
   * <code>this</code> matrix.
   *
   * @param L will be set to the lower-triangular LDL factor
   *        (not <code>null</code> and not equal to <code>this</code>)
   * @param D will be set to the diagonal elements of the diagonal LDL factor
   *        (not <code>null</code>)
   */
  public void choleskyLD(Matrix L, Vector D) {
    assert (rows() == cols());
    assert (L != null && L.rows() == rows() && L.cols() == cols());
    assert (D != null && D.length() == rows());
    final int n = rows();
    double Dj, Lij, v;
    for (int j = 0; j < n; ++j) {
      Dj = get(j,j);
      for (int k = 0; k < j; ++k) {
        v = L.get(j,k);
        Dj -= v * v * D.get(k);
      }
      D.set(j, Dj);
      L.set(j, j, 1.0);
            
      for (int i = j+1; i < n; ++i) {
        Lij = get(i,j);
        for (int k = 0; k < j; ++k) {
          Lij -= L.get(i,k) * L.get(j,k) * D.get(k);
        }
        L.set(i, j, Lij / Dj);
      }
    }
  }

  /**
   * LDL (Cholesky) decomposition of a (symmetric) positive-definite matrix.
   *
   * @return <code>{L,D}</code> matrices, where <code>L</code> is the
   * lower-triangular and <code>D</code> is the diagonal factor
   */
  public Matrix[] choleskyLD() {
    assert (rows() == cols());
    Matrix L = zero(rows(), rows());
    Vector D = Vector.create(rows());
    choleskyLD(L, D);
    return new Matrix[]{L, diag(D)};
  }

  //----------------------------------------------------------------------------

  /**
   * QR decomposition of an arbitrary matrix using Hauseholder transformations.
   *
   * Matrix <code>R</code> has to have the same size as <code>this</code>.
   * Vector <code>tmpV</code> is a temporary storage with length not smaller
   * than <code>min(rows()-1, cols())</code>.
   *
   * @param R will be set to the upper-triangular factor
   *          (not <code>null</code> with size <code>rows()</code> x <code>cols()</code>
   * @param tmpV temporary vector
   *             (not <code>null</code> with size >= <code>min(rows()-1, cols())</code>
   * @return matrix <code>R</code>
   */
  public Matrix QR(Matrix R, Vector tmpV) {
    QR(null, R, false, tmpV);
    return R;
  }

  /**
   * QR decomposition of an arbitrary matrix using Hauseholder transformations.
   *
   * Matrix <code>Q</code> has to be square of size <code>rows()</code> and
   * matrix <code>R</code> has to have the same size as <code>this</code>.
   * Vector <code>tmpV</code> is a temporary storage with length not smaller
   * than <code>min(rows()-1, cols())</code>.
   *
   * @param Q will be set to the orthogonal factor
   *          (not <code>null</code> with size <code>rows()</code> x <code>rows()</code>)
   * @param R will be set to the upper-triangular factor
   *          (not <code>null</code> with size <code>rows()</code> x <code>cols()</code>
   * @param tmpV temporary vector
   *             (not <code>null</code> with size >= <code>min(rows()-1, cols())</code>
   */
  public void QR(Matrix Q, Matrix R, Vector tmpV) {
    QR(Q, R, true, tmpV);
  }

  /**
   * QR decomposition of an arbitrary matrix using Hauseholder transformations.
   *
   * @return <code>{Q,R}</code> matrices, where <code>Q</code> is the orthogonal
   *         and <code>R</code> is the upper-triangular factor
   */
  public Matrix[] QR() {
    Matrix Q = create(rows(), rows());
    Matrix R = zero(rows(), cols());
    QR(Q, R, true, Vector.create(Math.min(rows()-1, cols())));
    return new Matrix[]{Q, R};
  }

  private void QR(Matrix Q, Matrix R, boolean isComputeQ, Vector tmpV) {
    assert (Q.rows() == rows() && Q.cols() == rows());
    assert (R.rows() == rows() && R.cols() == cols());

    final int rows = rows(), cols = cols();
    final int t = Math.min(rows-1, cols);
    assert (t <= tmpV.length());

    copy(R); // this -> R
    for (int k = 0; k < t; ++k) {
      double norm = 0.0;
      for (int i = k; i < rows; ++i) { norm = hypot(norm, R.get(i, k)); }

      if (norm != 0.0) {
        if (R.get(k, k) < 0) norm = -norm;

        for (int i = k; i < rows; ++i) { R.set(i, k, R.get(i, k) / norm); }
        R.set(k, k, R.get(k, k) + 1.0); // 1 <= R[k,k] <= 2

        for (int j = k+1; j < cols; ++j) {
          double s = 0.0;
          for (int i = k; i < rows; ++i) {
            s += R.get(i, k) * R.get(i, j);
          }
          s = -s / R.get(k, k);
          for (int i = k; i < rows; ++i) {
            R.set(i, j, R.get(i, j) + s*R.get(i, k));
          }
        }
      }
      tmpV.set(k, -norm);
    }

    if (isComputeQ) {
      Q.setToEye();
      for (int k = t-1; k >= 0; --k) {
        for (int j = k; j < rows; ++j) {
          if (R.get(k, k) != 0.0) {
            double s = 0.0;
            for (int i = k; i < rows; ++i) {
              s += R.get(i, k) * Q.get(i, j);
            }
            s = -s / R.get(k, k);
            for (int i = k; i < rows; ++i) {
              Q.set(i, j, Q.get(i, j) + s*R.get(i, k));
            }
          }
        }
      }
    }

    for (int k = 0; k < t; ++k) {
      R.set(k, k, tmpV.get(k));
      for (int i = k+1; i < rows; ++i) { R.set(i, k, 0.0); }
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

  // TODO: unify these to a general det method
  // use recursive formula for n <= 5, and QR above

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

  // TODO: unify these into a general inv method

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
    return result.div(det, result);
  }

  /**
   * @return inverse of a 2x2 matrix (placed into a new matrix)
   *         or null if the determinant is zero
   */
  public Matrix inv2x2() {
    return inv2x2(null);
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
    return result.div(det, result);
  }

  /**
   * @return inverse of a 3x3 matrix (placed into a new matrix)
   *         or null if the determinant is zero
   */    
  public Matrix inv3x3() {
    return inv3x3(null);
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
    invD.reciproc(invD);
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
   * Returns the matrix transpose.
   * Data is not moved, but the get/set methods are redefined,
   * so this operation cannot be done entirely in-place.
   *
   * @return transpose of the matrix (in new object)
   */
  abstract public Matrix T();

  //----------------------------------------------------------------------------

  @Override
  public String toString() {
    String str = "[";
    for (int i = 0; i < rows(); ++i) {
      if (0 != i) str += "; ";
      for (int j = 0; j < cols(); ++j) {
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
