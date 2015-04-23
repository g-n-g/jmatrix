package jvecmat;

import java.util.Random;

/**
 * Abstract matrix representation.
 */
public abstract class Matrix implements VecMat {

  /**
   * New row constant for the {@link Matrix#create(double...)} constructor.
   */
  public final static double NR = Double.NaN;

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
   * Creates a matrix of the given <code>values</code>.
   * New rows can be indicated by the {@link Matrix#NR} constant.
   *
   * The <code>values</code> should not contain any <code>Double.NaN</code> value.
   * Empty rows are silently ignored.
   *
   * @param values matrix elements and {@link Matrix#NR} new row indicators
   * @return new matrix with elements from <code>values</code>
   * @throws IllegalArgumentException if row lengths mismatch
   */
  public static Matrix create(double... values) {
    if (values.length == 0) { return create(0, 0); }

    int rows = 0;
    int cols = 0;
    boolean isFirstNR = true;
    for (int i = 1, iprev = 0; i <= values.length; ++i) {
      if (i == values.length || Double.isNaN(values[i])) {
        if (Double.isNaN(values[i-1])) { // empty row
          ++iprev;
        } else { // non-empty row
          ++rows;
          if (isFirstNR) {
            isFirstNR = false;
            cols = i;
          } else if (cols != i-iprev-1) {
            throw new IllegalArgumentException("Length mismatch of row " + rows + "!");
          }
          iprev = i;
        }
      }
    }
    Matrix result = create(rows, cols);
    int i = 0, j = 0;
    for (double value : values) {
      if (!Double.isNaN(value)) {
        result.set(i, j, value);
        if (++j == cols) { j = 0; ++i; } // new row
      }
    }
    return result;
  }

  /**
   * Creates a matrix of given size and sets all elements to <code>value</code>.
   *
   * @param rows number of rows
   * @param cols number of columns
   * @param value element initializer value
   * @return matrix with size <code>rows</code> x <code>cols</code>
   *         having elements set to <code>value</code>
   */
  public static Matrix scalars(int rows, int cols, double value) {
    Matrix m = create(rows, cols);
    m.setToScalars(value);
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
  public static Matrix zeros(int rows, int cols) {
    return scalars(rows, cols, 0.0);
  }

  /**
   * Creates a matrix of size <code>rows</code> x <code>cols</code>
   * and initializes its elements to zero.
   *
   * @param rows number of rows
   * @param cols number of columns
   * @return all one matrix of size <code>rows</code> x <code>cols</code>
   */
  public static Matrix ones(int rows, int cols) {
    return scalars(rows, cols, 1.0);
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
    Matrix m = Matrix.create(rows, cols);
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
    Matrix m = Matrix.create(rows, cols);
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
   * Returns a copy of the matrix placed into <code>result</code>.
   * If <code>result</code> is <code>null</code>, a new object is created.
   * The operation is skipped if <code>result</code> is equal to <code>this</code>.
   *
   * @param result appropriately sized storage for the copy (not <code>null</code>)
   * @return copy of the matrix
   */
  public Matrix copy(Matrix result) {
    assert(result != null && result.rows() == rows() && result.cols() == cols());
    if (result != this) {
      for (int i = 0; i < rows(); ++i) {
        for (int j = 0; j < cols(); ++j) {
          result.set(i, j, get(i,j));
        }
      }
    }
    return result;
  }

  @Override
  public Matrix copy() {
    return copy(create(rows(), cols()));
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
  public Matrix setToScalars(double c) {
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        set(i, j, c);
      }
    }
    return this;
  }

  @Override
  public Matrix setToZeros() {
    return setToScalars(0.0);
  }

  @Override
  public Matrix setToOnes() {
    return setToScalars(1.0);
  }

  /**
   * Set the diagonal elements to one and all other elements to zero.
   * If the matrix is square, this sets the identity matrix.
   *
   * @return <code>this</code> matrix
   */
  public Matrix setToEye() {
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        set(i, j, i==j ? 1.0 : 0.0);
      }
    }
    return this;
  }

  @Override
  public Matrix setToRand(Random rng) {
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        set(i, j, rng.nextDouble());
      }
    }
    return this;
  }

  @Override
  public Matrix setToRandN(Random rng) {
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
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

  //----------------------------------------------------------------------------

  /**
   * Returns the lower triangular part of the matrix (in <code>result</code>).
   * This includes all elements under and on the diagonal selected by
   * <code>offset</code>. For positive/negative values of <code>offset</code>,
   * the selected diagonal is above/under the main diagonal, respectively.
   * 
   * Matrix <code>result</code> has to have the same size as <code>this</code>
   * matrix.
   *
   * @param offset offset of the main diagonal
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> holding the lower triangular part
   */
  public Matrix getTriL(int offset, Matrix result) {
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    int jend = Math.min(offset+1, cols());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < jend; ++j) {
        result.set(i, j, get(i,j));
      }
      if (jend >= 0) {
        for (int j = jend; j < cols(); ++j) {
          result.set(i, j, 0.0);
        }
      }
      if (jend < cols()) { ++jend; }
    }
    return result;
  }

  /**
   * Returns the lower triangular part of the matrix (in new matrix).
   * This includes all elements under and on the diagonal selected by
   * <code>offset</code>. For positive/negative values of <code>offset</code>,
   * the selected diagonal is above/under the main diagonal, respectively.
   * 
   * @param offset offset of the main diagonal
   * @return new matrix holding the lower triangular part
   */
  public Matrix getTriL(int offset) {
    return getTriL(offset, create(rows(), cols()));
  }

  /**
   * Returns the lower triangular part of the matrix (in <code>result</code>).
   * This includes all elements under and on the main diagonal.
   * 
   * Matrix <code>result</code> has to have the same size as <code>this</code>
   * matrix.
   *
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> holding the lower triangular part
   */
  public Matrix getTriL(Matrix result) {
    return getTriL(0, result);
  }

  /**
   * Returns the lower triangular part of the matrix (in new matrix).
   * This includes all elements under and on the main diagonal.
   * 
   * @return new matrix holding the lower triangular part
   */
  public Matrix getTriL() {
    return getTriL(0);
  }

  /**
   * Returns the upper triangular part of the matrix (in <code>result</code>).
   * This includes all elements under and on the diagonal selected by
   * <code>offset</code>. For positive/negative values of <code>offset</code>,
   * the selected diagonal is above/under the main diagonal, respectively.
   * 
   * Matrix <code>result</code> has to have the same size as <code>this</code>
   * matrix.
   *
   * @param offset offset of the main diagonal
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> holding the upper triangular part
   */
  public Matrix getTriU(int offset, Matrix result) {
    return T().getTriL(-offset, result).T();
  }

  /**
   * Returns the upper triangular part of the matrix (in new matrix).
   * This includes all elements under and on the diagonal selected by
   * <code>offset</code>. For positive/negative values of <code>offset</code>,
   * the selected diagonal is above/under the main diagonal, respectively.
   * 
   * @param offset offset of the main diagonal
   * @return new matrix holding the upper triangular part
   */
  public Matrix getTriU(int offset) {
    return getTriU(offset, create(rows(), cols()));
  }

  /**
   * Returns the upper triangular part of the matrix (in <code>result</code>).
   * This includes all elements under and on the main diagonal.
   * 
   * Matrix <code>result</code> has to have the same size as <code>this</code>
   * matrix.
   *
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> holding the upper triangular part
   */
  public Matrix getTriU(Matrix result) {
    return getTriU(0, result);
  }

  /**
   * Returns the upper triangular part of the matrix (in new matrix).
   * This includes all elements under and on the main diagonal.
   * 
   * @return new matrix holding the upper triangular part
   */
  public Matrix getTriU() {
    return getTriU(0);
  }

  //----------------------------------------------------------------------------

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

  /**
   * Matrix-constant addition (in <code>result</code>).
   * Adds constant <code>c</code> to all elements of <code>this</code> matrix.
   *
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param c constant to add
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the values of
   *         <code>this</code> shifted by <code>c</code>
   */
  public Matrix add(double c, Matrix result) {
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < cols(); ++j) {
        result.set(i, j, get(i,j) + c);
      }
    }
    return result;
  }

  @Override
  public Matrix add(double c) {
    return add(c, create(rows(), cols()));
  }

  /**
   * Matrix-vector addition along rows (in <code>result</code>).
   * Adds vector <code>v</code> to all rows of <code>this</code> matrix.
   *
   * The length of vector <code>v</code> has to be equal to the column number
   * of <code>this</code> matrix.
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param v vector to add to each row (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the rows of
   *         <code>this</code> shifted by vector <code>v</code>
   */
  public Matrix addRow(Vector v, Matrix result) {
    assert (v != null && v.length() == cols());
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int j = 0; j < cols(); ++j) {
      double value = v.get(j);
      for (int i = 0; i < rows(); ++i) { result.set(i, j, get(i,j) + value); }
    }
    return result;
  }

  /**
   * Matrix-vector addition along rows (in new matrix).
   * Adds vector <code>v</code> to all rows of <code>this</code> matrix.
   *
   * The length of vector <code>v</code> has to be equal to the column number
   * of <code>this</code> matrix.
   *
   * @param v vector to add to each row (not <code>null</code>)
   * @return new matrix having the rows of
   *         <code>this</code> shifted by vector <code>v</code>
   */
  public Matrix addRow(Vector v) {
    return addRow(v, create(rows(), cols()));
  }

  /**
   * Matrix-vector addition along columns (in <code>result</code>).
   * Adds vector <code>v</code> to all columns of <code>this</code> matrix.
   *
   * The length of vector <code>v</code> has to be equal to the row number
   * of <code>this</code> matrix.
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param v vector to add to each column (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the columns of
   *         <code>this</code> shifted by vector <code>v</code>
   */
  public Matrix addCol(Vector v, Matrix result) {
    assert (v != null && v.length() == rows());
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      double value = v.get(i);
      for (int j = 0; j < cols(); ++j) { result.set(i, j, get(i,j) + value); }
    }
    return result;
  }

  /**
   * Matrix-vector addition along columns (in new matrix).
   * Adds vector <code>v</code> to all columns of <code>this</code> matrix.
   *
   * The length of vector <code>v</code> has to be equal to the row number
   * of <code>this</code> matrix.
   *
   * @param v vector to add to each column (not <code>null</code>)
   * @return new matrix having the columns of
   *         <code>this</code> shifted by vector <code>v</code>
   */
  public Matrix addCol(Vector v) {
    return addCol(v, create(rows(), cols()));
  }

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
    return sub(m, create(rows(), cols()));
  }

  /**
   * Matrix-constant subtraction (in <code>result</code>).
   * Subtracts constant <code>c</code> from all elements of <code>this</code>
   * matrix.
   *
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param c constant to subtract
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the values of
   *         <code>this</code> shifted by <code>-c</code>
   */
  public Matrix sub(double c, Matrix result) {
    return add(-c, result);
  }

  @Override
  public Matrix sub(double c) {
    return sub(c, create(rows(), cols()));
  }

  /**
   * Matrix-vector subtraction along rows (in <code>result</code>).
   * Subtracts vector <code>v</code> from all rows of <code>this</code> matrix.
   *
   * The length of vector <code>v</code> has to be equal to the column number
   * of <code>this</code> matrix.
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param v vector to subtract from each row (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the rows of
   *         <code>this</code> shifted by vector <code>-v</code>
   */
  public Matrix subRow(Vector v, Matrix result) {
    assert (v != null && v.length() == cols());
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int j = 0; j < cols(); ++j) {
      double value = v.get(j);
      for (int i = 0; i < rows(); ++i) { result.set(i, j, get(i,j) - value); }
    }
    return result;
  }

  /**
   * Matrix-vector subtraction along rows (in new matrix).
   * Subtracts vector <code>v</code> from all rows of <code>this</code> matrix.
   *
   * The length of vector <code>v</code> has to be equal to the column number
   * of <code>this</code> matrix.
   *
   * @param v vector to subtract from each row (not <code>null</code>)
   * @return new matrix having the rows of
   *         <code>this</code> shifted by vector <code>-v</code>
   */
  public Matrix subRow(Vector v) {
    return subRow(v, create(rows(), cols()));
  }

  /**
   * Matrix-vector subtraction along columns (in <code>result</code>).
   * Subtracts vector <code>v</code> from all columns of <code>this</code>
   * matrix.
   *
   * The length of vector <code>v</code> has to be equal to the row number
   * of <code>this</code> matrix.
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param v vector to subtract from each column (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the columns of
   *         <code>this</code> shifted by vector <code>-v</code>
   */
  public Matrix subCol(Vector v, Matrix result) {
    assert (v != null && v.length() == rows());
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      double value = v.get(i);
      for (int j = 0; j < cols(); ++j) { result.set(i, j, get(i,j) - value); }
    }
    return result;
  }

  /**
   * Matrix-vector subtraction along columns (in new matrix).
   * Subtracts vector <code>v</code> from all columns of <code>this</code>
   * matrix.
   *
   * The length of vector <code>v</code> has to be equal to the row number
   * of <code>this</code> matrix.
   *
   * @param v vector to subtract from each column (not <code>null</code>)
   * @return new matrix having the columns of
   *         <code>this</code> shifted by vector <code>-v</code>
   */
  public Matrix subCol(Vector v) {
    return subCol(v, create(rows(), cols()));
  }

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
    return div(c, create(rows(), cols()));
  }

  /**
   * Matrix-vector elementwise division along rows (in <code>result</code>).
   * Divides all rows of <code>this</code> matrix elementwise
   * by vector <code>v</code>.
   *
   * The length of vector <code>v</code> has to be equal to the column number
   * of <code>this</code> matrix.
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param v vector to divide each row (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the rows of
   *         <code>this</code> divided by the elements of vector <code>v</code>
   */
  public Matrix divRow(Vector v, Matrix result) {
    assert (v != null && v.length() == cols());
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int j = 0; j < cols(); ++j) {
      double value = v.get(j);
      for (int i = 0; i < rows(); ++i) { result.set(i, j, get(i,j) / value); }
    }
    return result;
  }

  /**
   * Matrix-vector elementwise division along rows (in new matrix).
   * Divides all rows of <code>this</code> matrix elementwise
   * by vector <code>v</code>.
   *
   * The length of vector <code>v</code> has to be equal to the column number
   * of <code>this</code> matrix.
   *
   * @param v vector to divide each row (not <code>null</code>)
   * @return new matrix having the rows of
   *         <code>this</code> divided by the elements of vector <code>v</code>
   */
  public Matrix divRow(Vector v) {
    return divRow(v, create(rows(), cols()));
  }

  /**
   * Matrix-vector elementwise division along columns (in <code>result</code>).
   * Divides all columns of <code>this</code> matrix elementwise
   * by vector <code>v</code>.
   *
   * The length of vector <code>v</code> has to be equal to the row number
   * of <code>this</code> matrix.
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param v vector to divide each column (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the columns of
   *         <code>this</code> divided by the elements of vector <code>v</code>
   */
  public Matrix divCol(Vector v, Matrix result) {
    assert (v != null && v.length() == rows());
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      double value = v.get(i);
      for (int j = 0; j < cols(); ++j) { result.set(i, j, get(i,j) / value); }
    }
    return result;
  }

  /**
   * Matrix-vector elementwise division along columns (in new matrix).
   * Divides all columns of <code>this</code> matrix elementwise
   * by vector <code>v</code>.
   *
   * The length of vector <code>v</code> has to be equal to the row number
   * of <code>this</code> matrix.
   *
   * @param v vector to divide each column (not <code>null</code>)
   * @return new matrix having the columns of
   *         <code>this</code> divided by the elements of vector <code>v</code>
   */
  public Matrix divCol(Vector v) {
    return divCol(v, create(rows(), cols()));
  }

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

  /**
   * Matrix-vector elementwise multiplication along rows (in <code>result</code>).
   * Multiplies all rows of <code>this</code> matrix elementwise
   * by vector <code>v</code>.
   *
   * The length of vector <code>v</code> has to be equal to the column number
   * of <code>this</code> matrix.
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param v vector to multiply each row (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the rows of
   *         <code>this</code> multiplied by the elements of vector <code>v</code>
   */
  public Matrix mulRow(Vector v, Matrix result) {
    assert (v != null && v.length() == cols());
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int j = 0; j < cols(); ++j) {
      double value = v.get(j);
      for (int i = 0; i < rows(); ++i) { result.set(i, j, get(i,j) * value); }
    }
    return result;
  }

  /**
   * Matrix-vector elementwise multiplication along rows (in new matrix).
   * Multiplies all rows of <code>this</code> matrix elementwise
   * by vector <code>v</code>.
   *
   * The length of vector <code>v</code> has to be equal to the column number
   * of <code>this</code> matrix.
   *
   * @param v vector to multiply each row (not <code>null</code>)
   * @return new matrix having the rows of
   *         <code>this</code> multiplied by the elements of vector <code>v</code>
   */
  public Matrix mulRow(Vector v) {
    return mulRow(v, create(rows(), cols()));
  }

  /**
   * Matrix-vector elementwise multiplication along columns (in <code>result</code>).
   * Multiplies all columns of <code>this</code> matrix elementwise
   * by vector <code>v</code>.
   *
   * The length of vector <code>v</code> has to be equal to the row number
   * of <code>this</code> matrix.
   * Matrix <code>result</code> has to have the same size as <code>this</code>.
   * The <code>result</code> parameter can be also set to <code>this</code>
   * providing in-place operation.
   *
   * @param v vector to multiply each column (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>)
   * @return <code>result</code> having the columns of
   *         <code>this</code> multiplied by the elements of vector <code>v</code>
   */
  public Matrix mulCol(Vector v, Matrix result) {
    assert (v != null && v.length() == rows());
    assert (result != null && result.rows() == rows() && result.cols() == cols());
    for (int i = 0; i < rows(); ++i) {
      double value = v.get(i);
      for (int j = 0; j < cols(); ++j) { result.set(i, j, get(i,j) * value); }
    }
    return result;
  }

  /**
   * Matrix-vector elementwise multiplication along columns (in new matrix).
   * Multiplies all columns of <code>this</code> matrix elementwise
   * by vector <code>v</code>.
   *
   * The length of vector <code>v</code> has to be equal to the row number
   * of <code>this</code> matrix.
   *
   * @param v vector to multiply each column (not <code>null</code>)
   * @return new matrix having the columns of
   *         <code>this</code> multiplied by the elements of vector <code>v</code>
   */
  public Matrix mulCol(Vector v) {
    return mulCol(v, create(rows(), cols()));
  }

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
   * @see Vector#mul(Matrix, Vector)
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
   * @see Vector#mul(Matrix)
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
   * @see Vector#mulD(Matrix, Matrix)
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
   * @see Vector#mulD(Matrix)
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

  /**
   * Permutes the columns of the matrix, i.e. multiplying the matrix from
   * the right with a permutation matrix represented by <code>p</code>
   * (in <code>result</code>).
   *
   * The length of permutation <code>p</code> has to be equal to the column
   * number of <code>this</code> matrix.
   * Furthermore, <code>result</code> has to have the same size
   * as <code>this</code> matrix.
   *
   * @param p permutation (not <code>null</code>)
   * @param result storage of the result (not <code>null</code>
   *                                      and not equal to <code>this</code>)
   * @return <code>result</code> having <code>this</code> matrix multiplied from
   *         the right with a permutation matrix represented by <code>p</code>
   * @see Permutation#mul(Matrix, Matrix)
   */
  public Matrix mul(Permutation p, Matrix result) {
    assert (p != null && p.length() == cols());
    assert (result != null && result != this &&
            result.rows() == rows() && result.cols() == cols());
    for (int j = 0; j < cols(); ++j) {
      int col = p.get(j);
      for (int i = 0; i < rows(); ++i) {
        result.set(i, j, get(i,col));
      }
    }
    return result;
  }

  /**
   * Permutes the columns of the matrix, i.e. multiplying the matrix from
   * the right with a permutation matrix represented by <code>p</code>
   * (in new matrix).
   *
   * The length of permutation <code>p</code> has to be equal to the column
   * number of <code>this</code> matrix.
   *
   * @param p permutation (not <code>null</code>)
   * @return new matrix being equal to <code>this</code> matrix multiplied from
   *         the right with a permutation matrix represented by <code>p</code>
   * @see Permutation#mul(Matrix)
   */
  public Matrix mul(Permutation p) {
    return mul(p, create(rows(), cols()));
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
    result.setToZeros();
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
    result.setToZeros();
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
   * Returns the result in matrix <code>L</code>
   * having the same size as <code>this</code>.
   *
   * @param L will be set to the lower triangular factor
   *        (not <code>null</code> and not equal to <code>this</code>)
   * @return <code>L</code> holding the lower triangular Cholesky factor
   * @throws UnsupportedOperationException if the matrix is not positive-definite
   */
  public Matrix choleskyL(Matrix L) {
    assert (rows() == cols());
    assert (L != null && L.rows() == rows() && L.cols() == cols());
    final int n = rows();
    double Ljj, Lij, v;
    for (int j = 0; j < n; ++j) {
      Ljj = get(j,j);
      for (int k = 0; k < j; ++k) {
        v = L.get(j,k);
        Ljj -= v * v;
      }
      if (Ljj <= 0.0) {
        throw new UnsupportedOperationException("Matrix has to be positive-definite.");
      }
      Ljj = Math.sqrt(Ljj);
      L.set(j, j, Ljj);
      Ljj = 1.0 / Ljj;
      if (Double.isInfinite(Ljj)) {
        throw new UnsupportedOperationException("Matrix has to be positive-definite.");
      }

      for (int i = 0; i < j; ++i) {
        L.set(i, j, 0.0); // zero above the diagonal
      }
      for (int i = j+1; i < n; ++i) {
        Lij = get(i,j);
        for (int k = 0; k < j; ++k) {
          Lij -= L.get(i,k) * L.get(j,k);
        }
        L.set(i, j, Lij * Ljj);
      }
    }
    return L;
  }

  /**
   * Cholesky decomposition of a (symmetric) positive-definite matrix.
   *
   * @return the lower triangular Cholesky factor in a new matrix
   * @throws UnsupportedOperationException if the matrix is not positive-definite
   */
  public Matrix choleskyL() {
    return choleskyL(zeros(rows(), cols()));
  }

  /**
   * LDL (Cholesky) decomposition of a (symmetric) positive-definite matrix.
   * Returns the result in matrix <code>L</code> and vector <code>D</code>.
   *
   * The <code>L</code> matrix should have the same size as <code>this</code>.
   * The diagonal elements of <code>L</code> will be set to one.
   * The length of vector D has to be equal to the row/column number of
   * <code>this</code> matrix.
   *
   * @param L will be set to the lower triangular LDL factor
   *        (not <code>null</code> and not equal to <code>this</code>)
   * @param D will be set to the diagonal elements of the diagonal LDL factor
   *        (not <code>null</code>)
   * @throws UnsupportedOperationException if the matrix is not positive-definite
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
      if (Dj <= 0.0) {
        throw new UnsupportedOperationException("Matrix has to be positive-definite.");
      }
      D.set(j, Dj);
      Dj = 1.0 / Dj;
      if (Double.isInfinite(Dj)) {
        throw new UnsupportedOperationException("Matrix has to be positive-definite.");
      }
      L.set(j, j, 1.0);
            
      for (int i = 0; i < j; ++i) {
        L.set(i, j, 0.0); // zero above the diagonal
      }
      for (int i = j+1; i < n; ++i) {
        Lij = get(i,j);
        for (int k = 0; k < j; ++k) {
          Lij -= L.get(i,k) * L.get(j,k) * D.get(k);
        }
        L.set(i, j, Lij * Dj);
      }
    }
  }

  /**
   * LDL (Cholesky) decomposition of a (symmetric) positive-definite matrix.
   *
   * @return <code>{L,D}</code> matrices, where <code>L</code> is the
   *         lower triangular and <code>D</code> is the diagonal factor
   * @throws UnsupportedOperationException if the matrix is not positive-definite
   */
  public Matrix[] choleskyLD() {
    assert (rows() == cols());
    Matrix L = zeros(rows(), rows());
    Vector D = Vector.create(rows());
    choleskyLD(L, D);
    return new Matrix[]{L, diag(D)};
  }

  //----------------------------------------------------------------------------

  /**
   * QR decomposition of an arbitrary matrix using Hauseholder transformations.
   *
   * Matrix <code>Q</code> can be <code>null</code> in which case its computation
   * is omitted. Otherwise, <code>Q</code> has to be a square matrix of size
   * <code>rows()</code>. Matrix <code>R</code> has to have the same size as
   * <code>this</code> and can be set to <code>this</code> providing in-place
   * operation. Vector <code>tmpV</code> is a temporary storage with length not
   * smaller than <code>min(rows()-1, cols())</code>.
   * Matrices <code>Q</code> and <code>R</code> cannot be the same.
   *
   * @param Q will be set to the orthogonal factor
   *          (<code>null</code> or with size <code>rows()</code> x <code>rows()</code>)
   * @param R will be set to the upper-triangular factor
   *          (not <code>null</code> with size <code>rows()</code> x <code>cols()</code>
   * @param tmpV temporary vector
   *             (not <code>null</code> with size >= <code>min(rows()-1, cols())</code>
   */
  public void QR(Matrix Q, Matrix R, Vector tmpV) {
    assert (Q == null || (Q.rows() == rows() && Q.cols() == rows()));
    assert (R != null && R != Q && R.rows() == rows() && R.cols() == cols());

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

    if (Q != null) { // compute Q only if requested
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

      // ensure R is upper triangular
      for (int i = k+1; i < rows; ++i) { R.set(i, k, 0.0); }
    }
  }

  /**
   * Computes the hypotenuse of a right triangle
   * having legs <code>x</code> and <code>y</code>.
   */
  private static double hypot(double x, double y) {
    final double absX = Math.abs(x);
    final double absY = Math.abs(y);

    double r = 0.0;
    if (absX > absY) {
      r = y/x;
      r = absX * Math.sqrt(1.0 + r*r);
    }
    else if (y != 0.0) {
      r = x/y;
      r = absY * Math.sqrt(1.0 + r*r);
    }
    return r;
  }

  /**
   * QR decomposition of an arbitrary matrix using Hauseholder transformations.
   *
   * @return <code>{Q,R}</code> matrices, where <code>Q</code> is the orthogonal
   *         and <code>R</code> is the upper-triangular factor
   */
  public Matrix[] QR() {
    Matrix Q = create(rows(), rows());
    Matrix R = zeros(rows(), cols());
    QR(Q, R, Vector.create(Math.min(rows()-1, cols())));
    return new Matrix[]{Q, R};
  }

  //----------------------------------------------------------------------------

  /**
   * P'LU decomposition of an arbitrary matrix
   * using Doolittle's method with partial pivoting.
   *
   * Matrix <code>L</code> will be set to the unit lower triangular factor
   * and matrix <code>U</code> will be set to the upper triangular factor.
   *
   * Matrix <code>L</code> can be set equal to <code>U</code> in which case
   * the unit diagonal elements of <code>L</code> will not be stored.
   * Both <code>L</code> and <code>U</code> can be set to <code>this</code>
   * supporting in-place operation.
   *
   * @param L the unit lower triangular factor, might be row permuted if
   *          <code>P</code> is set to <code>null</code> (not <code>null</code>
   *          of size <code>rows()</code> x <code>min(rows(), cols())</code>)
   * @param U the upper triangular factor (not <code>null</code>
   *          of size <code>min(rows(), cols())</code> x <code>cols()</code>)
   * @param P row permutations of <code>L</code>
   *          (not <code>null</code> with length <code>rows()</code>)
   * @return the number of row swaps (permutations)
   */
  public int LU(Matrix L, Matrix U, Permutation P) {
    assert (P != null);
    return LU(L, U, P, true);
  }

  private int LU(Matrix L, Matrix U, Permutation P, boolean isRowPivot) {
    if (rows() > cols()) {
      return T().LU(U.T(), L.T(), P, false);
    }

    assert (L != null && L.rows() == rows() && L.cols() == rows());
    assert (U != null && U.rows() == rows() && U.cols() == cols());

    copy(U);
    if (P != null) { P.setToEye(); }
    Matrix V = isRowPivot ? U : U.T(); // pivot matrix (sharing data with U)
    assert(P == null || P.length() == V.rows()); // (P == null) is for internal use only

    int nperms = 0;
    for (int i = 0; i < rows(); ++i) {
      // find the pivot
      int p = i;
      double pval = V.get(i,i);
      for (int k = i+1; k < rows(); ++k) {
        double val = V.get(k,i);
        if (val > pval) {
          p = k;
          pval = val;
        }
      }

      // skip the column if the pivot is zero
      if (pval == 0.0) { continue; }

      // swap rows if necessary
      if (p != i) {
        for (int j = 0; j < V.cols(); ++j) {
          double val = V.get(i,j);
          V.set(i, j, V.get(p,j));
          V.set(p, j, val);
        }
        ++nperms;
        if (P != null) { P.swap(i, p); }
      }

      // update below the diagonal: Lij = Aij - sum_k(Lik*Ukj)
      for (int j = 0; j < i; ++j) {
        double sum = U.get(i,j);
        for (int k = 0; k < j; ++k) {
          sum -= L.get(i,k) * U.get(k,j);
        }
        L.set(i, j, sum / U.get(j,j));
      }

      // update above the diagonal: Uij = Aij - sum_k(Lik*Ukj)
      for (int j = i; j < cols(); ++j) {
        double sum = U.get(i,j);
        for (int k = 0; k < i; ++k) {
          sum -= L.get(i,k) * U.get(k,j);
        }
        U.set(i, j, sum);
      }
    }

    if (L != U) { // ensure U is upper triangular and L is unit lower triangular
      for (int i = 0; i < U.rows(); ++i) {
        for (int j = 0; j < i; ++j) { U.set(i, j, 0.0); }
        L.set(i, i, 1.0);
        for (int k = 0; k < i; ++k) { L.set(k, i, 0.0); }
      }
    }

    return nperms;
  }

  /**
   * P'LU decomposition of an arbitrary matrix
   * using Doolittle's method with partial pivoting.
   *
   * Matrix <code>P</code> will be set tp a permutation matrix,
   * matrix <code>L</code> will be set to the unit lower triangular factor
   * and matrix <code>U</code> will be set to the upper triangular factor.
   *
   * @return <code>{L,U,P}</code> matrices, where <code>P</code> is a permutation
   *         matrix, <code>L</code> is the unit lower triangular factor
   *         and <code>U</code> is the upper triangular factor
   */
  public Matrix[] LU() {
    final int n = Math.min(rows(), cols());
    Matrix L = create(rows(), n);
    Matrix U = create(n, cols());
    Permutation P = Permutation.eye(rows()); // TODO: LU will set to eye...
    LU(L, U, P);
    return new Matrix[]{L, U, P.toMatrix()};
  }

  //----------------------------------------------------------------------------

  /**
   * Returns the determinant of a square matrix. An LU decomposition is performed
   * for size larger than <code>3x3</code> using only the provided
   * <code>tmpLU</code> storage, which has to have the same size as
   * <code>this</code> matrix.
   *
   * @param tmpLU storage for the LU decomposition
   *              (not <code>null</code> if size is larger than 3x3)
   * @return determinant
   */
  public double det(Matrix tmpLU) {
    assert (rows() == cols());

    switch (rows()) {
    case 0: return 1.0;
    case 1: return get(0,0);
    case 2: return get(0,0)*get(1,1) - get(0,1)*get(1,0);
    case 3: return get(0,0) * (get(1,1)*get(2,2) - get(1,2)*get(2,1))
                 + get(0,1) * (get(1,2)*get(2,0) - get(1,0)*get(2,2))
                 + get(0,2) * (get(1,0)*get(2,1) - get(1,1)*get(2,0));
    default:
      assert (tmpLU.rows() == rows() && tmpLU.cols() == rows());
      int nperms = LU(tmpLU, tmpLU, null, true);
      double det = tmpLU.prodDiag();
      if ((nperms & 1) != 0) { det = -det; }
      return det;
    }
  }

  /**
   * Returns the determinant of a square matrix. An LU decomposition is performed
   * for size larger than <code>3x3</code>.
   *
   * @param result storage of the result
   *        (not <code>null</code> and not equal to <code>this</code>)
   * @return determinant
   */
  public double det() {
    assert (rows() == cols());
    return det((rows() > 3) ? create(rows(), rows()) : null);
  }

  //----------------------------------------------------------------------------

  /**
   * Computes the inverse of a square matrix (in <code>result</code>).
   *
   * Matrix <code>result</code> and <code>tmpM</code> should be square matrices
   * of the same size as <code>this</code>. The length of permutation <code>tmpP</code>
   * should match the side length of the square matrices.
   *
   * @param result storage of the result
   *               (not <code>null</code> and not equal to <code>this</code>)
   * @param tmpM temporary matrix (not <code>null</code> and not equal to
   *             <code>this</code>, having size of <code>rows() x cols()</code>)
   * @param tmpP temporary permutation
   *             (not <code>null</code> with size of <code>rows()</code>)
   * @return <code>result</code> holding the matrix inverse
   * @throws UnsupportedOperationException if the matrix is singular
   */
  public Matrix inv(Matrix result, Matrix tmpM, Permutation tmpP) {
    assert (rows() == cols());
    assert (result != null && result.rows() == rows() && result.cols() == cols());

    switch (rows()) {
    case 0:
      break;
    case 1: {
      double rdet = 1.0 / get(0,0);
      if (Double.isInfinite(rdet)) {
        throw new UnsupportedOperationException("Matrix is singular.");
      }
      result.set(0, 0, rdet);
      break;
    }
    case 2: {
      double a = get(0,0), b = get(0,1), c = get(1,0), d = get(1,1);
      double rdet = 1.0 / (a*d - b*c);
      if (Double.isInfinite(rdet)) {
        throw new UnsupportedOperationException("Matrix is singular.");
      }
      result.set(0, 0,  d); result.set(0, 1, -b);
      result.set(1, 0, -c); result.set(1, 1,  a);
      result.mul(rdet, result);
      break;
    }
    case 3: {
      double a = get(0,0), b = get(0,1), c = get(0,2);
      double d = get(1,0), e = get(1,1), f = get(1,2);
      double g = get(2,0), h = get(2,1), k = get(2,2);
      double A = e*k-f*h, D = c*h-b*k, G = b*f-c*e;
      double B = f*g-d*k, E = a*k-c*g, H = c*d-a*f;
      double C = d*h-e*g, F = g*b-a*h, K = a*e-b*d;
      double rdet = 1.0 / (a*A + b*B + c*C);
      if (Double.isInfinite(rdet)) {
        throw new UnsupportedOperationException("Matrix is singular.");
      }
      result.set(0, 0, A); result.set(0, 1, D); result.set(0, 2, G);
      result.set(1, 0, B); result.set(1, 1, E); result.set(1, 2, H);
      result.set(2, 0, C); result.set(2, 1, F); result.set(2, 2, K);
      result.mul(rdet, result);
      break;
    }
    default:
      LU(tmpM, tmpM, tmpP);
      tmpM.invLT(result, false); // inv(L)
      tmpM.T().invLT(result.T(), true); // inv(U)
      result.mulUL(tmpM); // tmpM = inv(U)*inv(L)
      tmpM.mul(tmpP, result); // result = tmpM*tmpP
    }
    return result;
  }

  /**
   * Inverts a square lower triangular matrix <code>L</code> in
   * </code>result</code>. If <code>mirror</code> is set, the result is also
   * copied into the upper half, otherwise <code>result</code> is kept untouched
   * above the diagonal.
   */
  private Matrix invLT(Matrix result, boolean useDiag) {
    for (int i = 0; i < rows(); ++i) {
      double Rii = 1.0;
      if (useDiag) { Rii /= get(i,i); }
      if (Double.isInfinite(Rii)) {
        throw new UnsupportedOperationException("Matrix is singular.");
      }
      result.set(i, i, Rii);
      for (int j = 0; j < i; ++j) {
        double Rij = 0.0;
        for (int k = j; k < i; ++k) { Rij -= get(i,k) * result.get(k,j); }
        result.set(i, j, Rii * Rij);
      }
    }
    return result;
  }

  /**
   * Multiplication of a upper triangular (U) and a unit lower triangular
   * matrix (L) as U*L, both stored in <code>this</code> square matrix as
   * <code>[L \ U]</code>.
   */
  private Matrix mulUL(Matrix result) {
    final int n = rows(); // == cols()
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        int k = Math.max(i, j);
        double value = 0.0;
        if (k == j) {
          value = get(i,j);
          ++k;
        }
        while (k < n) {
          value += get(i,k) * get(k,j);
          ++k;
        }
        result.set(i, j, value);
      }
    }
    return result;
  }

  /**
   * Computes the inverse of a square matrix (in new matrix).
   *
   * @return the matrix inverse
   * @throws UnsupportedOperationException if the matrix is singular
   */
  public Matrix inv() {
    return inv(create(rows(), cols()),
               create(rows(), cols()),
               Permutation.eye(rows()));
  }

  /**
   * Computes the inverse of a (symmatric) positive-definite matrix
   * (in <code>result</code>). Symmetry is not verified, violating this condition
   * might silently produce an invalid result.
   *
   * Matrix <code>result</code> and <code>tmpM</code> should be square matrices
   * of the same size as <code>this</code>. They should be different from each
   * other and from <code>this</code>.
   *
   * @param result storage of the result
   *               (not <code>null</code> and not equal to
   *                <code>tmpM</code> and <code>this</code>)
   * @param tmpM temporary matrix (not <code>null</code> and not equal to
   *                               <code>tmpM</code> and <code>this</code>,
   *                               having size of <code>rows() x cols()</code>)
   * @return <code>result</code> holding the matrix inverse
   * @throws UnsupportedOperationException if the matrix is singular
   */
  public Matrix invPsd(Matrix result, Matrix tmpM) {
    assert (rows() == cols());
    assert (result != null && result != this && result != tmpM &&
            result.rows() == rows() && result.cols() == cols());
    assert (tmpM != null && tmpM != this &&
            tmpM.rows() == rows() && tmpM.cols() == cols());
    choleskyL(result);
    tmpM.setToZeros();
    result.invLT(tmpM, true);
    return tmpM.T().mul(tmpM, result); // TODO: optimize (private method mulLL).
  }

  /**
   * Computes the inverse of a (symmatric) positive-definite matrix
   * (in new matrix). Symmetry is not verified, violating this condition
   * might silently produce an invalid result.
   *
   * @return the matrix inverse
   * @throws UnsupportedOperationException if the matrix is singular
   */
  public Matrix invPsd() {
    return invPsd(create(rows(), cols()), create(rows(), cols()));
  }

  //----------------------------------------------------------------------------

  /**
   * Returns the matrix transpose (lazy operation).
   *
   * Data is not moved, but shared between <code>this</code> and the returned
   * matrix object. Only the access to the data is changed by the redefinition
   * of the get/set methods. Use <code>copy().T()</code> to avoid data sharing.
   *
   * @return transpose of the matrix (new object with shared data)
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
