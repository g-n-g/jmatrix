package jmatrix;

import static org.junit.Assert.assertEquals;

/**
 * Matrix assertions.
 */
public class MatrixAssert {

  /**
   * Asserts that two matrices are equal.
   */
  public static void assertMatrixEquals(Matrix a, Matrix b, double tol) {
    assertEquals(a.rows(), b.rows());
    assertEquals(a.cols(), b.cols());
    for (int i = 0; i < a.rows(); ++i) {
      for (int j = 0; j < a.cols(); ++j) {
        assertEquals(a.get(i,j), b.get(i,j), tol);
      }
    }
  }

  /**
   * Asserts that matrix is lower triangular.
   */
  public static void assertMatrixLT(Matrix m, double tol) {
    for (int i = 0; i < m.rows(); ++i) {
      for (int j = i+1; j < m.cols(); ++j) {
        assertEquals(0.0, m.get(i,j), tol);
      }
    }
  }

  /**
   * Asserts that matrix is unit lower triangular.
   */
  public static void assertMatrixUnitLT(Matrix m, double tol) {
    assertMatrixLT(m, tol);
    for (int i = 0; i < Math.min(m.rows(), m.cols()); ++i) {
      assertEquals(1.0, m.get(i,i), tol);
    }
  }

  /**
   * Asserts that matrix is upper triangular.
   */
  public static void assertMatrixUT(Matrix m, double tol) {
    assertMatrixLT(m.T(), tol);
  }

  /**
   * Asserts that matrix is unit upper triangular.
   */
  public static void assertMatrixUnitUT(Matrix m, double tol) {
    assertMatrixUnitLT(m.T(), tol);
  }

  /**
   * Asserts that matrix is diagonal.
   */
  public static void assertMatrixDiag(Matrix m, double tol) {
    assertMatrixLT(m, tol);
    assertMatrixUT(m, tol);
  }

  /**
   * Asserts that matrix is orthogonal.
   */
  public static void assertMatrixOrtho(Matrix m, double tol) {
    assertEquals(m.rows(), m.cols());
    Matrix I = Matrix.eye(m.rows());
    assertMatrixEquals(I, m.mul(m.T()), tol);
    assertMatrixEquals(I, m.T().mul(m), tol);
  }

  /**
   * Asserts that matrix has orthogonal rows.
   */
  public static void assertMatrixOrthoRows(Matrix m, double tol) {
    Matrix I = Matrix.eye(m.rows());
    assertMatrixEquals(I, m.mul(m.T()), tol);
  }

  /**
   * Asserts that matrix has orthogonal columns.
   */
  public static void assertMatrixOrthoCols(Matrix m, double tol) {
    Matrix I = Matrix.eye(m.cols());
    assertMatrixEquals(I, m.T().mul(m), tol);
  }
}
