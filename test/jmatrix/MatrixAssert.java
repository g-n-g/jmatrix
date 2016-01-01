package jmatrix;

import static org.junit.Assert.assertEquals;

/**
 * Matrix assertions.
 */
public class MatrixAssert {

  /**
   * Asserts that two matrices are equal.
   */
  public static double assertMatrixEquals(Matrix A, Matrix B, double tol) {
    final int rows = A.rows(), cols = A.cols();
    assertEquals(rows, B.rows());
    assertEquals(cols, B.cols());
    double err = 0.0;
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        double Aij = A.get(i,j);
        double Bij = B.get(i,j);
        double diff = Math.abs(Aij - Bij);
        if (diff > err) { err = diff; }
        assertEquals(Aij, Bij, tol);
      }
    }
    return err;
  }

  /**
   * Asserts that matrix is lower triangular.
   */
  public static double assertMatrixLT(Matrix M, double tol) {
    final int rows = M.rows(), cols = M.cols();
    double err = 0.0;
    for (int i = 0; i < rows; ++i) {
      for (int j = i+1; j < cols; ++j) {
        double Mij = M.get(i,j);
        double diff = Math.abs(Mij);
        if (diff > err) { err = diff; }
        assertEquals(0.0, Mij, tol);
      }
    }
    return err;
  }

  /**
   * Asserts that matrix is unit lower triangular.
   */
  public static double assertMatrixUnitLT(Matrix M, double tol) {
    double err = assertMatrixLT(M, tol);
    for (int i = 0; i < Math.min(M.rows(), M.cols()); ++i) {
      double Mii = M.get(i,i);
      double diff = Math.abs(1.0 - Mii);
      if (diff > err) { err = diff; }
      assertEquals(1.0, Mii, tol);
    }
    return err;
  }

  /**
   * Asserts that matrix is upper triangular.
   */
  public static double assertMatrixUT(Matrix M, double tol) {
    return assertMatrixLT(M.T(), tol);
  }

  /**
   * Asserts that matrix is unit upper triangular.
   */
  public static double assertMatrixUnitUT(Matrix M, double tol) {
    return assertMatrixUnitLT(M.T(), tol);
  }

  /**
   * Asserts that matrix is diagonal.
   */
  public static double assertMatrixDiag(Matrix M, double tol) {
    double errLT = assertMatrixLT(M, tol);
    double errUT = assertMatrixUT(M, tol);
    return Math.max(errLT, errUT);
  }

  /**
   * Asserts that matrix is identity (might allow zeros on the diagonal).
   */
  public static double assertMatrixEye(Matrix M, boolean allowZero, double tol) {
    assertEquals(M.rows(), M.cols());
    double err = assertMatrixDiag(M, tol);
    for (int i = 0; i < M.rows(); ++i) {
      double Mii = M.get(i,i);
      double diff = Math.abs(1.0 - Mii);
      if (allowZero) {
        diff = Math.min(diff, Mii);
        assertEquals(0.0, diff, tol);
      }
      else {
        assertEquals(1.0, Mii, tol);
      }
      if (diff > err) { err = diff; }
    }
    return err;
  }

  /**
   * Asserts that matrix is identity.
   */
  public static double assertMatrixEye(Matrix M, double tol) {
    return assertMatrixEye(M, false, tol);
  }

  /**
   * Asserts that matrix is orthogonal (might allow zeros on the diagonal).
   */
  public static double assertMatrixOrtho(Matrix M, boolean allowZero, double tol) {
    assertEquals(M.rows(), M.cols());
    double errMMt = assertMatrixEye(M.mul(M.T()), allowZero, tol);
    double errMtM = assertMatrixEye(M.T().mul(M), allowZero, tol);
    return Math.max(errMMt, errMtM);
  }

  /**
   * Asserts that matrix is orthogonal.
   */
  public static double assertMatrixOrtho(Matrix M, double tol) {
    return assertMatrixOrtho(M, false, tol);
  }

  /**
   * Asserts that matrix has orthogonal rows (might allow zeros on the diagonal).
   */
  public static double assertMatrixOrthoRows(Matrix M, boolean allowZero, double tol) {
    return assertMatrixEye(M.mul(M.T()), allowZero, tol);
  }

  /**
   * Asserts that matrix has orthogonal rows.
   */
  public static double assertMatrixOrthoRows(Matrix M, double tol) {
    return assertMatrixOrthoRows(M, false, tol);
  }

  /**
   * Asserts that matrix has orthogonal columns (might allow zeros on the diagonal).
   */
  public static double assertMatrixOrthoCols(Matrix M, boolean allowZero, double tol) {
    return assertMatrixEye(M.T().mul(M), allowZero, tol);
  }

  /**
   * Asserts that matrix has orthogonal columns.
   */
  public static double assertMatrixOrthoCols(Matrix M, double tol) {
    return assertMatrixOrthoCols(M, false, tol);
  }

  /**
   * Asserts that matrix is symmetric.
   */
  public static double assertMatrixSymmetric(Matrix M, double tol) {
    assertEquals(M.rows(), M.cols());
    double err = 0.0;
    for (int i = 0; i < M.rows(); ++i) {
      for (int j = 0; j < i; ++j) {
        double Mij = M.get(i,j);
        double Mji = M.get(j,i);
        double diff = Math.abs(Mij - Mji);
        if (diff > err) { err = diff; }
        assertEquals(Mij, Mji, tol);
      }
    }
    return err;
  }
}
