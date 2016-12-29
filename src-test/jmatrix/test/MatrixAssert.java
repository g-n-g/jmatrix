package jmatrix.test;

import jmatrix.Matrix;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

/**
 * Matrix assertions.
 */
public class MatrixAssert {

  /** Asserts that two matrices have the same size. */
  public static void assertMatrixSizeEquals(Matrix A, Matrix B) {
    assertEquals(A.rows(), B.rows());
    assertEquals(A.cols(), B.cols());
  }

  /** Asserts that matrix is square. */
  public static void assertMatrixSquareSize(Matrix A) {
    assertEquals(A.rows(), A.cols());
  }

  /** Asserts that matrix has no NaN values. */
  public static void assertMatrixHasNoNaNs(Matrix A) {
    for (int i = 0; i < A.rows(); ++i) {
      for (int j = 0; j < A.cols(); ++j) {
        assertFalse(Double.isNaN(A.get(i,j)));
      }
    }
  }

  /**
   * Asserts that two matrices are equal.
   */
  public static void assertMatrixEquals(Matrix A, Matrix B, double tol) {
    assertNotNull(A);
    assertNotNull(B);
    assertMatrixSizeEquals(A, B);
    for (int i = 0; i < A.rows(); ++i) {
      for (int j = 0; j < A.cols(); ++j) {
        assertEquals(A.get(i,j), B.get(i,j), tol);
      }
    }
  }

  /** Asserts that matrix is lower triangular. */
  public static void assertMatrixLT(Matrix A, double tol) {
    for (int i = 0; i < A.rows(); ++i) {
      for (int j = i+1; j < A.cols(); ++j) {
        assertEquals(0.0, A.get(i,j), tol);
      }
    }
  }

  /** Asserts that matrix is unit lower triangular. */
  public static void assertMatrixUnitLT(Matrix A, double tol) {
    assertMatrixLT(A, tol);
    for (int i = 0; i < Math.min(A.rows(), A.cols()); ++i) {
      assertEquals(1.0, A.get(i,i), tol);
    }
  }

  /** Asserts that matrix is upper triangular. */
  public static void assertMatrixUT(Matrix A, double tol) {
    assertMatrixLT(A.T(), tol);
  }

  /** Asserts that matrix is unit upper triangular. */
  public static void assertMatrixUnitUT(Matrix A, double tol) {
    assertMatrixUnitLT(A.T(), tol);
  }

  /** Asserts that matrix is diagonal. */
  public static void assertMatrixDiag(Matrix A, double tol) {
    assertMatrixLT(A, tol);
    assertMatrixUT(A, tol);
  }

  /** Asserts that matrix is zero. */
  public static void assertMatrixZero(Matrix A, double tol) {
    for (int i = 0; i < A.rows(); ++i) {
      for (int j = 0; j < A.cols(); ++j) {
        assertEquals(0.0, A.get(i,j), tol);
      }
    }
  }

  /** Asserts that matrix is identity. */
  public static void assertMatrixEye(Matrix A, double tol) {
    assertMatrixSquareSize(A);
    assertMatrixUnitLT(A, tol);
    assertMatrixUT(A, tol);
  }

  /** Checks that matrix is diagonal with zero or one elements only. */
  public static void assertMatrixEyeWithZeros(Matrix A, double tol) {
    assertMatrixSquareSize(A);
    assertMatrixDiag(A, tol);
    for (int i = 0; i < A.rows(); ++i) {
      double Aii = A.get(i,i);
      double diff = Math.min(Math.abs(Aii), Math.abs(1.0-Aii));
      assertEquals(0.0, diff, tol);
    }
  }

  /** Asserts that matrix has orthogonal rows. */
  public static void assertMatrixOrthoRows(Matrix A, double tol) {
    assertMatrixEye(A.mul(A.T()), tol);
  }

  /** Asserts that matrix has orthogonal columns. */
  public static void assertMatrixOrthoCols(Matrix A, double tol) {
    assertMatrixEye(A.T().mul(A), tol);
  }

  /** Asserts that matrix is orthogonal. */
  public static void assertMatrixOrtho(Matrix A, double tol) {
    assertMatrixOrthoRows(A, tol);
    assertMatrixOrthoCols(A, tol);
  }

  /** Asserts that matrix has orthogonal or zero rows. */
  public static void assertMatrixOrthoOrZeroRows(Matrix A, double tol) {
    assertMatrixEyeWithZeros(A.mul(A.T()), tol);
  }

  /** Asserts that matrix has orthogonal or zero columns. */
  public static void assertMatrixOrthoOrZeroCols(Matrix A, double tol) {
    assertMatrixEyeWithZeros(A.T().mul(A), tol);
  }

  /** Asserts that matrix has orthogonal or zero rows and columns. */
  public static void assertMatrixOrthoOrZero(Matrix A, double tol) {
    assertMatrixOrthoOrZeroRows(A, tol);
    assertMatrixOrthoOrZeroCols(A, tol);
  }

  /** Asserts that matrix is symmetric. */
  public static void assertMatrixSymmetric(Matrix A, double tol) {
    assertMatrixSquareSize(A);
    for (int i = 0; i < A.rows(); ++i) {
      for (int j = 0; j < i; ++j) {
        assertEquals(A.get(i,j), A.get(j,i), tol);
      }
    }
  }
}
