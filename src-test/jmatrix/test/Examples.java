package jmatrix.test;

import jmatrix.Matrix;
import static jmatrix.Matrix.NR;
import static jmatrix.Matrix.TOL;

import org.junit.Test;
import static jmatrix.test.MatrixAssert.assertMatrixEquals;

public class Examples {

  @Test public void solveEqnFR() {
    // A has full rank.
    Matrix A = Matrix.create(2., -3., -1., 2., NR,
                             4., -5., -1., 4., NR,
                             2., -5., -2., 2., NR,
                             0.,  2.,  1., 3.);
    Matrix b = Matrix.create(4., 4., 9., -5.).T();
    Matrix[] LU = A.LU();
    Matrix x = LU[1].backsU(LU[0].backsL(LU[2].mul(b)));
    assertMatrixEquals(Matrix.create(-1., -1., -3., 0.).T(), x, TOL);
  }

  @Test public void solveLSFCR() {
    // A has full column rank.
    Matrix A = Matrix.create(1., 1., NR, 1., 2., NR, 1., 3., NR, 1., 4.);
    Matrix b = Matrix.create(6., 5., 7., 10.).T();

    // by QR decomposition
    Matrix[] QR = A.QR(false); // without computing Q
    Matrix xQR = QR[1].backsU(QR[1].T().backsL(A.T().mul(b)));
    assertMatrixEquals(Matrix.create(3.5, 1.4).T(), xQR, TOL);

    // by Cholesky decomposition
    Matrix L = A.T().mul(A).choleskyL();
    Matrix xLL = L.T().backsU(L.backsL(A.T().mul(b)));
    assertMatrixEquals(Matrix.create(3.5, 1.4).T(), xLL, TOL);
  }

  // @Test public void solveLeastSquaresFR() {
  //   // A has full rank.
  //   Matrix A = Matrix.create(2., -3., -1., 2., NR,
  //                            4., -5., -1., 4., NR,
  //                            2., -5., -2., 2., NR,
  //                            0.,  2.,  1., 3.);
  //   Matrix b = Matrix.create(4., 4., 9., -5.).T();
  //   //Matrix x = A.solve(b);
  //   Matrix[] QR = A.QR(false); // without computing Q
  //   Matrix x = QR[1].backsU(QR[1].T().backsL(A.T().mul(b)));
  //   assertMatrixEquals(b, A.mul(x), TOL);
  // }

  // @Test public void solveLeastSquaresFCR() {
  //   // A has full column rank.
  //   Matrix A = Matrix.create(1., 1., NR,
  //                            1., 2., NR,
  //                            1., 3., NR,
  //                            1., 4.);
  //   Matrix b = Matrix.create(6., 5., 7., 10.).T();
  //   // Matrix x = A.solve(b);
  //   Matrix[] QR = A.QR(false); // without computing Q
  //   Matrix x = QR[1].backsU(QR[1].T().backsL(A.T().mul(b)));
  //   assertMatrixEquals(Matrix.create(3.5, 1.4).T(), x, TOL);
  // }

  // @Test public void solveLeastSquaresFRR() {
  //   // A has full row rank.
  //   Matrix A = Matrix.create(1., 1., NR,
  //                            1., 2., NR,
  //                            1., 3., NR,
  //                            1., 4.).T();
  //   Matrix b = Matrix.create(6., 5.).T();
  //   // Matrix x = A.solve(b);
  //   Matrix[] QR = A.T().reducedQR();
  //   Matrix x = QR[0].mul(QR[1].T().backsL(b));
  //   assertMatrixEquals(Matrix.create(4.5, 2.5, 0.5, -1.5).T(), x, TOL);
  // }

  // @Test public void solveLeastSquaresLCR() {
  //   // A has low (column) rank.
  //   Matrix A = Matrix.create(1., 1., 3., NR,
  //                            1., 2., 5., NR,
  //                            1., 3., 7., NR,
  //                            1., 4., 9.);
  //   Matrix b = Matrix.create(6., 5., 7., 10.).T();
  //   Matrix[] QR = A.reducedQR();
  //   Matrix x = QR[1].backsU(QR[0].T().mul(b));
  //   assertMatrixEquals(Matrix.create(2.45, -0.7, 1.05).T(), x, TOL);
  // }

  // @Test public void solveLeastSquaresLRR() {
  //   // A has low (row) rank.
  //   Matrix A = Matrix.create(1., 1., 1., 1., NR,
  //                            1., 2., 3., 4., NR,
  //                            3., 5., 7., 9.);
  //   Matrix b = Matrix.create(6., 5., 7.).T();
  //   Matrix[] QR = A.T().reducedQR();
  //   Matrix x = QR[0].mul(QR[1].T().backsL(b));
  //   assertMatrixEquals(Matrix.create(3.9, 2.05, 0.2, -1.65).T(), x, TOL);
  // }
}
