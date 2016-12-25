package jmatrix;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static jmatrix.MatrixAssert.assertMatrixEquals;

import static jmatrix.Matrix.NR;
import static jmatrix.Matrix.TOL;

public class Examples {

  @Test public void solveEquationLU() {
    Matrix A = Matrix.create(2., -3., -1., 2., NR,
                             4., -4., -1., 4., NR,
                             2., -5., -2., 2., NR,
                             0.,  2.,  1., 0.);
    Matrix b = Matrix.create(4., 4., 9., -5.).T();
    Matrix[] LU = A.LU();
    Matrix x = LU[0].backsL(LU[1].backsU(b), true);
    Matrix r = A.mul(x).sub(b);
    assertEquals(0.0, r.normI(), TOL);
  }

  @Test public void solveLeastSquaresFR() { // full rank
    Matrix A = Matrix.create(2., -3., -1., 2., NR,
                             4., -4., -1., 4., NR,
                             2., -5., -2., 2., NR,
                             0.,  2.,  1., 0.);
    Matrix b = Matrix.create(4., 4., 9., -5.).T();
    Matrix[] QR = A.QR(false); // without computing Q
    Matrix x = QR[1].backsU(QR[1].T().backsL(A.T().mul(b)));
    assertEquals(0.0, A.mul(x).sub(b).normI(), TOL);
  }

  @Test public void solveLeastSquaresFCR() { // full column rank
    Matrix A = Matrix.create(1., 1., NR,
                             1., 2., NR,
                             1., 3., NR,
                             1., 4.);
    Matrix b = Matrix.create(6., 5., 7., 10.).T();
    Matrix[] LUP = A.T().mul(A).LU();
    Matrix x = LUP[1].backsU(LUP[0].backsL(LUP[2].mul(A.T().mul(b))));
    assertMatrixEquals(Matrix.create(3.5, 1.4).T(), x, TOL);
  }
}
