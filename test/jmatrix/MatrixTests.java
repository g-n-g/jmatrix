package jmatrix;

import java.util.Random;

import org.junit.Test;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertEquals;
import static jmatrix.MatrixAssert.assertMatrixEquals;
import static jmatrix.MatrixAssert.assertMatrixLT;
import static jmatrix.MatrixAssert.assertMatrixUnitLT;
import static jmatrix.MatrixAssert.assertMatrixUT;
import static jmatrix.MatrixAssert.assertMatrixOrtho;
import static jmatrix.MatrixAssert.assertMatrixOrthoCols;

import static jmatrix.Matrix.NR;
import static jmatrix.Matrix.TOL;
import static jmatrix.BasicUnaryOperation.*;
import static jmatrix.BasicBinaryOperation.*;

/**
 * Tests for the Matrix class.
 */
public class MatrixTests {

  public static final double PREC = 1e-8;
  public static final Random RNG = new Random();

  //---------------------------------------------------------------------------

  @Test public void isEmpty() {
    Matrix e = Matrix.create(0,0);
    assertTrue(e.isEmpty());
    assertEquals(0, e.rows());
    assertEquals(0, e.cols());

    e = Matrix.create(0,1);
    assertTrue(e.isEmpty());
    assertEquals(0, e.rows());
    assertEquals(1, e.cols());

    e = Matrix.create(1,0);
    assertTrue(e.isEmpty());
    assertEquals(1, e.rows());
    assertEquals(0, e.cols());

    e = Matrix.create(1,1);
    assertFalse(e.isEmpty());
    assertEquals(1, e.rows());
    assertEquals(1, e.cols());
  }

  @Test public void norms() {
    Matrix m = Matrix.create(1.0, 1.5, -0.5, 1.5, NR,
                             2.0, 4.0, -4.0, 0.0, NR,
                             0.5, 3.0, -3.0, 2.0);
    assertEquals(8.5, m.norm1(), PREC);
    assertEquals(10.0, m.normI(), PREC);
    assertEquals(8.0, m.normF(), PREC);

    m = Matrix.create(new double[][]{ // test this form of create
        new double[]{1.0, 1.5, -0.5, 1.5},
        new double[]{2.0, 4.0, -4.0, 0.0},
        new double[]{0.5, 3.0, -3.0, 2.0}
      });
    assertEquals(8.5, m.norm1(), PREC);
    assertEquals(10.0, m.normI(), PREC);
    assertEquals(8.0, m.normF(), PREC);
  }

  @Test public void checkNaNandInf() {
    Matrix m = Matrix.create(1.0, 1.5, -0.5, 1.5, NR,
                             2.0, 4.0, -4.0, 0.0, NR, NR, NR,
                             0.5, 3.0, -3.0, 2.0, NR);
    assertFalse(m.hasNaN());
    assertFalse(m.hasInf());

    m.set(1, 1, Double.NaN);
    assertTrue(m.hasNaN());
    assertFalse(m.hasInf());

    m.set(1, 2, Double.POSITIVE_INFINITY);
    assertTrue(m.hasNaN());
    assertTrue(m.hasInf());

    m.set(1, 1, 4.4);
    m.set(1, 2, Double.NEGATIVE_INFINITY);
    assertFalse(m.hasNaN());
    assertTrue(m.hasInf());

    m.set(1, 1, Double.NaN);
    m.set(2, 2, Double.POSITIVE_INFINITY);
    m.replaceNaNandInf(1.1, 2.2, 3.3);
    assertFalse(m.hasNaN());
    assertFalse(m.hasInf());
  }

  @Test public void basicLinearOps() {
    final int n = 8, m = 17;
    Matrix m1 = Matrix.zeros(n, m);
    Matrix m2 = Matrix.ones(n, m);
    Matrix m3 = Matrix.scalars(n, m, 42.42);
    Matrix m4 = Matrix.eye(n);

    assertMatrixEquals(m2, m2.add(m1), PREC);
    assertMatrixEquals(m3, m2.mul(42.42), PREC);
    assertMatrixEquals(m2, m3.div(42.42), PREC);

    assertMatrixEquals(m2, m1.add(1.0), PREC);
    assertMatrixEquals(m1, m2.sub(1.0), PREC);

    m2.add(m2, m2);
    assertMatrixEquals(m2, m3.div(42.42).mul(2), PREC);
    m2.setToOnes();

    m2.sub(m2.div(2), m2);
    assertMatrixEquals(m2, m3.div(42.42).mul(0.5), PREC);
    m2.setToOnes();

    m2.mul(2, m2);
    assertMatrixEquals(m2, m3.div(42.42).mul(2), PREC);
    m2.setToOnes();

    m2.div(8, m2);
    assertMatrixEquals(m2, m3.div(8*42.42), PREC);
    m2.setToOnes();
  }

  @Test public void ewu() {
    // ABS, NEG, SIGN
    Matrix m = Matrix.randN(4, 6, RNG);
    assertMatrixEquals(m, m.ewu(SIGN).emul(m.ewu(ABS)), PREC);
    assertMatrixEquals(m.ewu(NEG), m.ewu(NEG).ewu(SIGN).emul(m.ewu(ABS)), PREC);

    m.ewu(NEG);
    assertMatrixEquals(m, m.ewu(SIGN).emul(m.ewu(ABS)), PREC);

    Matrix s = Matrix.zeros(4, 6);
    Matrix a = Matrix.zeros(4, 6);
    m.copy(s).ewu(SIGN, s);
    m.copy(a).ewu(ABS, m);
    assertMatrixEquals(m, s.emul(a), PREC);

    s.copy(a).ewu(SIGN, s);
    assertMatrixEquals(s, a, PREC);

    s.set(2, 2, 0.0);
    a.set(2, 2, 0.0);
    assertMatrixEquals(s, a, PREC);
  }

  @Test public void ewb12() {
    // RECIPROC
    int nr = 2, nc = 5;
    Matrix v = Matrix.scalars(nr, nc, 0.1);
    Matrix o = Matrix.ones(nr, nc);
    Matrix r = Matrix.zeros(nr, nc);

    assertMatrixEquals(o, v.emul(v.ewb1(DIV, 1.0)), PREC);

    v.ewb1(DIV, 1.0, r);
    assertMatrixEquals(o, v.emul(r), PREC);

    r.ewb1(DIV, 1.0, r);
    assertMatrixEquals(v, r, PREC);

    // POW
    Matrix m1 = Matrix.create(1, 2, 3, NR,
                              4, 5, 6);
    Matrix m2 = Matrix.create(1, 4, 9, NR,
                              16, 25, 36);
    assertMatrixEquals(m2, m1.ewb2(POW, 2.0), PREC);

    Matrix mr = Matrix.randN(4, 5, RNG);
    assertMatrixEquals(mr.emul(mr), mr.ewb2(POW, 2.0), PREC);
    assertMatrixEquals(mr.emul(mr).emul(mr), mr.ewb2(POW, 3.0), PREC);

    // MOD
    double m = 2.3;
    Matrix oo = Matrix.scalars(5, 4, m);
    assertMatrixEquals(Matrix.zeros(5, 4), oo.ewb2(MOD, m), PREC);

    oo.ewb2(MOD, m, oo);
    assertMatrixEquals(Matrix.zeros(5, 4), oo.ewb2(MOD, m), PREC);

    Matrix mm = Matrix.rand(5, 4, RNG);
    mm.mul(8.0, mm);
    assertMatrixEquals(mm.add(oo.mul(2.0)).ewb2(MOD, m), mm.ewb2(MOD, m), PREC);

    mm.mul(-1.0, mm);
    assertMatrixEquals(mm.add(oo.mul(-2.0)).ewb2(MOD, m), mm.ewb2(MOD, m), PREC);
  }

  @Test public void testewb() {
    Matrix m = Matrix.create(1, 2, 3, NR,
                             4, 5, 6);

    Matrix r = Matrix.create(1, 2, 3);
    Matrix t1 = Matrix.create(1,  4,  9, NR,
                              4, 10, 18);
    assertMatrixEquals(t1, m.ewb(MUL, r), PREC);
    m.ewb(MUL, r, m);
    assertMatrixEquals(t1, m, PREC);

    m = Matrix.create(1, 2, 3, NR,
                      4, 5, 6);

    Matrix c = Matrix.create(1.0, 2.0).T();
    Matrix t2 = Matrix.create(1,  2,  3, NR,
                              8, 10, 12);
    assertMatrixEquals(t2, m.ewb(MUL, c), PREC);
    m.ewb(MUL, c, m);
    assertMatrixEquals(t2, m, PREC);
  }

  @Test public void copy() {
    Matrix m = Matrix.rand(5, 7, RNG);

    assertMatrixEquals(m, m.copy(m), PREC);

    Matrix c = m.copy();
    assertMatrixEquals(m, c, PREC);

    m.add(Matrix.ones(5, 7), m);
    c.copy(m);
    assertTrue(PREC > m.sub(c).norm1());
  }

  @Test public void blkdiag() {
    Matrix m1 = Matrix.create(1.1, 2.2);
    Matrix m2 = Matrix.create(3.3, 4.4, 5.5, NR,
                              6.6, 7.7, 8.8);
    Matrix m3 = Matrix.create(9.9);

    Matrix m = Matrix.blkdiag(m1, m2, m3);
    Matrix t = Matrix.create(1.1, 2.2, 0.0, 0.0, 0.0, 0.0, NR,
                             0.0, 0.0, 3.3, 4.4, 5.5, 0.0, NR,
                             0.0, 0.0, 6.6, 7.7, 8.8, 0.0, NR,
                             0.0, 0.0, 0.0, 0.0, 0.0, 9.9);
    assertMatrixEquals(t, m, PREC);
  }

  @Test public void cat() {
    Matrix m1 = Matrix.create(0.0, 1.1, NR,
                              0.0, 2.2);
    Matrix m2 = Matrix.create(3.3, 4.4, 5.5, NR,
                              6.6, 7.7, 8.8);
    Matrix m3 = Matrix.create(9.9, NR,
                              9.9);

    Matrix m = Matrix.horzcat(m1, m2, m3);
    Matrix t = Matrix.create(0.0, 1.1, 3.3, 4.4, 5.5, 9.9, NR,
                             0.0, 2.2, 6.6, 7.7, 8.8, 9.9);
    assertMatrixEquals(t, m, PREC);

    m1 = m1.T();
    m2 = m2.T();
    m3 = m3.T();
    m = Matrix.vertcat(m1, m2, m3);
    t = t.T();
    assertMatrixEquals(t, m, PREC);
  }

  @Test public void setTo() {
    Matrix m = Matrix.ones(4, 7);
    m.setToZeros();
    assertMatrixEquals(Matrix.zeros(4, 7), m, PREC);
    m.setToOnes();
    assertMatrixEquals(Matrix.ones(4, 7), m, PREC);

    Matrix I = Matrix.eye(4);
    I.setToScalars(2);
    assertMatrixEquals(Matrix.scalars(4, 4, 2), I, PREC);
    I.setToEye();
    assertMatrixEquals(Matrix.eye(4), I, PREC);
  }

  @Test public void getSetMat() {
    Matrix A = Matrix.rand(7, 7, RNG);
    Matrix B = Matrix.rand(7, 7, RNG);

    Matrix AB = Matrix.scalars(10, 10, Double.NaN);
    AB.setMat(2, 8, 2, 8, A.mul(B));

    assertMatrixEquals(AB.getMat(2, 2, 4, 5),
                       A.getMat(0, 0, 0, 6).mul(B.getMat(0, 6, 2, 3)),
                       PREC);
    assertMatrixEquals(AB.getMat(6, 8, 4, 5),
                       A.getMat(4, 6, 0, 6).mul(B.getMat(0, 6, 2, 3)),
                       PREC);
    assertMatrixEquals(AB.getMat(3, 4, 4, 5),
                       A.getMat(1, 2, 0, 6).mul(B.getMat(0, 6, 2, 3)),
                       PREC);
  }

  @Test public void getDiag() {
    Matrix A = Matrix.create(1, 2, 3, NR,
                             4, 5, 6, NR,
                             7, 8, 9);
    Matrix a = Matrix.create(1, 5, 9).T();
    assertEquals(0.0, A.getDiag().sub(a).normI(), PREC);
    assertEquals(0.0, A.T().getDiag().sub(a).normI(), PREC);

    Matrix B = Matrix.create(1, 2, 3, 4, 5, NR,
                             6, 7, 8, 9, 0);
    Matrix b = Matrix.create(1.0, 7.0).T();
    assertEquals(0.0, B.getDiag().sub(b).normI(), PREC);
    assertEquals(0.0, B.T().getDiag().sub(b).normI(), PREC);
  }

  @Test public void getTrilLU() {
    Matrix A = Matrix.create(1, 2, 3, NR,
                             4, 5, 6, NR,
                             7, 8, 9);

    Matrix L1 = Matrix.create(1, 0, 0, NR,
                              4, 5, 0, NR,
                              7, 8, 9);
    assertEquals(0.0, A.getTriL().sub(L1).normF(), PREC);

    Matrix L2 = Matrix.create(0, 0, 0, NR,
                              4, 0, 0, NR,
                              7, 8, 0);
    assertEquals(0.0, A.getTriL(-1).sub(L2).normF(), PREC);

    Matrix L3 = Matrix.create(1, 2, 0, NR,
                              4, 5, 6, NR,
                              7, 8, 9);
    assertEquals(0.0, A.getTriL(1).sub(L3).normF(), PREC);

    assertEquals(0.0, A.getTriL(-3).normF(), PREC);
    assertEquals(0.0, A.getTriL(-5).normF(), PREC);
    assertEquals(0.0, A.getTriL(2).sub(A).normF(), PREC);
    assertEquals(0.0, A.getTriL(5).sub(A).normF(), PREC);

    Matrix U1 = Matrix.create(1, 2, 3, NR,
                              0, 5, 6, NR,
                              0, 0, 9);
    assertEquals(0.0, A.getTriU().sub(U1).normF(), PREC);

    Matrix U2 = Matrix.create(1, 2, 3, NR,
                              4, 5, 6, NR,
                              0, 8, 9);
    assertEquals(0.0, A.getTriU(-1).sub(U2).normF(), PREC);

    Matrix U3 = Matrix.create(0, 2, 3, NR,
                              0, 0, 6, NR,
                              0, 0, 0);
    assertEquals(0.0, A.getTriU(1).sub(U3).normF(), PREC);

    assertEquals(0.0, A.getTriU(-2).sub(A).normF(), PREC);
    assertEquals(0.0, A.getTriU(-5).sub(A).normF(), PREC);
    assertEquals(0.0, A.getTriU(3).normF(), PREC);
    assertEquals(0.0, A.getTriU(5).normF(), PREC);
  }

  @Test public void matrixVectorProduct() {
    final int n = 5, m = 8;
    Matrix v1 = Matrix.ones(m, 1);
    Matrix v2 = Matrix.zeros(m, 1); v2.set(m/2, 0, 1.0); // unit vector
    Matrix one = Matrix.ones(n, 1);
    Matrix m1 = Matrix.ones(n, m);
    Matrix m2 = Matrix.eye(m);

    assertMatrixEquals(v1, m2.mul(v1), PREC);
    assertMatrixEquals(one, m1.mul(v2), PREC);

    for (int i = 0; i < m; ++i) v2.set(i, 0, i+1);
    double sum = m * (m+1) / 2.0;
    assertMatrixEquals(one.mul(sum), m1.mul(v2), PREC);
    assertMatrixEquals(v2, m2.mul(v2), PREC);

    assertMatrixEquals(v2.emul(v2), Matrix.diag(v2).mul(v2), PREC);
    assertMatrixEquals(v2.emul(v2), Matrix.diag(v2.T()).mul(v2), PREC);
    assertMatrixEquals(m1.mul(v1.add(v2)), m1.mul(v1).add(m1.mul(v2)), PREC);
  }

  @Test public void matrixProduct() {
    final int n = 5, m = 8;
    Matrix m1 = Matrix.ones(n, m).mul(2);
    Matrix m2 = Matrix.eye(n);
    Matrix m3 = Matrix.eye(m);

    assertMatrixEquals(m1, m1.mul(m3), PREC);
    assertMatrixEquals(m1, m2.mul(m1), PREC);

    m2.mul(3, m2);
    m3.mul(5, m3);

    assertMatrixEquals(m1.mul(m3), m3.T().mul(m1.T()).T(), PREC);
    assertMatrixEquals(m2.mul(m1), m1.T().mul(m2.T()).T(), PREC);

    Matrix m4 = Matrix.ones(n, m).mul(7);

    assertMatrixEquals(m1.add(m4).mul(m3),
                       m1.mul(m3).add(m4.mul(m3)),
                       PREC);
    assertMatrixEquals(m2.mul(m1.add(m4)),
                       m2.mul(m1).add(m2.mul(m4)),
                       PREC);

    Matrix m5 = Matrix.ones(n, m).mul(-8);

    m1.mul(m3, m5);
    assertMatrixEquals(m5, m1.mul(m3), PREC);

    m2.mul(m1, m5);
    assertMatrixEquals(m5, m2.mul(m1), PREC);
  }

  @Test public void entrywiseMultiplication() {
    final int n = 5, m = 7;
    Matrix m1 = Matrix.scalars(n, m, 3);
    Matrix m2 = Matrix.scalars(n, m, 9);
    Matrix m3 = Matrix.zeros(n, m);

    assertMatrixEquals(m2, m1.emul(m1), PREC);
    assertMatrixEquals(m3, m1.emul(m3), PREC);
    assertMatrixEquals(m3, m3.emul(m2), PREC);

    m1.emul(m2, m1);
    assertMatrixEquals(m1, m2.mul(3), PREC);
    m1.setToScalars(3.0);
  }

  @Test public void dotProduct() {
    Matrix v1 = Matrix.create(1, 2, 0, 3);
    Matrix v2 = Matrix.create(5, 2, 9, -1);

    assertEquals(6.0, v1.dot(v2), PREC);
    assertEquals(6.0, v1.T().dot(v2), PREC);
    assertEquals(6.0, v1.dot(v2.T()), PREC);
    assertEquals(6.0, v1.T().dot(v2.T()), PREC);

    assertEquals(6.0, v2.dot(v1), PREC);
    assertEquals(6.0, v2.T().dot(v1), PREC);
    assertEquals(6.0, v2.dot(v1.T()), PREC);
    assertEquals(6.0, v2.T().dot(v1.T()), PREC);
  }

  @Test public void trace() {
    Matrix m = Matrix.create(1.0, 1.5, -0.5, NR,
                             2.0, 4.0, -4.0, NR,
                             0.5, 3.0, -3.0, NR);
    assertEquals(2.0, m.trace(), PREC);
    assertEquals(-12.0, m.prodDiag(), PREC);
    assertEquals(4.0, Matrix.eye(4).trace(), PREC);
    assertEquals(8.0, Matrix.scalars(4, 4, 2).trace(), PREC);

    Matrix m3x2 = Matrix.rand(3, 2, RNG);
    Matrix m2x3 = Matrix.rand(2, 3, RNG);
    assertEquals(m3x2.mul(m2x3).trace(), m3x2.traceMul(m2x3), PREC);
  }

  @Test public void rowColNorms() {
    Matrix m = Matrix.create(17, 24, -1,  8, 15, 1, NR,
                             23, -5, -7, 14, 16, 1, NR,
                             -4,  6, 31, 20, 22, 1, NR,
                             10, 12, 19, 21,  3, 1, NR,
                             11, 18, 25, -2,  9, 1);

    Matrix n1 = Matrix.create(66, 66, 84, 66, 66).T();
    assertMatrixEquals(n1, m.rowNorms1(), PREC);
    n1 = Matrix.create(65, 65, 83, 65, 65, 5);
    assertMatrixEquals(n1, m.colNorms1(), PREC);

    Matrix n2 = Matrix.create(34.0,
                              32.4961536185438,
                              43.5660418215839,
                              32.4961536185438,
                              34.0).T();
    assertMatrixEquals(n2, m.rowNorms2(), PREC);
    n2 = Matrix.create(32.48076353782343,
                       33.24154027718932,
                       44.68780594300866,
                       33.24154027718932,
                       32.48076353782343,
                       2.23606797749979);
    assertMatrixEquals(n2, m.colNorms2(), PREC);

    Matrix nI = Matrix.create(24, 23, 31, 21, 25).T();
    assertMatrixEquals(nI, m.rowNormsI(), PREC);
    nI = Matrix.create(23, 24, 31, 21, 22, 1);
    assertMatrixEquals(nI, m.colNormsI(), PREC);
  }

  @Test public void rowColSums() {
    Matrix m = Matrix.create(17, 24,  1,  8, 15, 1, NR,
                             23,  5,  7, 14, 16, 1, NR,
                              4,  6, 31, 20, 22, 1, NR,
                             10, 12, 19, 21,  3, 1, NR,
                             11, 18, 25,  2,  9, 1);

    Matrix vRow = Matrix.create(65, 65, 83, 65, 65, 5);
    assertMatrixEquals(vRow, m.sumRows(), PREC);

    Matrix vCol = Matrix.create(66, 66, 84, 66, 66).T();
    assertMatrixEquals(vCol, m.sumCols(), PREC);
  }

  //---------------------------------------------------------------------------
  // Cholesky decomposition tests

  @Test public void CholeskyDecomposition() {
    Matrix m = Matrix.create(2.0, 1.0, 1.0, NR,
                             1.0, 2.0, 1.0, NR,
                             1.0, 1.0, 2.0);
    Matrix L = m.choleskyL();
    assertMatrixEquals(m, L.mul(L.T()), PREC);

    m.set(0, 0, 3.33333); m.set(1, 1, 4.44444); m.set(2, 2, 5.55555);
    m.choleskyL(L);
    assertMatrixEquals(m, L.mul(L.T()), PREC);
    assertMatrixLT(L, PREC);

    L = m.copy();
    L.choleskyL(L); // in-place
    assertMatrixEquals(m, L.mul(L.T()), PREC);
    assertMatrixLT(L, PREC);
  }

  @Test public void CholeskyDecompositionLD() {
    Matrix m = Matrix.create(2.0, 1.0, 1.0, NR,
                             1.0, 3.0, 1.0, NR,
                             1.0, 1.0, 2.0);
    Matrix LD[] = m.choleskyLD();
    Matrix L = LD[0];
    Matrix D = LD[1];
    assertMatrixEquals(m, L.mul(D).mul(L.T()), PREC);

    L = m.copy();
    D = Matrix.create(3, 1);
    L.choleskyLD(L, D); // "in-place"
    assertMatrixEquals(m, L.mul(Matrix.diag(D)).mul(L.T()), PREC);
  }

  //---------------------------------------------------------------------------
  // QR decomposition tests

  @Test public void basicQR() {
    Matrix M1 = Matrix.zeros(2, 2);
    Matrix[] QR = M1.QR();
    Matrix Q = QR[0];
    Matrix R = QR[1];
    assertMatrixEquals(M1, Q.mul(R), PREC);
    assertMatrixOrtho(Q, PREC);
    assertMatrixUT(R, PREC);

    Matrix M2 = Matrix.eye(2);
    QR = M1.QR();
    Q = QR[0];
    R = QR[1];
    assertMatrixEquals(M1, Q.mul(R), PREC);
    assertMatrixOrtho(Q, PREC);
    assertMatrixUT(R, PREC);

    Matrix M3 = Matrix.create(-1.0, 2.0, NR,
                              2.0, -1.0);
    QR = M3.QR();
    Q = QR[0];
    R = QR[1];
    assertMatrixEquals(M3, Q.mul(R), PREC);
    assertMatrixOrtho(Q, PREC);
    assertMatrixUT(R, PREC);
  }

  @Test public void fullrankQR() {
    Matrix M1 = Matrix.create(1, 2, NR,
                              4, 1);
    Matrix[] QR = M1.QR();
    Matrix Q = QR[0];
    Matrix R = QR[1];
    assertMatrixEquals(M1, Q.mul(R), PREC);
    assertMatrixOrtho(Q, PREC);
    assertMatrixUT(R, PREC);

    Matrix M2 = Matrix.create(17, 24,  1,  8, 15, NR,
                              23,  5,  7, 14, 16, NR,
                               4,  6, 13, 20, 22, NR,
                              10, 12, 19, 21,  3, NR,
                              11, 18, 25,  2,  9);
    QR = M2.QR();
    Q = QR[0];
    R = QR[1];
    assertMatrixEquals(M2, Q.mul(R), PREC);
    assertMatrixOrtho(Q, PREC);
    assertMatrixUT(R, PREC);
  }

  @Test public void lowrankQR() {
    Matrix M3 = Matrix.create(+1, +2, NR,
                              -3, +4, NR,
                              +5, -6, NR,
                              -7, -8);
    Matrix Q = Matrix.create(4, 4);
    Matrix R = Matrix.create(4, 2);
    M3.QR(Q, R, Matrix.create(3, 1));
    assertMatrixEquals(M3, Q.mul(R), PREC);
    assertMatrixOrtho(Q, PREC);
    assertMatrixUT(R, PREC);

    Matrix M4 = Matrix.create(+1, +2, NR,
                              -0, +0, NR,
                              +0, -6, NR,
                              -7, -8);
    Q = Matrix.rand(4, 4, RNG);
    R = Matrix.rand(4, 2, RNG);
    M4.QR(Q, R, Matrix.rand(3, 1, RNG));
    assertMatrixEquals(M4, Q.mul(R), PREC);
    assertMatrixOrtho(Q, PREC);
    assertMatrixUT(R, PREC);

    Matrix M5 = Matrix.create(  0,   0, NR, // size: 3x2
                               10,  20, NR,
                              200, 100);
    Q = Matrix.create(3, 3);
    R = Matrix.create(3, 2);
    M5.QR(Q, R, Matrix.create(10, 1));
    assertMatrixEquals(M5, Q.mul(R), PREC);
    assertMatrixOrtho(Q, PREC);
    assertMatrixUT(R, PREC);

    Matrix M6 = Matrix.create(0, 10, NR, // size: 3x2
                              0, 20, NR,
                              0, 30);
    Q = Matrix.create(3, 3);
    R = Matrix.create(3, 2);
    M6.QR(Q, R, Matrix.create(10, 1));
    assertMatrixEquals(M6, Q.mul(R), PREC);
    assertMatrixOrtho(Q, PREC);
    assertMatrixUT(R, PREC);

    Matrix M7 = Matrix.create(1, 2, 3, 4, NR, // size: 3x4
                              2, 4, 6, 8, NR,
                              8, 7, 6, 5);
    Q = Matrix.rand(3, 3, RNG);
    R = Matrix.rand(3, 4, RNG);
    M7.QR(Q, R, Matrix.rand(1, 2, RNG));
    assertMatrixEquals(M7, Q.mul(R), PREC);
    assertMatrixOrtho(Q, PREC);
    assertMatrixUT(R, PREC);
  }

  //---------------------------------------------------------------------------
  // LU decomposition tests

  @Test public void basicLU() {
    Matrix M, L, U, P;
    Matrix[] LUP;

    M = Matrix.zeros(2, 2);
    LUP = M.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M, P.T().mul(L).mul(U), PREC);
    assertMatrixEquals(Matrix.eye(2), L, PREC);
    assertMatrixEquals(Matrix.zeros(2,2), U, PREC);

    M = Matrix.eye(2);
    LUP = M.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M, P.T().mul(L).mul(U), PREC);
    assertMatrixEquals(Matrix.eye(2), L, PREC);
    assertMatrixEquals(Matrix.eye(2), U, PREC);

    M = M.mul(3.0);
    LUP = M.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M, P.T().mul(L).mul(U), PREC);
    assertMatrixEquals(Matrix.eye(2), L, PREC);
    assertMatrixEquals(Matrix.eye(2).mul(3), U, PREC);

    M = Matrix.create(2.0, 1.0, NR,
                      1.0, 2.0);
    LUP = M.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M, P.T().mul(L).mul(U), PREC);
    assertMatrixEquals(Matrix.create(1.0, 0.0, NR, 0.5, 1.0), L, PREC);
    assertMatrixEquals(Matrix.create(2.0, 1.0, NR, 0.0, 1.5), U, PREC);
    assertMatrixEquals(Matrix.eye(2), P, PREC);
  }

  @Test public void fullrankLU() {
    Matrix L, U, P;
    Matrix[] LUP;

    Matrix M1 = Matrix.create(1, 2, NR,
                              4, 1);
    LUP = M1.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M1, P.T().mul(L).mul(U), PREC);
    assertMatrixUnitLT(L, PREC);
    assertMatrixUT(U, PREC);

    M1 = M1.T();
    LUP = M1.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M1, P.T().mul(L).mul(U), PREC);
    assertMatrixUnitLT(L, PREC);
    assertMatrixUT(U, PREC);

    Matrix M3 = Matrix.create(17, 24,  1,  8, 15, NR,
                              23,  5,  7, 14, 16, NR,
                               4,  6, 13, 20, 22, NR,
                              10, 12, 19, 21,  3, NR,
                              11, 18, 25,  2,  9);
    L = Matrix.zeros(5,5);
    U = Matrix.zeros(5,5);
    Permutation p = Permutation.eye(5);
    M3.LU(L, U, p);
    assertMatrixEquals(M3, p.inv().mul(L).mul(U), PREC);
    assertMatrixUnitLT(L, PREC);
    assertMatrixUT(U, PREC);

    M3 = M3.T();
    M3.LU(L, U, p);
    assertMatrixEquals(M3, p.inv().mul(L).mul(U), PREC);
    assertMatrixUnitLT(L, PREC);
    assertMatrixUT(U, PREC);
  }

  @Test public void inplaceLU() {
    Matrix A = Matrix.create(17, 24,  1,  8, 15, NR,
                             23,  5,  7, 14, 16, NR,
                             4,  6, 13, 20, 22, NR,
                             10, 12, 19, 21,  3, NR,
                             11, 18, 25,  2,  9);
    Matrix B = Matrix.zeros(5,5);
    Permutation p = Permutation.eye(5);
    A.LU(B, B, p);
    Matrix L = B.getTriL();
    for (int i = 0; i < 5; ++i) { L.set(i, i, 1.0); }
    Matrix U = B.getTriU();
    assertMatrixEquals(A, p.inv().mul(L).mul(U), PREC);

    A = A.T();
    A.LU(B, B, p);
    L = B.getTriL();
    for (int i = 0; i < 5; ++i) { L.set(i, i, 1.0); }
    U = B.getTriU();
    assertMatrixEquals(A, p.inv().mul(L).mul(U), PREC);
  }

  @Test public void lowrankLU() {
    Matrix L, U, P;
    Matrix[] LUP;

    Matrix M2 = Matrix.create(1, 0, 2, NR, // rank: 2
                              0, 0, 0, NR,
                              4, 0, 1);
    LUP = M2.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M2, P.T().mul(L).mul(U), PREC);
    assertMatrixUnitLT(L, PREC);
    assertMatrixUT(U, PREC);

    M2 = M2.T();
    LUP = M2.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M2, P.T().mul(L).mul(U), PREC);
    assertMatrixUnitLT(L, PREC);
    assertMatrixUT(U, PREC);

    Matrix M4 = Matrix.create(1, 2, 3, 4, NR, // rank: 3
                              2, 4, 6, 8, NR,
                              8, 7, 6, 5);
    LUP = M4.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M4, P.T().mul(L.mul(U)), PREC);
    assertMatrixUnitLT(L, PREC);
    assertMatrixUT(U, PREC);

    M4 = M4.T();
    LUP = M4.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M4, P.T().mul(L.mul(U)), PREC);
    assertMatrixUnitLT(L, PREC);
    assertMatrixUT(U, PREC);

    Matrix M5 = Matrix.create(+1, +2, NR, // rank: 2
                              -3, +4, NR,
                              +5, -6, NR,
                              -7, -8);
    LUP = M5.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M5, P.T().mul(L.mul(U)), PREC);
    assertMatrixUnitLT(L, PREC);
    assertMatrixUT(U, PREC);

    M5 = M5.T();
    LUP = M5.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M5, P.T().mul(L.mul(U)), PREC);
    assertMatrixUnitLT(L, PREC);
    assertMatrixUT(U, PREC);

    Matrix M6 = Matrix.create( // rank: 5
                              96, 42, 20, 15,  0,  5, 50, 87, NR,
                              16, 19, 30,  8, 60, 28,  6, 27, NR,
                              28, 73, 36, 14, 78, 18, 89, 62, NR,
                              99, 80, 30, 53, 27,  7,  9, 28, NR,
                              49,  6, 31, 32, 94, 93, 14, 44);
    LUP = M6.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M6, P.T().mul(L.mul(U)), PREC);
    assertMatrixUnitLT(L, PREC);
    assertMatrixUT(U, PREC);

    M6 = M6.T();
    LUP = M6.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M6, P.T().mul(L.mul(U)), PREC);
    assertMatrixUnitLT(L, PREC);
    assertMatrixUT(U, PREC);

    Matrix M7 = Matrix.create( // rank: 3
                              61,   88,   2, 36,  32,  73,   99, 25,  80, NR,
                              83,   36,  36, 98,   2,  67,   51, 46,  33, NR,
                              89,  168, -14, 13, -17,  92,  153, 17, 121, NR,
                              -6, -132,  50, 85,  19, -25, -102, 29, -88, NR,
                              33,    8,  18, 59,  81,  54,   45, 33,  39);
    LUP = M7.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M7, P.T().mul(L.mul(U)), PREC);
    assertMatrixUnitLT(L, PREC);
    assertMatrixUT(U, PREC);

    M7 = M7.T();
    LUP = M7.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertMatrixEquals(M7, P.T().mul(L.mul(U)), PREC);
    assertMatrixUnitLT(L, PREC);
    assertMatrixUT(U, PREC);
  }

  //---------------------------------------------------------------------------
  // back substitution tests

  @Test public void backs() {
    Matrix A = Matrix.create(0.275186, 0.492146, 0.967304, 0.027302, NR,
                             0.889619, 0.760571, 0.496819, 0.499948, NR,
                             0.913949, 0.246029, 0.510040, 0.411547, NR,
                             0.769612, 0.280185, 0.554000, 0.086159);
    Matrix b = Matrix.create(0.71099, 0.78331, 0.53646, 0.28704).T();
    Matrix[] LU = A.LU();
    Matrix y = LU[0].backsL(LU[2].mul(b));
    Matrix x = LU[1].backsU(y);
    assertMatrixEquals(b, A.mul(x), PREC);

    Matrix L = Matrix.create(4,4);
    Matrix U = Matrix.create(4,4);
    Permutation P = Permutation.eye(4);
    A.LU(L, U, P);
    y = L.backsL(P.mul(b), true);
    x = U.backsU(y, false);
    assertMatrixEquals(b, A.mul(x), PREC);
  }

  //---------------------------------------------------------------------------
  // matrix inverse tests

  @Test public void inverseSmallMatrix() {
    Matrix i2 = Matrix.eye(2);
    Matrix m2 = Matrix.create(2.0, 1.0, NR,
                              1.0, 3.0);
    assertMatrixEquals(i2, m2.mul(m2.inv()), PREC);
    assertMatrixEquals(i2, m2.inv().mul(m2), PREC);

    Matrix i3 = Matrix.eye(3);
    Matrix m3 = Matrix.create(2.0, 1.0, 1.0, NR,
                              1.0, 3.0, 1.0, NR,
                              1.0, 1.0, 4.0);
    assertMatrixEquals(i3, m3.mul(m3.inv()), PREC);
    assertMatrixEquals(i3, m3.inv().mul(m3), PREC);
  }

  @Test public void inverseLargeMatrix() {
    Matrix i5 = Matrix.eye(5);
    Matrix m5 = Matrix.create(1.7, 2.4, 0.1, 0.8, 1.5, NR,
                              2.3, 0.5, 0.7, 1.4, 1.6, NR,
                              0.4, 0.6, 1.3, 2.0, 2.2, NR,
                              1.0, 1.2, 1.9, 2.1, 0.3, NR,
                              1.1, 1.8, 2.5, 0.2, 0.9);
    assertMatrixEquals(i5, m5.mul(m5.inv()), PREC);
    assertMatrixEquals(i5, m5.inv().mul(m5), PREC);

    Matrix m5b = m5.mul(m5.T());
    Matrix r1 = Matrix.rand(m5b.rows(), m5b.cols(), RNG);
    Matrix r2 = Matrix.rand(m5b.rows(), m5b.cols(), RNG);
    Permutation pr = Permutation.rand(m5b.rows(), RNG);
    assertMatrixEquals(i5, m5b.mul(m5b.inv(r1, r2, pr)), PREC);
    pr = Permutation.rand(m5b.rows(), RNG);
    assertMatrixEquals(i5, m5b.inv(r2, r1, pr).mul(m5b), PREC);
  }

  @Test public void inverseSmallMatrixPsd() {
    Matrix i2 = Matrix.eye(2);
    Matrix m2 = Matrix.create(2.0, 1.0, NR,
                              1.0, 3.0);
    assertMatrixEquals(i2, m2.mul(m2.invPsd()), PREC);
    assertMatrixEquals(i2, m2.invPsd().mul(m2), PREC);

    Matrix i3 = Matrix.eye(3);
    Matrix m3 = Matrix.create(2.0, 1.0, 1.0, NR,
                              1.0, 3.0, 1.0, NR,
                              1.0, 1.0, 4.0);
    assertMatrixEquals(i3, m3.mul(m3.invPsd()), PREC);
    assertMatrixEquals(i3, m3.invPsd().mul(m3), PREC);

    Matrix ip = m3.copy();
    ip.invPsd(ip); // in-place
    assertMatrixEquals(i3, m3.mul(ip), PREC);
    assertMatrixEquals(i3, ip.mul(m3), PREC);
  }

  @Test public void inverseLargeMatrixPsd() {
    Matrix i5 = Matrix.eye(5);
    Matrix m5 = Matrix.create(1.7, 2.4, 0.1, 0.8, 1.5, NR,
                              2.3, 0.5, 0.7, 1.4, 1.6, NR,
                              0.4, 0.6, 1.3, 2.0, 2.2, NR,
                              1.0, 1.2, 1.9, 2.1, 0.3, NR,
                              1.1, 1.8, 2.5, 0.2, 0.9);
    Matrix m5a = m5.T().mul(m5);
    assertMatrixEquals(i5, m5a.mul(m5a.invPsd()), PREC);
    assertMatrixEquals(i5, m5a.invPsd().mul(m5a), PREC);

    Matrix m5b = m5.mul(m5.T());
    Matrix r1 = Matrix.rand(m5.rows(), m5.cols(), RNG);
    assertMatrixEquals(i5, m5b.mul(m5b.invPsd(r1)), PREC);
    Matrix r2 = Matrix.rand(m5.rows(), m5.cols(), RNG);
    assertMatrixEquals(i5, m5b.invPsd(r2).mul(m5b), PREC);

    Matrix ip = m5a.copy();
    ip.invPsd(ip); // in-place
    assertMatrixEquals(i5, m5a.mul(ip), PREC);
    assertMatrixEquals(i5, ip.mul(m5a), PREC);
  }

  //---------------------------------------------------------------------------
  // determinant tests

  @Test public void determinant() {
    Matrix m2x2 = Matrix.create(1.0, 0.0, NR,
                                2.0, 3.0);
    assertEquals(3.0, m2x2.prodDiag(), PREC);
    assertEquals(3.0, m2x2.T().prodDiag(), PREC);

    assertEquals(3.0, m2x2.det(), PREC);
    assertEquals(3.0, m2x2.T().det(), PREC);
    assertEquals(9.0, m2x2.mul(m2x2.T()).det(), PREC);

    Matrix m3x3 = Matrix.create(1.0, 0.0, 0.0, NR,
                                2.0, 3.0, 0.0, NR,
                                4.0, 5.0, 6.0);
    assertEquals(18.0, m3x3.prodDiag(), PREC);
    assertEquals(18.0, m3x3.T().prodDiag(), PREC);

    assertEquals(18.0, m3x3.det(), PREC);
    assertEquals(18.0, m3x3.T().det(), PREC);
    assertEquals(324.0, m3x3.mul(m3x3.T()).det(), PREC);

    Matrix m5x5 = Matrix.create(1.7, 2.4, 0.1, 0.8, 1.5, NR,
                                2.3, 0.5, 0.7, 1.4, 1.6, NR,
                                0.4, 0.6, 1.3, 2.0, 2.2, NR,
                                1.0, 1.2, 1.9, 2.1, 0.3, NR,
                                1.1, 1.8, 2.5, 0.2, 0.9);
    assertEquals(50.7, m5x5.det(), PREC);
    assertEquals(50.7, m5x5.T().det(), PREC);

    m5x5 = Matrix.create(1.7, 2.4, 0.1, 0.8, 1.5, NR,
                         2.3, 0.5, 0.7, 1.4, 1.6, NR,
                         2.7, 3.6, 2.0, 2.9, 1.8, NR,
                         1.0, 1.2, 1.9, 2.1, 0.3, NR,
                         1.1, 1.8, 2.5, 0.2, 0.9);
    Matrix r = Matrix.rand(m5x5.rows(), m5x5.cols(), RNG);
    assertEquals(0.0, m5x5.det(r), PREC);
    assertEquals(0.0, m5x5.T().det(r.setToRand(RNG)), PREC);
  }

  //---------------------------------------------------------------------------
  // singular value decomposition tests

  @Test public void basicReducedSVD() {
    Matrix m2x2 = Matrix.zeros(2,2);
    Matrix[] USV = m2x2.reducedSVD();
    assertEquals(3, USV.length);
    Matrix U = USV[0], S = USV[1], V = USV[2];
    assertTrue(U.isEmpty());
    assertTrue(S.isEmpty());
    assertTrue(V.isEmpty());

    m2x2.setToEye();
    USV = m2x2.reducedSVD();
    assertEquals(3, USV.length);
    U = USV[0]; S = USV[1]; V = USV[2];
    assertMatrixEquals(Matrix.eye(2), U, PREC);
    assertMatrixEquals(Matrix.create(1.0, 1.0).T(), S, PREC);
    assertMatrixEquals(Matrix.eye(2), V, PREC);

    m2x2.set(1, 1, 22.0);
    USV = m2x2.reducedSVD();
    assertEquals(3, USV.length);
    U = USV[0]; S = USV[1]; V = USV[2];
    assertMatrixEquals(Matrix.eye(2), U, PREC);
    assertMatrixEquals(Matrix.create(1.0, 22.0).T(), S, PREC);
    assertMatrixEquals(Matrix.eye(2), V, PREC);

    m2x2 = Matrix.create(1.0, 2.0, NR,
                         2.0, 1.0);
    USV = m2x2.reducedSVD();
    assertEquals(3, USV.length);
    U = USV[0]; S = USV[1]; V = USV[2];
    assertMatrixEquals(m2x2, U.mul(Matrix.diag(S)).mul(V.T()), PREC);
    assertMatrixEquals(Matrix.create(1.0, 3.0).T(), S, PREC);
    assertMatrixOrtho(U, PREC);
    assertMatrixOrtho(V, PREC);

    m2x2 = Matrix.create(4.0, 4.0, NR,
                         -3.0, 3.0);
    USV = m2x2.reducedSVD();
    assertEquals(3, USV.length);
    U = USV[0]; S = USV[1]; V = USV[2];
    assertMatrixEquals(m2x2, U.mul(Matrix.diag(S)).mul(V.T()), PREC);
    assertMatrixEquals(Matrix.create(3.0, 4.0).T().mul(Math.sqrt(2.0)), S, PREC);
    assertMatrixOrtho(U, PREC);
    assertMatrixOrtho(V, PREC);
  }

  @Test public void fullrankReducedSVD() {
    Matrix m2x2 = Matrix.create(1.0, 0.0, NR,
                                2.0, 3.0);
    Matrix U = Matrix.create(2,2);
    Matrix S = Matrix.create(2,1);
    Matrix V = Matrix.create(2,2);
    assertEquals(2, m2x2.reducedSVD(U, S, V));
    assertMatrixEquals(m2x2, U.mul(Matrix.diag(S)).mul(V.T()), PREC);
    assertMatrixOrtho(U, PREC);
    assertMatrixOrtho(V, PREC);

    Matrix m4x2 = Matrix.create(+1, +2, NR,
                                -3, +4, NR,
                                +5, -6, NR,
                                -7, -8);
    U = Matrix.create(4,2);
    S = Matrix.create(2,1);
    V = Matrix.create(2,2);
    assertEquals(2, m4x2.reducedSVD(U, S, V));
    assertMatrixEquals(m4x2, U.mul(Matrix.diag(S)).mul(V.T()), PREC);
    assertMatrixOrthoCols(U, PREC);
    assertMatrixOrtho(V, PREC);

    Matrix m3x4 = Matrix.create(1, 2, 3, 4, NR,
                                2, 4, 7, 8, NR,
                                8, 7, 6, 5);
    U = Matrix.create(3,3);
    S = Matrix.create(3,1);
    V = Matrix.create(4,3);
    assertEquals(3, m3x4.reducedSVD(U, S, V));
    assertMatrixEquals(m3x4, U.mul(Matrix.diag(S)).mul(V.T()), PREC);
    assertMatrixOrtho(U, PREC);
    assertMatrixOrthoCols(V, PREC);
  }

  @Test public void lowrankReducedSVD() {
    Matrix m4x3 = Matrix.create(1, 2, 3, 4, NR,
                                2, 4, 6, 8, NR,
                                8, 7, 6, 5).T();
    Matrix U = Matrix.create(4,3);
    Matrix S = Matrix.create(3,1);
    Matrix V = Matrix.create(3,3);
    assertEquals(2, m4x3.reducedSVD(U, S, V));
    assertMatrixEquals(m4x3, U.mul(Matrix.diag(S)).mul(V.T()), PREC);
    U = U.getMat(0, 3, 0, 1);
    S = S.getMat(0, 1, 0, 0);
    V = V.getMat(0, 2, 0, 1);
    assertMatrixEquals(m4x3, U.mul(Matrix.diag(S)).mul(V.T()), PREC);
    assertMatrixOrthoCols(U, PREC);
    assertMatrixOrthoCols(V, PREC);

    Matrix m5x9 = Matrix.create(61,   88,   2, 36,  32,  73,   99, 25,  80, NR,
                                83,   36,  36, 98,   2,  67,   51, 46,  33, NR,
                                89,  168, -14, 13, -17,  92,  153, 17, 121, NR,
                                -6, -132,  50, 85,  19, -25, -102, 29, -88, NR,
                                33,    8,  18, 59,  81,  54,   45, 33,  39);

    Matrix[] USV = m5x9.reducedSVD();
    assertEquals(3, USV.length);
    U = USV[0]; S = USV[1]; V = USV[2];
    assertMatrixEquals(m5x9, U.mul(Matrix.diag(S)).mul(V.T()), PREC);
    assertMatrixOrthoCols(U, PREC);
    assertMatrixOrthoCols(V, PREC);

    Matrix m5x9save = m5x9.copy();
    U = Matrix.create(5,5);
    S = Matrix.create(5,1);
    V = m5x9.T();
    assertEquals(3, m5x9.reducedSVD(U, S, V));
    assertMatrixEquals(m5x9save, U.mul(Matrix.diag(S)).mul(V.T()), PREC);
    U = U.getMat(0, 4, 0, 2);
    S = S.getMat(0, 2, 0, 0);
    V = V.getMat(0, 8, 0, 2);
    assertMatrixEquals(m5x9save, U.mul(Matrix.diag(S)).mul(V.T()), PREC);
    assertMatrixOrthoCols(U, PREC);
    assertMatrixOrthoCols(V, PREC);
  }
}
