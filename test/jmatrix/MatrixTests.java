package jmatrix;

import java.util.Random;
import static jmatrix.Matrix.NR;
import static jmatrix.BasicUnaryOperation.*;
import static jmatrix.BasicBinaryOperation.*;

/**
 * Tests for the Matrix class.
 */
public class MatrixTests extends AssertionBaseTest {

  public static final double PREC = 1e-8;
  public static final Random RNG = new Random();

  //---------------------------------------------------------------------------

  public MatrixTests(String name) { super(name); }

  //---------------------------------------------------------------------------

  public void testIsEmpty() {
    Matrix e = Matrix.create(0,0);
    assertTrue(e.isEmpty());
    assertTrue(e.rows() == 0);
    assertTrue(e.cols() == 0);

    e = Matrix.create(0,1);
    assertTrue(e.isEmpty());
    assertTrue(e.rows() == 0);
    assertTrue(e.cols() == 1);

    e = Matrix.create(1,0);
    assertTrue(e.isEmpty());
    assertTrue(e.rows() == 1);
    assertTrue(e.cols() == 0);

    e = Matrix.create(1,1);
    assertFalse(e.isEmpty());
    assertTrue(e.rows() == 1);
    assertTrue(e.cols() == 1);
  }

  public void testNorms() {
    Matrix m = Matrix.create(1.0, 1.5, -0.5, 1.5, NR,
                             2.0, 4.0, -4.0, 0.0, NR,
                             0.5, 3.0, -3.0, 2.0);
    assertTrue(PREC > Math.abs(8.5 - m.norm1()));
    assertTrue(PREC > Math.abs(10.0 - m.normI()));
    assertTrue(PREC > Math.abs(8.0 - m.normF()));

    m = Matrix.create(new double[][]{ // test this form of create
        new double[]{1.0, 1.5, -0.5, 1.5},
        new double[]{2.0, 4.0, -4.0, 0.0},
        new double[]{0.5, 3.0, -3.0, 2.0}
      });
    assertTrue(PREC > Math.abs(8.5 - m.norm1()));
    assertTrue(PREC > Math.abs(10.0 - m.normI()));
    assertTrue(PREC > Math.abs(8.0 - m.normF()));
  }

  public void testNaNandInfChecks() {
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

  public void testBasicLinearOps() {
    final int n = 8, m = 17;
    Matrix m1 = Matrix.zeros(n, m);
    Matrix m2 = Matrix.ones(n, m);
    Matrix m3 = Matrix.scalars(n, m, 42.42);
    Matrix m4 = Matrix.eye(n);

    assertTrue(PREC > m2.sub(m2.add(m1)).normF());
    assertTrue(PREC > m3.sub(m2.mul(42.42)).normF());
    assertTrue(PREC > m2.sub(m3.div(42.42)).norm1());
    assertTrue(PREC > m4.sub(m4).normI());

    assertTrue(PREC > m1.add(1.0).sub(m2).normF());
    assertTrue(PREC > m2.sub(1.0).normF());

    m2.add(m2, m2);
    assertTrue(PREC > m2.sub(m3.div(42.42).mul(2)).norm1());
    m2.setToOnes();

    m2.sub(m2.div(2), m2);
    assertTrue(PREC > m2.sub(m3.div(42.42).mul(0.5)).normI());
    m2.setToOnes();

    m2.mul(2, m2);
    assertTrue(PREC > m2.sub(m3.div(42.42).mul(2)).normF());
    m2.setToOnes();

    m2.div(8, m2);
    assertTrue(PREC > m2.sub(m3.div(8*42.42)).normI());
    m2.setToOnes();
  }

  public void testEwu() {
    // ABS, NEG, SIGN
    Matrix m = Matrix.randN(4, 6, RNG);
    assertTrue(PREC > m.sub(m.ewu(SIGN).emul(m.ewu(ABS))).norm1());
    assertTrue(PREC > m.ewu(NEG).sub(m.ewu(NEG).ewu(SIGN).emul(m.ewu(ABS))).norm1());

    m.ewu(NEG);
    assertTrue(PREC > m.sub(m.ewu(SIGN).emul(m.ewu(ABS))).norm1());

    Matrix s = Matrix.zeros(4, 6);
    Matrix a = Matrix.zeros(4, 6);
    m.copy(s).ewu(SIGN, s);
    m.copy(a).ewu(ABS, m);
    assertTrue(PREC > m.sub(s.emul(a)).norm1());

    s.copy(a).ewu(SIGN, s);
    assertTrue(PREC > s.sub(a).norm1());

    s.set(2, 2, 0.0);
    a.set(2, 2, 0.0);
    assertTrue(PREC > s.sub(a).norm1());
  }

  public void testEwb12() {
    // RECIPROC
    int nr = 2, nc = 5;
    Matrix v = Matrix.scalars(nr, nc, 0.1);
    Matrix o = Matrix.ones(nr, nc);
    Matrix r = Matrix.zeros(nr, nc);

    assertTrue(PREC > o.sub(v.emul(v.ewb1(DIV, 1.0))).normF());

    v.ewb1(DIV, 1.0, r);
    assertTrue(PREC > o.sub(v.emul(r)).normF());

    r.ewb1(DIV, 1.0, r);
    assertTrue(PREC > v.sub(r).normF());

    // POW
    Matrix m1 = Matrix.create(1, 2, 3, NR,
                              4, 5, 6);
    Matrix m2 = Matrix.create(1, 4, 9, NR,
                              16, 25, 36);
    assertTrue(PREC > m2.sub(m1.ewb2(POW, 2.0)).norm1());

    Matrix mr = Matrix.randN(4, 5, RNG);
    assertTrue(PREC > mr.emul(mr).sub(mr.ewb2(POW, 2.0)).norm1());
    assertTrue(PREC > mr.emul(mr).emul(mr).sub(mr.ewb2(POW, 3.0)).norm1());

    // MOD
    double m = 2.3;
    Matrix mm = Matrix.rand(5, 4, RNG);
    mm.mul(8.0, mm);
    Matrix oo = Matrix.scalars(5, 4, m);

    assertTrue(PREC > oo.ewb2(MOD, m).norm1());
    assertTrue(PREC > mm.add(oo.mul(2.0)).ewb2(MOD, m).sub(mm.ewb2(MOD, m)).norm1());

    mm.mul(-1.0, mm);
    assertTrue(PREC > mm.add(oo.mul(-2.0)).ewb2(MOD, m).sub(mm.ewb2(MOD, m)).norm1());

    oo.ewb2(MOD, m, oo);
    assertTrue(PREC > oo.norm1());
  }

  public void testEwb() {
    Matrix m = Matrix.create(1, 2, 3, NR,
                             4, 5, 6);

    Matrix r = Matrix.create(1, 2, 3);
    Matrix t1 = Matrix.create(1,  4,  9, NR,
                              4, 10, 18);
    assertTrue(PREC > m.ewb(MUL, r).sub(t1).normI());
    m.ewb(MUL, r, m);
    assertTrue(PREC > m.sub(t1).normI());

    m = Matrix.create(1, 2, 3, NR,
                      4, 5, 6);

    Matrix c = Matrix.create(1.0, 2.0).T();
    Matrix t2 = Matrix.create(1,  2,  3, NR,
                              8, 10, 12);
    assertTrue(PREC > m.ewb(MUL, c).sub(t2).normI());
    m.ewb(MUL, c, m);
    assertTrue(PREC > m.sub(t2).normI());
  }

  public void testCopy() {
    Matrix m = Matrix.rand(5, 7, RNG);
    Matrix c = m.copy();

    m.add(Matrix.randN(5, 7, RNG), m);
    assertTrue(PREC > m.sub(m).norm1());
    assertTrue(PREC < c.sub(m).norm1());

    c.copy(m);
    assertTrue(PREC > m.sub(c).norm1());
  }

  public void testBlkdiag() {
    Matrix m1 = Matrix.create(1.1, 2.2);
    Matrix m2 = Matrix.create(3.3, 4.4, 5.5, NR,
                              6.6, 7.7, 8.8);
    Matrix m3 = Matrix.create(9.9);

    Matrix m = Matrix.blkdiag(m1, m2, m3);
    Matrix t = Matrix.create(1.1, 2.2, 0.0, 0.0, 0.0, 0.0, NR,
                             0.0, 0.0, 3.3, 4.4, 5.5, 0.0, NR,
                             0.0, 0.0, 6.6, 7.7, 8.8, 0.0, NR,
                             0.0, 0.0, 0.0, 0.0, 0.0, 9.9);
    assertTrue(PREC > t.sub(m).normI());
  }

  public void testCat() {
    Matrix m1 = Matrix.create(0.0, 1.1, NR,
                              0.0, 2.2);
    Matrix m2 = Matrix.create(3.3, 4.4, 5.5, NR,
                              6.6, 7.7, 8.8);
    Matrix m3 = Matrix.create(9.9, NR,
                              9.9);

    Matrix m = Matrix.horzcat(m1, m2, m3);
    Matrix t = Matrix.create(0.0, 1.1, 3.3, 4.4, 5.5, 9.9, NR,
                             0.0, 2.2, 6.6, 7.7, 8.8, 9.9);
    assertTrue(PREC > t.sub(m).normI());

    m1 = m1.T();
    m2 = m2.T();
    m3 = m3.T();
    m = Matrix.vertcat(m1, m2, m3);
    t = t.T();
    assertTrue(PREC > t.sub(m).normI());
  }

  public void testSetTo() {
    Matrix m = Matrix.ones(4, 7);
    m.setToZeros();
    assertTrue(PREC > Matrix.zeros(4, 7).sub(m).norm1());
    m.setToOnes();
    assertTrue(PREC > Matrix.ones(4, 7).sub(m).norm1());

    Matrix I = Matrix.eye(4);
    I.setToScalars(2);
    assertTrue(PREC > Matrix.scalars(4, 4, 2).sub(I).norm1());
    I.setToEye();
    assertTrue(PREC > Matrix.eye(4).sub(I).norm1());
  }

  public void testGetSetMat() {
    Matrix A = Matrix.rand(7, 7, RNG);
    Matrix B = Matrix.rand(7, 7, RNG);

    Matrix AB = Matrix.scalars(10, 10, Double.NaN);
    AB.setMat(2, 8, 2, 8, A.mul(B));

    assertTrue(PREC > A.getMat(0, 0, 0, 6).mul(B.getMat(0, 6, 2, 3))
                       .sub(AB.getMat(2, 2, 4, 5)).normF());
    assertTrue(PREC > A.getMat(4, 6, 0, 6).mul(B.getMat(0, 6, 2, 3))
                       .sub(AB.getMat(6, 8, 4, 5)).normF());
    assertTrue(PREC > A.getMat(1, 2, 0, 6).mul(B.getMat(0, 6, 2, 3))
                       .sub(AB.getMat(3, 4, 4, 5)).normF());
  }

  public void testGetDiag() {
    Matrix A = Matrix.create(1, 2, 3, NR,
                             4, 5, 6, NR,
                             7, 8, 9);
    Matrix a = Matrix.create(1, 5, 9).T();
    assertEquals(0.0, A.getDiag().sub(a).normI());
    assertEquals(0.0, A.T().getDiag().sub(a).normI());

    Matrix B = Matrix.create(1, 2, 3, 4, 5, NR,
                             6, 7, 8, 9, 0);
    Matrix b = Matrix.create(1.0, 7.0).T();
    assertEquals(0.0, B.getDiag().sub(b).normI());
    assertEquals(0.0, B.T().getDiag().sub(b).normI());
  }

  public void testGetTrilLU() {
    Matrix A = Matrix.create(1, 2, 3, NR,
                             4, 5, 6, NR,
                             7, 8, 9);

    Matrix L1 = Matrix.create(1, 0, 0, NR,
                              4, 5, 0, NR,
                              7, 8, 9);
    assertEquals(0.0, A.getTriL().sub(L1).normF());

    Matrix L2 = Matrix.create(0, 0, 0, NR,
                              4, 0, 0, NR,
                              7, 8, 0);
    assertEquals(0.0, A.getTriL(-1).sub(L2).normF());

    Matrix L3 = Matrix.create(1, 2, 0, NR,
                              4, 5, 6, NR,
                              7, 8, 9);
    assertEquals(0.0, A.getTriL(1).sub(L3).normF());

    assertEquals(0.0, A.getTriL(-3).normF());
    assertEquals(0.0, A.getTriL(-5).normF());
    assertEquals(0.0, A.getTriL(2).sub(A).normF());
    assertEquals(0.0, A.getTriL(5).sub(A).normF());

    Matrix U1 = Matrix.create(1, 2, 3, NR,
                              0, 5, 6, NR,
                              0, 0, 9);
    assertEquals(0.0, A.getTriU().sub(U1).normF());

    Matrix U2 = Matrix.create(1, 2, 3, NR,
                              4, 5, 6, NR,
                              0, 8, 9);
    assertEquals(0.0, A.getTriU(-1).sub(U2).normF());

    Matrix U3 = Matrix.create(0, 2, 3, NR,
                              0, 0, 6, NR,
                              0, 0, 0);
    assertEquals(0.0, A.getTriU(1).sub(U3).normF());

    assertEquals(0.0, A.getTriU(-2).sub(A).normF());
    assertEquals(0.0, A.getTriU(-5).sub(A).normF());
    assertEquals(0.0, A.getTriU(3).normF());
    assertEquals(0.0, A.getTriU(5).normF());
  }

  public void testMatrixVectorProduct() {
    final int n = 5, m = 8;
    Matrix v1 = Matrix.ones(m, 1);
    Matrix v2 = Matrix.zeros(m, 1); v2.set(m/2, 0, 1.0); // unit vector
    Matrix one = Matrix.ones(n, 1);
    Matrix m1 = Matrix.ones(n, m);
    Matrix m2 = Matrix.eye(m);

    assertTrue(PREC > v1.sub(m2.mul(v1)).norm1());
    assertTrue(PREC > one.sub(m1.mul(v2)).norm1());

    for (int i = 0; i < m; ++i) v2.set(i, 0, i+1);
    double sum = m * (m+1) / 2.0;
    assertTrue(PREC > one.mul(sum).sub(m1.mul(v2)).norm1());
    assertTrue(PREC > v2.sub(m2.mul(v2)).normI());
    assertTrue(PREC < v2.sub(Matrix.ones(m, m).mul(v2)).normI());

    assertTrue(PREC > v2.emul(v2).sub(Matrix.diag(v2).mul(v2)).norm1());
    assertTrue(PREC > v2.emul(v2).sub(Matrix.diag(v2.T()).mul(v2)).norm1());
    assertTrue(PREC > m1.mul(v1.add(v2)).sub(m1.mul(v1).add(m1.mul(v2))).norm1());
  }

  public void testMatrixProduct() {
    final int n = 5, m = 8;
    Matrix m1 = Matrix.ones(n, m).mul(2);
    Matrix m2 = Matrix.eye(n);
    Matrix m3 = Matrix.eye(m);

    assertTrue(PREC > m1.sub(m1.mul(m3)).norm1());
    assertTrue(PREC > m1.sub(m2.mul(m1)).norm1());

    m2.mul(3, m2);
    m3.mul(5, m3);

    assertTrue(PREC > m1.mul(m3).sub(m3.T().mul(m1.T()).T()).normF());
    assertTrue(PREC > m2.mul(m1).sub(m1.T().mul(m2.T()).T()).normF());

    Matrix m4 = Matrix.ones(n, m).mul(7);

    assertTrue(PREC > m1.add(m4).mul(m3)
                        .sub(m1.mul(m3).add(m4.mul(m3))).normI());
    assertTrue(PREC > m2.mul(m1.add(m4))
                        .sub(m2.mul(m1).add(m2.mul(m4))).normI());

    Matrix m5 = Matrix.ones(n, m).mul(-8);

    m1.mul(m3, m5);
    assertTrue(PREC > m5.sub(m1.mul(m3)).norm1());

    m2.mul(m1, m5);
    assertTrue(PREC > m5.sub(m2.mul(m1)).norm1());
  }

  public void testEntrywiseMultiplication() {
    final int n = 5, m = 7;
    Matrix m1 = Matrix.scalars(n, m, 3);
    Matrix m2 = Matrix.scalars(n, m, 9);
    Matrix m3 = Matrix.zeros(n, m);

    assertTrue(PREC > m1.emul(m1).sub(m2).norm1());
    assertTrue(PREC > m3.sub(m1.emul(m3)).norm1());
    assertTrue(PREC > m3.sub(m3.emul(m2)).norm1());

    m1.emul(m2, m1);
    assertTrue(PREC > m1.sub(m2.mul(3)).normI());
    m1.setToScalars(3.0);
  }

  public void testDotProduct() {
    Matrix v1 = Matrix.create(1, 2, 0, 3);
    Matrix v2 = Matrix.create(5, 2, 9, -1);

    assertEquals(6.0, v1.dot(v2));
    assertEquals(6.0, v1.T().dot(v2));
    assertEquals(6.0, v1.dot(v2.T()));
    assertEquals(6.0, v1.T().dot(v2.T()));

    assertEquals(6.0, v2.dot(v1));
    assertEquals(6.0, v2.T().dot(v1));
    assertEquals(6.0, v2.dot(v1.T()));
    assertEquals(6.0, v2.T().dot(v1.T()));
  }

  public void testTrace() {
    Matrix m = Matrix.create(1.0, 1.5, -0.5, NR,
                             2.0, 4.0, -4.0, NR,
                             0.5, 3.0, -3.0, NR);
    assertEquals(2.0, m.trace());
    assertEquals(-12.0, m.prodDiag());
    assertEquals(4.0, Matrix.eye(4).trace());
    assertEquals(8.0, Matrix.scalars(4, 4, 2).trace());

    Matrix m3x2 = Matrix.rand(3, 2, RNG);
    Matrix m2x3 = Matrix.rand(2, 3, RNG);
    assertEquals(m3x2.mul(m2x3).trace(), m3x2.traceMul(m2x3), PREC);
  }

  public void testRowColNorms() {
    Matrix m = Matrix.create(17, 24, -1,  8, 15, 1, NR,
                             23, -5, -7, 14, 16, 1, NR,
                             -4,  6, 31, 20, 22, 1, NR,
                             10, 12, 19, 21,  3, 1, NR,
                             11, 18, 25, -2,  9, 1);

    Matrix n1 = Matrix.create(66, 66, 84, 66, 66).T();
    assertTrue(PREC > m.rowNorms1().sub(n1).normI());
    n1 = Matrix.create(65, 65, 83, 65, 65, 5);
    assertTrue(PREC > m.colNorms1().sub(n1).normI());

    Matrix n2 = Matrix.create(34.0,
                              32.4961536185438,
                              43.5660418215839,
                              32.4961536185438,
                              34.0).T();
    assertTrue(PREC > m.rowNorms2().sub(n2).normI());
    n2 = Matrix.create(32.48076353782343,
                       33.24154027718932,
                       44.68780594300866,
                       33.24154027718932,
                       32.48076353782343,
                       2.23606797749979);
    assertTrue(PREC > m.colNorms2().sub(n2).normI());

    Matrix nI = Matrix.create(24, 23, 31, 21, 25).T();
    assertTrue(PREC > m.rowNormsI().sub(nI).normI());
    nI = Matrix.create(23, 24, 31, 21, 22, 1);
    assertTrue(PREC > m.colNormsI().sub(nI).normI());
  }

  public void testRowColSums() {
    Matrix m = Matrix.create(17, 24,  1,  8, 15, 1, NR,
                             23,  5,  7, 14, 16, 1, NR,
                              4,  6, 31, 20, 22, 1, NR,
                             10, 12, 19, 21,  3, 1, NR,
                             11, 18, 25,  2,  9, 1);

    Matrix vRow = Matrix.create(65, 65, 83, 65, 65, 5);
    assertTrue(PREC > m.sumRows().sub(vRow).norm1());

    Matrix vCol = Matrix.create(66, 66, 84, 66, 66).T();
    assertTrue(PREC > m.sumCols().sub(vCol).norm1());
  }

  public void testCholeskyDecomposition() {
    Matrix m = Matrix.create(2.0, 1.0, 1.0, NR,
                             1.0, 2.0, 1.0, NR,
                             1.0, 1.0, 2.0);
    Matrix L = m.choleskyL();
    assertTrue(PREC > m.sub(L.mul(L.T())).norm1());

    m.set(0, 0, 3.33333); m.set(1, 1, 4.44444); m.set(2, 2, 5.55555);
    m.choleskyL(L);
    assertTrue(PREC > m.sub(L.mul(L.T())).norm1());

    L = m.copy();
    L.choleskyL(L); // in-place
    assertTrue(PREC > m.sub(L.mul(L.T())).norm1());
  }

  public void testCholeskyDecompositionLD() {
    Matrix m = Matrix.create(2.0, 1.0, 1.0, NR,
                             1.0, 3.0, 1.0, NR,
                             1.0, 1.0, 2.0);
    Matrix LD[] = m.choleskyLD();
    Matrix L = LD[0];
    Matrix D = LD[1];
    assertTrue(PREC > m.sub(L.mul(D).mul(L.T())).norm1());

    L = m.copy();
    D = Matrix.create(3, 1);
    L.choleskyLD(L, D); // "in-place"
    assertTrue(PREC > m.sub(L.mul(Matrix.diag(D)).mul(L.T())).norm1());
  }

  public void testQR() {
    Matrix M1 = Matrix.create(1, 2, NR, // size: 2x2
                              4, 1);
    Matrix[] QR = M1.QR();
    Matrix Q = QR[0];
    Matrix R = QR[1];
    assertTrue(PREC > M1.sub(Q.mul(R)).norm1());
    assertTrue(PREC > Matrix.eye(2).sub(Q.T().mul(Q)).norm1());

    Matrix M2 = Matrix.create(17, 24,  1,  8, 15, NR, // size: 5x5
                              23,  5,  7, 14, 16, NR,
                               4,  6, 13, 20, 22, NR,
                              10, 12, 19, 21,  3, NR,
                              11, 18, 25,  2,  9);
    QR = M2.QR();
    Q = QR[0];
    R = QR[1];
    assertTrue(PREC > M2.sub(Q.mul(R)).norm1());
    assertTrue(PREC > Matrix.eye(5).sub(Q.T().mul(Q)).norm1());

    Matrix M3 = Matrix.create(+1, +2, NR, // size: 4x2
                              -3, +4, NR,
                              +5, -6, NR,
                              -7, -8);
    Q = Matrix.create(4, 4);
    R = Matrix.create(4, 2);
    M3.QR(Q, R, Matrix.create(3, 1));
    assertTrue(PREC > M3.sub(Q.mul(R)).norm1());
    assertTrue(PREC > Matrix.eye(4).sub(Q.T().mul(Q)).norm1());

    Matrix M4 = Matrix.create(+1, +2, NR, // size: 4x2
                              -0, +0, NR,
                              +0, -6, NR,
                              -7, -8);
    Q = Matrix.rand(4, 4, RNG);
    R = Matrix.rand(4, 2, RNG);
    M4.QR(Q, R, Matrix.rand(3, 1, RNG));
    assertTrue(PREC > M4.sub(Q.mul(R)).norm1());
    assertTrue(PREC > Matrix.eye(4).sub(Q.T().mul(Q)).norm1());

    Matrix M5 = Matrix.create(  0,   0, NR, // size: 3x2
                               10,  20, NR,
                              200, 100);
    Q = Matrix.create(3, 3);
    R = Matrix.create(3, 2);
    M5.QR(Q, R, Matrix.create(10, 1));
    assertTrue(PREC > M5.sub(Q.mul(R)).norm1());
    assertTrue(PREC > Matrix.eye(3).sub(Q.T().mul(Q)).norm1());

    Matrix M6 = Matrix.create(0, 10, NR, // size: 3x2
                              0, 20, NR,
                              0, 30);
    Q = Matrix.create(3, 3);
    R = Matrix.create(3, 2);
    M6.QR(Q, R, Matrix.create(10, 1));
    assertTrue(PREC > M6.sub(Q.mul(R)).norm1());
    assertTrue(PREC > Matrix.eye(3).sub(Q.T().mul(Q)).norm1());

    Matrix M7 = Matrix.create(1, 2, 3, 4, NR, // size: 3x4
                              2, 4, 6, 8, NR,
                              8, 7, 6, 5);
    Q = Matrix.rand(3, 3, RNG);
    R = Matrix.rand(3, 4, RNG);
    M7.QR(Q, R, Matrix.rand(1, 2, RNG));
    assertTrue(PREC > M7.sub(Q.mul(R)).norm1());
    assertTrue(PREC > Matrix.eye(3).sub(Q.T().mul(Q)).norm1());
  }

  public void testLU() {
    Matrix L, U, P;

    Matrix M1 = Matrix.create(1, 2, NR, // size: 2x2
                              4, 1);
    Matrix[] LUP = M1.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertTrue(PREC > M1.sub(P.T().mul(L).mul(U)).norm1());

    Matrix M2 = Matrix.create(1, 0, 2, NR, // size: 3x3, singular
                              0, 0, 0, NR,
                              4, 0, 1);
    LUP = M2.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertTrue(PREC > M2.sub(P.T().mul(L).mul(U)).norm1());

    Matrix M3 = Matrix.create(17, 24,  1,  8, 15, NR, // size: 5x5
                              23,  5,  7, 14, 16, NR,
                               4,  6, 13, 20, 22, NR,
                              10, 12, 19, 21,  3, NR,
                              11, 18, 25,  2,  9);
    L = Matrix.rand(M3.rows(), Math.min(M3.rows(), M3.cols()), RNG);
    U = Matrix.rand(Math.min(M3.rows(), M3.cols()), M3.cols(), RNG);
    Permutation p = Permutation.rand(M3.rows(), RNG);
    M3.LU(L, U, p);
    assertTrue(PREC > M3.sub(p.inv().mul(L).mul(U)).norm1());

    Matrix M4 = Matrix.create(1, 2, 3, 4, NR, // size: 3x4
                              2, 4, 6, 8, NR,
                              8, 7, 6, 5);
    LUP = M4.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertTrue(PREC > M4.sub(P.T().mul(L.mul(U))).norm1());

    Matrix M5 = Matrix.create(+1, +2, NR, // size: 4x2
                              -3, +4, NR,
                              +5, -6, NR,
                              -7, -8);
    LUP = M5.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertTrue(PREC > M5.sub(P.T().mul(L.mul(U))).norm1());

    Matrix M6 = Matrix.create(96, 42, 20, 15,  0,  5, 50, 87, NR, // size: 5x8, rank: 5
                              16, 19, 30,  8, 60, 28,  6, 27, NR,
                              28, 73, 36, 14, 78, 18, 89, 62, NR,
                              99, 80, 30, 53, 27,  7,  9, 28, NR,
                              49,  6, 31, 32, 94, 93, 14, 44);
    LUP = M6.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertTrue(PREC > M6.sub(P.T().mul(L).mul(U)).norm1());

    M6 = M6.T(); // size: 8x5
    LUP = M6.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertTrue(PREC > M6.sub(P.T().mul(L).mul(U)).norm1());

    Matrix M7 = Matrix.create( // size: 5x9, rank: 3
                              61,   88,   2, 36,  32,  73,   99, 25,  80, NR,
                              83,   36,  36, 98,   2,  67,   51, 46,  33, NR,
                              89,  168, -14, 13, -17,  92,  153, 17, 121, NR,
                              -6, -132,  50, 85,  19, -25, -102, 29, -88, NR,
                              33,    8,  18, 59,  81,  54,   45, 33,  39);
    LUP = M7.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertTrue(PREC > M7.sub(P.T().mul(L).mul(U)).norm1());

    M7 = M7.T(); // size: 9x5
    LUP = M7.LU(); L = LUP[0]; U = LUP[1]; P = LUP[2];
    assertTrue(PREC > M7.sub(P.T().mul(L).mul(U)).norm1());
  }

  public void testBacks() {
    Matrix A = Matrix.create(0.275186, 0.492146, 0.967304, 0.027302, NR,
                             0.889619, 0.760571, 0.496819, 0.499948, NR,
                             0.913949, 0.246029, 0.510040, 0.411547, NR,
                             0.769612, 0.280185, 0.554000, 0.086159);
    Matrix b = Matrix.create(0.71099, 0.78331, 0.53646, 0.28704).T();
    Matrix[] LU = A.LU();
    Matrix y = LU[0].backsL(LU[2].mul(b));
    assertTrue(!y.hasNaN());
    Matrix x = LU[1].backsU(y);
    assertTrue(!x.hasNaN());
    assertTrue(PREC > A.mul(x).sub(b).normI());

    Matrix L = Matrix.create(4,4);
    Matrix U = Matrix.create(4,4);
    Permutation P = Permutation.eye(4);
    A.LU(L, U, P);
    y = L.backsL(P.mul(b), true);
    assertTrue(!y.hasNaN());
    x = U.backsU(y, false);
    assertTrue(!x.hasNaN());
    assertTrue(PREC > A.mul(x).sub(b).normI());
  }

  public void testInverseSmallMatrix() {
    Matrix i2 = Matrix.eye(2);
    Matrix m2 = Matrix.create(2.0, 1.0, NR,
                              1.0, 3.0);
    assertTrue(PREC > i2.sub(m2.mul(m2.inv())).norm1());
    assertTrue(PREC > i2.sub(m2.inv().mul(m2)).norm1());

    Matrix i3 = Matrix.eye(3);
    Matrix m3 = Matrix.create(2.0, 1.0, 1.0, NR,
                              1.0, 3.0, 1.0, NR,
                              1.0, 1.0, 4.0);
    assertTrue(PREC > i3.sub(m3.mul(m3.inv())).norm1());
    assertTrue(PREC > i3.sub(m3.inv().mul(m3)).norm1());
  }

  public void testInverseLargeMatrix() {
    Matrix i5 = Matrix.eye(5);
    Matrix m5 = Matrix.create(1.7, 2.4, 0.1, 0.8, 1.5, NR,
                              2.3, 0.5, 0.7, 1.4, 1.6, NR,
                              0.4, 0.6, 1.3, 2.0, 2.2, NR,
                              1.0, 1.2, 1.9, 2.1, 0.3, NR,
                              1.1, 1.8, 2.5, 0.2, 0.9);
    assertTrue(PREC > i5.sub(m5.mul(m5.inv())).norm1());
    assertTrue(PREC > i5.sub(m5.inv().mul(m5)).norm1());

    Matrix m5b = m5.mul(m5.T());
    Matrix r1 = Matrix.rand(m5b.rows(), m5b.cols(), RNG);
    Matrix r2 = Matrix.rand(m5b.rows(), m5b.cols(), RNG);
    Permutation pr = Permutation.rand(m5b.rows(), RNG);
    assertTrue(PREC > i5.sub(m5b.mul(m5b.inv(r1, r2, pr))).norm1());
    pr = Permutation.rand(m5b.rows(), RNG);
    assertTrue(PREC > i5.sub(m5b.inv(r2, r1, pr).mul(m5b)).norm1());
  }

  public void testInverseSmallMatrixPsd() {
    Matrix i2 = Matrix.eye(2);
    Matrix m2 = Matrix.create(2.0, 1.0, NR,
                              1.0, 3.0);
    assertTrue(PREC > i2.sub(m2.mul(m2.invPsd())).norm1());
    assertTrue(PREC > i2.sub(m2.invPsd().mul(m2)).norm1());

    Matrix i3 = Matrix.eye(3);
    Matrix m3 = Matrix.create(2.0, 1.0, 1.0, NR,
                              1.0, 3.0, 1.0, NR,
                              1.0, 1.0, 4.0);
    assertTrue(PREC > i3.sub(m3.mul(m3.invPsd())).norm1());
    assertTrue(PREC > i3.sub(m3.invPsd().mul(m3)).norm1());

    Matrix ip = m3.copy();
    ip.invPsd(ip); // in-place
    assertTrue(PREC > i3.sub(m3.mul(ip)).norm1());
    assertTrue(PREC > i3.sub(ip.mul(m3)).norm1());
  }

  public void testInverseLargeMatrixPsd() {
    Matrix i5 = Matrix.eye(5);
    Matrix m5 = Matrix.create(1.7, 2.4, 0.1, 0.8, 1.5, NR,
                              2.3, 0.5, 0.7, 1.4, 1.6, NR,
                              0.4, 0.6, 1.3, 2.0, 2.2, NR,
                              1.0, 1.2, 1.9, 2.1, 0.3, NR,
                              1.1, 1.8, 2.5, 0.2, 0.9);
    Matrix m5a = m5.T().mul(m5);
    assertTrue(PREC > i5.sub(m5a.mul(m5a.invPsd())).norm1());
    assertTrue(PREC > i5.sub(m5a.invPsd().mul(m5a)).norm1());

    Matrix m5b = m5.mul(m5.T());
    Matrix r1 = Matrix.rand(m5.rows(), m5.cols(), RNG);
    assertTrue(PREC > i5.sub(m5b.mul(m5b.invPsd(r1))).norm1());
    Matrix r2 = Matrix.rand(m5.rows(), m5.cols(), RNG);
    assertTrue(PREC > i5.sub(m5b.invPsd(r2).mul(m5b)).norm1());

    Matrix ip = m5a.copy();
    ip.invPsd(ip); // in-place
    assertTrue(PREC > i5.sub(m5a.mul(ip)).norm1());
    assertTrue(PREC > i5.sub(ip.mul(m5a)).norm1());
  }

  public void testDeterminant() {
    Matrix m2x2 = Matrix.create(1.0, 0.0, NR,
                                2.0, 3.0);
    assertTrue(PREC > Math.abs(3.0 - m2x2.prodDiag()));
    assertTrue(PREC > Math.abs(3.0 - m2x2.T().prodDiag()));

    assertTrue(PREC > Math.abs(3.0 - m2x2.det()));
    assertTrue(PREC > Math.abs(3.0 - m2x2.T().det()));
    assertTrue(PREC > Math.abs(9.0 - m2x2.mul(m2x2.T()).det()));

    Matrix m3x3 = Matrix.create(1.0, 0.0, 0.0, NR,
                                2.0, 3.0, 0.0, NR,
                                4.0, 5.0, 6.0);
    assertTrue(PREC > Math.abs(18.0 - m3x3.prodDiag()));
    assertTrue(PREC > Math.abs(18.0 - m3x3.T().prodDiag()));

    assertTrue(PREC > Math.abs(18.0 - m3x3.det()));
    assertTrue(PREC > Math.abs(18.0 - m3x3.T().det()));
    assertTrue(PREC > Math.abs(324.0 - m3x3.mul(m3x3.T()).det()));

    Matrix m5x5 = Matrix.create(1.7, 2.4, 0.1, 0.8, 1.5, NR,
                                2.3, 0.5, 0.7, 1.4, 1.6, NR,
                                0.4, 0.6, 1.3, 2.0, 2.2, NR,
                                1.0, 1.2, 1.9, 2.1, 0.3, NR,
                                1.1, 1.8, 2.5, 0.2, 0.9);
    assertTrue(PREC > Math.abs(50.7 - m5x5.det()));
    assertTrue(PREC > Math.abs(50.7 - m5x5.T().det()));

    m5x5 = Matrix.create(1.7, 2.4, 0.1, 0.8, 1.5, NR,
                         2.3, 0.5, 0.7, 1.4, 1.6, NR,
                         2.7, 3.6, 2.0, 2.9, 1.8, NR,
                         1.0, 1.2, 1.9, 2.1, 0.3, NR,
                         1.1, 1.8, 2.5, 0.2, 0.9);
    Matrix r = Matrix.rand(m5x5.rows(), m5x5.cols(), RNG);
    assertTrue(PREC > Math.abs(0.0 - m5x5.det(r)));
    assertTrue(PREC > Math.abs(0.0 - m5x5.T().det(r.setToRand(RNG))));
  }

  public void testSvd() {
    Matrix m2x2 = Matrix.create(1.0, 0.0, NR,
                                2.0, 3.0);
    Matrix U = Matrix.create(2,2);
    Matrix S = Matrix.create(2,1);
    Matrix V = Matrix.create(2,2);
    assertEquals(2, m2x2.compactSVD(U, S, V));
    assertTrue(PREC > U.mul(Matrix.diag(S)).mul(V.T()).sub(m2x2).normF());
    assertTrue(PREC > U.T().mul(U).sub(Matrix.eye(2)).normF());
    assertTrue(PREC > U.mul(U.T()).sub(Matrix.eye(2)).normF());
    assertTrue(PREC > V.T().mul(V).sub(Matrix.eye(2)).normF());
    assertTrue(PREC > V.mul(V.T()).sub(Matrix.eye(2)).normF());

    Matrix m4x2 = Matrix.create(+1, +2, NR,
                                -3, +4, NR,
                                +5, -6, NR,
                                -7, -8);
    U = Matrix.create(4,2);
    S = Matrix.create(2,1);
    V = Matrix.create(2,2);
    assertEquals(2, m4x2.compactSVD(U, S, V));
    assertTrue(PREC > U.mul(Matrix.diag(S)).mul(V.T()).sub(m4x2).normF());
    assertTrue(PREC > U.T().mul(U).sub(Matrix.eye(2)).normF());
    assertTrue(PREC > V.T().mul(V).sub(Matrix.eye(2)).normF());
    assertTrue(PREC > V.mul(V.T()).sub(Matrix.eye(2)).normF());

    Matrix m3x4 = Matrix.create(1, 2, 3, 4, NR,
                                2, 4, 7, 8, NR,
                                8, 7, 6, 5);
    U = Matrix.create(3,3);
    S = Matrix.create(3,1);
    V = Matrix.create(4,3);
    assertEquals(3, m3x4.compactSVD(U, S, V));
    assertTrue(PREC > U.mul(Matrix.diag(S)).mul(V.T()).sub(m3x4).normF());
    assertTrue(PREC > U.T().mul(U).sub(Matrix.eye(3)).normF());
    assertTrue(PREC > U.mul(U.T()).sub(Matrix.eye(3)).normF());
    assertTrue(PREC > V.T().mul(V).sub(Matrix.eye(3)).normF());

    Matrix m4x3 = Matrix.create(1, 2, 3, 4, NR,
                                2, 4, 6, 8, NR,
                                8, 7, 6, 5).T();
    U = Matrix.create(4,3);
    S = Matrix.create(3,1);
    V = Matrix.create(3,3);
    assertEquals(2, m4x3.compactSVD(U, S, V));
    assertTrue(PREC > U.mul(Matrix.diag(S)).mul(V.T()).sub(m4x3).normF());
    U = U.getMat(0, 3, 0, 1);
    S = S.getMat(0, 1, 0, 0);
    V = V.getMat(0, 2, 0, 1);
    assertTrue(PREC > U.mul(Matrix.diag(S)).mul(V.T()).sub(m4x3).normF());
    assertTrue(PREC > U.T().mul(U).sub(Matrix.eye(2)).normF());
    assertTrue(PREC > V.T().mul(V).sub(Matrix.eye(2)).normF());

    Matrix m5x9 = Matrix.create(61,   88,   2, 36,  32,  73,   99, 25,  80, NR,
                                83,   36,  36, 98,   2,  67,   51, 46,  33, NR,
                                89,  168, -14, 13, -17,  92,  153, 17, 121, NR,
                                -6, -132,  50, 85,  19, -25, -102, 29, -88, NR,
                                33,    8,  18, 59,  81,  54,   45, 33,  39);

    Matrix[] USV = m5x9.compactSVD();
    assertEquals(3, USV.length);
    U = USV[0]; S = USV[1]; V = USV[2];
    assertTrue(PREC > U.mul(Matrix.diag(S)).mul(V.T()).sub(m5x9).normF());
    assertTrue(PREC > U.mul(Matrix.diag(S)).mul(V.T()).sub(m5x9).normF());
    assertTrue(PREC > U.T().mul(U).sub(Matrix.eye(3)).normF());
    assertTrue(PREC > V.T().mul(V).sub(Matrix.eye(3)).normF());

    Matrix m5x9save = m5x9.copy();
    U = Matrix.create(5,5);
    S = Matrix.create(5,1);
    V = m5x9.T();
    assertEquals(3, m5x9.compactSVD(U, S, V));
    assertTrue(PREC > U.mul(Matrix.diag(S)).mul(V.T()).sub(m5x9save).normF());
    U = U.getMat(0, 4, 0, 2);
    S = S.getMat(0, 2, 0, 0);
    V = V.getMat(0, 8, 0, 2);
    assertTrue(PREC > U.mul(Matrix.diag(S)).mul(V.T()).sub(m5x9save).normF());
    assertTrue(PREC > U.T().mul(U).sub(Matrix.eye(3)).normF());
    assertTrue(PREC > V.T().mul(V).sub(Matrix.eye(3)).normF());
  }

  //---------------------------------------------------------------------------

  public static void main(String[] args) {
    junit.textui.TestRunner.run(MatrixTests.class);
  }
}
