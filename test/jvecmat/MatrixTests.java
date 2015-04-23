package jvecmat;

import java.util.Random;
import static jvecmat.Matrix.NR;

/**
 * Tests for the Matrix class.
 */
public class MatrixTests extends AssertionBaseTest {

  public static final double PREC = 1e-8;
  public static final Random RNG = new Random();

  //---------------------------------------------------------------------------

  public MatrixTests(String name) { super(name); }

  //---------------------------------------------------------------------------

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

    Vector vr = Vector.ones(m);
    assertTrue(PREC > m1.addRow(vr).sub(m2).normF());
    assertTrue(PREC > m2.subRow(vr).normF());

    Vector vc = Vector.ones(n);
    assertTrue(PREC > m1.addCol(vc).sub(m2).normF());
    assertTrue(PREC > m2.subCol(vc).normF());

    vr = Vector.scalars(m, 2.0);
    vc = Vector.scalars(n, 4.0);
    assertTrue(PREC > m2.divRow(vr).sub(0.50).normF());
    assertTrue(PREC > m2.divCol(vc).sub(0.25).normF());
    assertTrue(PREC > m2.mulRow(vr).sub(2.0).normF());
    assertTrue(PREC > m2.mulCol(vc).sub(4.0).normF());
  }

  public void testAbsAndSignAndNeg() {
    Matrix m = Matrix.randN(4, 6, RNG);
    assertTrue(PREC > m.sub(m.sign().emul(m.abs())).norm1());
    assertTrue(PREC > m.neg().sub(m.neg().sign().emul(m.abs())).norm1());

    m.neg(m);
    assertTrue(PREC > m.sub(m.sign().emul(m.abs())).norm1());

    Matrix s = Matrix.zeros(4, 6);
    Matrix a = Matrix.zeros(4, 6);
    m.copy(s).sign(s);
    m.copy(a).abs(m);
    assertTrue(PREC > m.sub(s.emul(a)).norm1());

    s.copy(a).sign(s);
    assertTrue(PREC > s.sub(a).norm1());

    s.set(2, 2, 0.0);
    a.set(2, 2, 0.0);
    assertTrue(PREC > s.sub(a).norm1());
    s.sign(8.8, s);
    assertFalse(PREC > s.sub(a).norm1());
  }

  public void testMod() {
    double m = 2.3;
    Matrix mm = Matrix.rand(5, 4, RNG);
    mm.mul(8.0, mm);
    Matrix oo = Matrix.scalars(5, 4, m);

    assertTrue(PREC > oo.mod(m).norm1());
    assertTrue(PREC > mm.add(oo.mul(2.0)).mod(m).sub(mm.mod(m)).norm1());

    mm.mul(-1.0, mm);
    assertTrue(PREC > mm.add(oo.mul(-2.0)).mod(m).sub(mm.mod(m)).norm1());

    oo.mod(m, oo);
    assertTrue(PREC > oo.norm1());
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

  public void testTransposeAndSomeGetSet() {
    Vector v1 = new Vector(new double[]{1.0, 2.0, 3.0});
    Vector v2 = new Vector(new double[]{4.5, 5.5, 6.5});
    Vector v3 = new Vector(new double[]{6.1, 2.2, -5.1});
    Vector v4 = new Vector(new double[]{2.9, -1.1, 0.3});

    Matrix m1 = Matrix.createByRows(new Vector[]{v1, v2, v3, v4});
    assertEquals(4, m1.rows());
    assertEquals(v1.length(), m1.cols());

    Matrix m2 = Matrix.createByCols(new Vector[]{v1, v2, v3, v4});
    assertEquals(v1.length(), m2.rows());
    assertEquals(4, m2.cols());

    assertTrue(PREC > m1.sub(m2.T()).norm1());
    assertTrue(PREC > m2.sub(m1.T()).norm1());

    assertTrue(PREC > m1.getRow(3).sub(v4).norm2());
    assertTrue(PREC > m2.getCol(3).sub(v4).norm2());

    Matrix m3 = m1.getMat(2, 3, 1, 2);
    Matrix m4 = m2.getMat(1, 2, 2, 3);
    assertTrue(PREC > m3.sub(m4.T()).normF());

    Matrix A = Matrix.rand(7, 7, RNG);
    Matrix B = Matrix.rand(7, 7, RNG);

    assertEquals(A.mul(B).get(2, 5),
                 A.getRow(2).inner(B.getCol(5)),
                 PREC);

    Vector row = Vector.create(7);
    Vector col = Vector.create(7);
    assertEquals(A.mul(B).get(2, 5),
                 A.getRow(2, row).inner(B.getCol(5, col)),
                 PREC);

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
    Vector a = Vector.create(new double[]{1, 5, 9});
    assertEquals(0.0, A.getDiag().sub(a).normI());
    assertEquals(0.0, A.T().getDiag().sub(a).normI());

    Matrix B = Matrix.create(1, 2, 3, 4, 5, NR,
                             6, 7, 8, 9, 0);
    Vector b = Vector.create(new double[]{1, 7});
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
    Vector v1 = Vector.ones(m);
    Vector v2 = Vector.unit(m, m/2);
    Matrix m1 = Matrix.ones(n, m);
    Matrix m2 = Matrix.eye(m);

    assertTrue(PREC > v1.sub(m2.mul(v1)).norm1());
    assertTrue(PREC > Vector.ones(n).sub(m1.mul(v2)).norm1());

    for (int i = 0; i < m; ++i) v2.set(i, i+1);
    double sum = m * (m+1) / 2.0;
    assertTrue(PREC > Vector.ones(n).mul(sum).sub(m1.mul(v2)).norm1());
    assertTrue(PREC > v2.sub(m2.mul(v2)).normI());
    assertTrue(PREC < v2.sub(Matrix.ones(m, m).mul(v2)).normI());

    assertTrue(PREC > v2.emul(v2).sub(Matrix.diag(v2).mul(v2)).norm2());

    assertTrue(PREC > m1.mul(v1.add(v2)).sub(m1.mul(v1).add(m1.mul(v2)))
                        .norm1());
  }

  public void testMatrixDiagProduct() {
    final int n = 3, m = 4;
    Matrix m1 = Matrix.ones(n, m);
    Vector v = new Vector(new double[]{2.0, 3.0, 5.0, 7.0});
    Matrix m2 = Matrix.createByRows(new Vector[]{v, v, v});

    Matrix m3 = m1.mulD(v);
    assertTrue(PREC > m2.sub(m3).norm1());

    m3.mulD(v, m3);
    assertTrue(PREC > m2.emul(m2).sub(m3).norm1());
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

  public void testRowColSums() {
    Matrix m = Matrix.create(17, 24,  1,  8, 15, 1, NR,
                             23,  5,  7, 14, 16, 1, NR,
                              4,  6, 31, 20, 22, 1, NR,
                             10, 12, 19, 21,  3, 1, NR,
                             11, 18, 25,  2,  9, 1);

    Vector vRow = Vector.create(new double[]{65, 65, 83, 65, 65, 5});
    assertTrue(PREC > m.rowSum().sub(vRow).norm1());

    Vector vCol = Vector.create(new double[]{66, 66, 84, 66, 66});
    assertTrue(PREC > m.colSum().sub(vCol).norm1());
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
  }

  public void testCholeskyDecompositionLD() {
    Matrix m = Matrix.create(2.0, 1.0, 1.0, NR,
                             1.0, 3.0, 1.0, NR,
                             1.0, 1.0, 2.0);
    Matrix LD[] = m.choleskyLD();
    Matrix L = LD[0];
    Vector D = LD[1].getDiag();
    assertTrue(PREC > m.sub(L.mulD(D).mul(L.T())).norm1());
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
    M3.QR(Q, R, Vector.create(3));
    assertTrue(PREC > M3.sub(Q.mul(R)).norm1());
    assertTrue(PREC > Matrix.eye(4).sub(Q.T().mul(Q)).norm1());

    Matrix M4 = Matrix.create(+1, +2, NR, // size: 4x2
                              -0, +0, NR,
                              +0, -6, NR,
                              -7, -8);
    Q = Matrix.rand(4, 4, RNG);
    R = Matrix.rand(4, 2, RNG);
    M4.QR(Q, R, Vector.rand(3, RNG));
    assertTrue(PREC > M4.sub(Q.mul(R)).norm1());
    assertTrue(PREC > Matrix.eye(4).sub(Q.T().mul(Q)).norm1());

    Matrix M5 = Matrix.create(  0,   0, NR, // size: 3x2
                               10,  20, NR,
                              200, 100);
    Q = Matrix.create(3, 3);
    R = Matrix.create(3, 2);
    M5.QR(Q, R, Vector.create(10));
    assertTrue(PREC > M5.sub(Q.mul(R)).norm1());
    assertTrue(PREC > Matrix.eye(3).sub(Q.T().mul(Q)).norm1());

    Matrix M6 = Matrix.create(0, 10, NR, // size: 3x2
                              0, 20, NR,
                              0, 30);
    Q = Matrix.create(3, 3);
    R = Matrix.create(3, 2);
    M6.QR(Q, R, Vector.create(10));
    assertTrue(PREC > M6.sub(Q.mul(R)).norm1());
    assertTrue(PREC > Matrix.eye(3).sub(Q.T().mul(Q)).norm1());

    Matrix M7 = Matrix.create(1, 2, 3, 4, NR, // size: 3x4
                              2, 4, 6, 8, NR,
                              8, 7, 6, 5);
    Q = Matrix.rand(3, 3, RNG);
    R = Matrix.rand(3, 4, RNG);
    M7.QR(Q, R, Vector.rand(2, RNG));
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
    Matrix r2 = Matrix.rand(m5.rows(), m5.cols(), RNG);
    assertTrue(PREC > i5.sub(m5b.mul(m5b.invPsd(r1, r2))).norm1());
    assertTrue(PREC > i5.sub(m5b.invPsd(r2, r1).mul(m5b)).norm1());
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

  public void testReciproc() {
    final int n = 2, m = 5;
    Matrix v = Matrix.scalars(n, m, 0.1);
    Matrix o = Matrix.ones(n, m);
    Matrix r = Matrix.zeros(n, m);

    assertTrue(PREC > o.sub(v.emul(v.reciproc())).normF());

    v.reciproc(r);
    assertTrue(PREC > o.sub(v.emul(r)).normF());

    r.reciproc(r);
    assertTrue(PREC > v.sub(r).normF());
  }

  //---------------------------------------------------------------------------

  public static void main(String[] args) {
    junit.textui.TestRunner.run(MatrixTests.class);
  }
}
