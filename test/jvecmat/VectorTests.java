package jvecmat;

import java.util.Random;

/**
 * Tests for the Vector class.
 */
public class VectorTests extends AssertionBaseTest {

  public static final double PREC = 1e-8;
  public static final Random RNG = new Random();

  //---------------------------------------------------------------------------

  public VectorTests(String name) {
    super(name);
  }

  //---------------------------------------------------------------------------

  public void testNorms() {
    Vector v = Vector.create(new double[]{1.0, 3.0, -5.0, 1.0});
    assertTrue(PREC > Math.abs(10 - v.norm1()));
    assertTrue(PREC > Math.abs(6.0 - v.norm2()));
    assertTrue(PREC > Math.abs(5.0 - v.normI()));
  }

  public void testNaNandInfChecks() {
    Vector v = Vector.create(new double[]{1.0, 3.0, -5.0, 1.0});
    assertFalse(v.hasNaN());
    assertFalse(v.hasInf());

    v.set(1, Double.NaN);
    assertTrue(v.hasNaN());
    assertFalse(v.hasInf());

    v.set(2, Double.POSITIVE_INFINITY);
    assertTrue(v.hasNaN());
    assertTrue(v.hasInf());

    v.set(1, 4.4);
    v.set(2, Double.NEGATIVE_INFINITY);
    assertFalse(v.hasNaN());
    assertTrue(v.hasInf());

    v.set(1, Double.NaN);
    v.set(3, Double.POSITIVE_INFINITY);
    v.replaceNaNandInf(1.1, 2.2, 3.3);
    assertFalse(v.hasNaN());
    assertFalse(v.hasInf());
  }

  public void testBasicLinearOps() {
    final int len = 10;
    Vector v1 = Vector.zero(len);
    Vector v2 = Vector.one(len);
    Vector v3 = Vector.constant(len, 42.42);
    Vector v4 = Vector.unit(len, len/2);

    assertTrue(PREC > v2.sub(v2.add(v1)).norm2());
    assertTrue(PREC > v3.sub(v2.mul(42.42)).norm2());
    assertTrue(PREC > v2.sub(v3.div(42.42)).norm1());
    assertTrue(PREC > v1.sub(v4.sub(v4)).normI());

    assertTrue(PREC > v2.add(-1.0).norm1());
    assertTrue(PREC > v2.sub(1.0).norm2());

    v2.add(v2, v2);
    assertTrue(PREC > v2.sub(v3.div(42.42).mul(2)).norm2());
    v2.setToOne();

    v2.sub(v2.div(2), v2);
    assertTrue(PREC > v2.sub(v3.div(42.42).mul(0.5)).norm2());
    v2.setToOne();

    v2.mul(2, v2);
    assertTrue(PREC > v2.sub(v3.div(42.42).mul(2)).norm2());
    v2.setToOne();

    v2.div(8, v2);
    assertTrue(PREC > v2.sub(v3.div(8*42.42)).norm2());
    v2.setToOne();

    v2.setToConstant(-2.1);
    v2.add(2.1, v2);
    assertTrue(PREC > v2.normI());
    v2.sub(-1.0, v2);
    assertTrue(PREC > v2.sub(Vector.one(len)).normI());
    v2.setToOne();
  }

  public void testMod() {
    double m = 2.3;
    Vector v = Vector.rand(10, RNG).mul(8.0);
    Vector o = Vector.constant(10, m);

    assertTrue(PREC > o.mod(m).norm1());
    assertTrue(PREC > v.add(o.mul(2.0)).mod(m).sub(v.mod(m)).norm1());

    v.mul(-1.0, v);
    assertTrue(PREC > v.add(o.mul(-2.0)).mod(m).sub(v.mod(m)).norm1());

    o.mod(m, o);
    assertTrue(PREC > o.norm1());
  }

  public void testCopy() {
    Vector v = Vector.rand(5, RNG);
    Vector c = v.copy();

    v.add(Vector.randN(5, RNG), v);
    assertTrue(PREC > v.sub(v).norm1());
    assertTrue(PREC < c.sub(v).norm1());

    c.copy(v);
    assertTrue(PREC > v.sub(c).norm1());
  }

  public void testGetVec() {
    Vector v = Vector.create(new double[]{1.0, 2.0, 3.0, 4.0});

    Vector v1 = v.getVec(0, 1);
    assertEquals(2, v1.length());
    assertEquals(1.0, v1.get(0));
    assertEquals(2.0, v1.get(1));

    Vector v2 = v.getVec(1, 3);
    assertEquals(3, v2.length());
    assertEquals(2.0, v2.get(0));
    assertEquals(3.0, v2.get(1));
    assertEquals(4.0, v2.get(2));

    Vector u = Vector.rand(5, RNG);
    v = Vector.rand(5, RNG);

    assertTrue(PREC > u.emul(v).getVec(0, 3)
                       .sub(u.getVec(0, 3).emul(v.getVec(0, 3))).norm1());
    assertTrue(PREC > u.emul(v).getVec(1, 4)
                       .sub(u.getVec(1, 4).emul(v.getVec(1, 4))).norm1());
    assertTrue(PREC > u.emul(v).getVec(1, 3)
                       .sub(u.getVec(1, 3).emul(v.getVec(1, 3))).norm1());
  }

  public void testSetToZeroOne() {
    Vector v = Vector.one(4);
    v.setToZero();
    assertTrue(PREC > Vector.zero(4).sub(v).norm2());
    v.setToOne();
    assertTrue(PREC > Vector.one(4).sub(v).norm2());
  }

  public void testNormalize() {
    Vector v = Vector.zero(4);

    v.normalize1();
    assertTrue(PREC > v.norm1());

    v.normalize2();
    assertTrue(PREC > v.norm2());

    v.normalizeI();
    assertTrue(PREC > v.normI());

    v.setToOne();
    v.normalize1();
    assertTrue(PREC > Math.abs(1.0 - v.norm1()));

    v.setToRandN(RNG);
    v.normalize2();
    assertTrue(PREC > Math.abs(1.0 - v.norm2()));

    v.setToRand(RNG);
    v.normalizeI();
    assertTrue(PREC > Math.abs(1.0 - v.normI()));

    v = Vector.unit(4, 2);
    assertTrue(PREC > Math.abs(1.0 - v.norm1()));
    assertTrue(PREC > Math.abs(1.0 - v.norm2()));
    assertTrue(PREC > Math.abs(1.0 - v.normI()));
  }

  public void testAbsAndSignAndNeg() {
    Vector v = Vector.randN(10, RNG);        
    assertTrue(PREC > v.sub(v.sign().emul(v.abs())).norm1());
    assertTrue(PREC > v.neg().sub(v.neg().sign().emul(v.abs())).norm1());

    v.neg(v);
    assertTrue(PREC > v.sub(v.sign().emul(v.abs())).norm1());

    Vector s = Vector.zero(10);
    Vector a = Vector.zero(10);
    v.copy(s).sign(s);
    v.copy(a).abs(v);
    assertTrue(PREC > v.sub(s.emul(a)).norm1());

    s.copy(a).sign(s);
    assertTrue(PREC > s.sub(a).norm1());

    s.set(2, 0.0);
    a.set(2, 0.0);
    assertTrue(PREC > s.sub(a).norm1());
    s.sign(8.8, s);
    assertFalse(PREC > s.sub(a).norm1());
  }

  public void testInnerProduct() {
    final int len = 10;
    Vector v1 = Vector.one(len);
    Vector v2 = Vector.unit(len, len / 3);

    assertTrue(PREC > Math.abs(len - v1.inner(v1)));
    assertTrue(PREC > Math.abs(1.0 - v1.inner(v2)));
    assertTrue(PREC > Math.abs(2.0 * len - v1.mul(2).inner(v1)));
    assertTrue(PREC > Math.abs(2.0 - v2.inner(v1.mul(2))));
  }

  public void testEntrywiseMultiplication() {
    final int len = 10;
    Vector v1 = Vector.constant(len, 2);
    Vector v2 = Vector.unit(len, 2*len/3);

    assertTrue(PREC > v1.mul(2).sub(v1.emul(v1)).norm1());
    assertTrue(PREC > v2.mul(2).sub(v1.emul(v2)).norm2());

    v1.emul(v1.mul(2), v1);
    assertTrue(PREC > Vector.one(len).mul(8).sub(v1).normI());
    v1.setToConstant(2.0);
  }

  public void testVectorMatrixMultiplication() {
    final int n = 5, m = 10;
    Matrix m1 = Matrix.eye(n);
    Matrix m2 = Matrix.one(n, m);
    Vector v1 = Vector.unit(n, 3*n/5);
    Vector v2 = Vector.constant(n, 1.42);
    Vector v3 = Vector.one(m);

    assertTrue(PREC > v2.sub(v2.mul(m1)).norm1());
    assertTrue(PREC > v3.sub(v1.mul(m2)).norm1());

    v1.mul(m2, v3);
    assertTrue(PREC > Vector.one(m).sub(v3).norm1());
  }

  public void testOuterProduct() {
    final int n = 5;
    Matrix m = Matrix.constant(n, n, 9);
    Vector v = Vector.constant(n, 3);

    assertTrue(PREC > m.sub(v.outer(v)).norm1());

    v.mul(2).outer(v, m);
    assertTrue(PREC > Matrix.constant(n, n, 18).sub(m).norm1());
  }

  public void testQuadraticProduct() {
    Vector v = Vector.create(new double[]{1.1, 2.2, 3.3});
    Matrix m = Matrix.createByCols(new Vector[]{v, v, v});
    assertTrue(PREC > Math.abs(v.mulQ(m) - v.mul(m).inner(v)));

    Vector rv = Vector.rand(4, RNG);
    Matrix rm = Matrix.randN(4, 4, RNG);
    assertTrue(PREC > Math.abs(rv.mulQ(rm) - rv.mul(rm).inner(rv)));
  }

  public void testMatrixDiagProduct() {
    final int n = 3, m = 4;
    Matrix m1 = Matrix.one(n, m);
    Vector v = Vector.create(new double[]{2.0, 3.0, 5.0});
    Matrix m2 = Matrix.createByCols(new Vector[]{v, v, v, v});

    Matrix m3 = v.mulD(m1);
    assertTrue(PREC > m2.sub(m3).norm1());

    v.mulD(m3, m3);
    assertTrue(PREC > m2.emul(m2).sub(m3).norm1());
  }

  public void testReciproc() {
    final int n = 5;
    Vector v = Vector.constant(n, 0.1);
    Vector o = Vector.one(n);
    Vector r = Vector.zero(n);

    assertTrue(PREC > o.sub(v.emul(v.reciproc())).norm2());

    v.reciproc(r);
    assertTrue(PREC > o.sub(v.emul(r)).norm2());

    r.reciproc(r);
    assertTrue(PREC > v.sub(r).norm2());
  }

  //---------------------------------------------------------------------------

  public static void main(String[] args) {
    junit.textui.TestRunner.run(VectorTests.class);
  }
}
