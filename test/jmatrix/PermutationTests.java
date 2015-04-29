package jmatrix;

import java.util.Random;

/**
 * Tests for permuting vectors and matrices.
 */
public class PermutationTests extends AssertionBaseTest {

  public static final Random RNG = new Random();

  //---------------------------------------------------------------------------

  public PermutationTests(String name) {
    super(name);
  }

  //---------------------------------------------------------------------------

  public void testEmptyPerm() {
    Permutation.create(new int[]{});
  }

  public void testVectorPerms() {
    Matrix v1 = Matrix.create(1.0, 2.0, 3.0);
    Matrix v2 = Matrix.create(3.0, 1.0, 2.0);

    // test this form of create
    Permutation p = Permutation.create(new int[]{2, 0, 1});
    assertEquals(0.0, v1.mul(p).sub(v2).norm1());
    assertEquals(0.0, p.mul(v1).sub(v2).norm1());
  }

  public void testPermInv() {
    Matrix v1 = Matrix.create(1.0, 2.0, 3.0, 4.0);
    Matrix v2 = Matrix.create(3.0, 1.0, 4.0, 2.0);

    Permutation p = Permutation.create(2, 0, 3, 1);
    assertEquals(0.0, v1.mul(p).sub(v2).norm1());
    assertEquals(0.0, p.mul(v1).sub(v2).norm1());

    Permutation invp = p.inv();
    assertEquals(0.0, v2.mul(invp).sub(v1).norm1());
    assertEquals(0.0, invp.mul(v2).sub(v1).norm1());
  }

  public void testMatrixRowPerms() {
    Matrix m1 = Matrix.create(new double[][]{
        new double[]{1.0, 2.0},
        new double[]{3.0, 4.0},
        new double[]{5.0, 6.0},
        new double[]{7.0, 8.0}
      });
    Matrix m2 = Matrix.create(new double[][]{
        new double[]{5.0, 6.0},
        new double[]{3.0, 4.0},
        new double[]{7.0, 8.0},
        new double[]{1.0, 2.0}
      });

    Permutation p = Permutation.create(new int[]{2, 1, 3, 0});
    assertEquals(0.0, p.mul(m1).sub(m2).norm1());
    assertEquals(0.0, p.toMatrix().mul(m1).sub(m2).norm1());
  }

  public void testMatrixColumnPerms() {
    Matrix m1 = Matrix.create(new double[][]{
        new double[]{1.1, 2.2, 3.3, 4.4, 5.5},
        new double[]{6.6, 7.7, 8.8, 9.9, 0.0}
      });
    Matrix m2 = Matrix.create(new double[][]{
        new double[]{2.2, 3.3, 5.5, 4.4, 1.1},
        new double[]{7.7, 8.8, 0.0, 9.9, 6.6}
      });

    Permutation p = Permutation.create(new int[]{1, 2, 4, 3, 0});
    assertEquals(0.0, m1.mul(p).sub(m2).norm1());
    assertEquals(0.0, m1.mul(p.toMatrix().T()).sub(m2).norm1());
  }

  public void testPermSwaps() {
    Matrix v1 = Matrix.create(1.0, 2.0, 3.0);
    Matrix v2 = Matrix.create(3.0, 1.0, 2.0);

    Permutation p = Permutation.eye(3).swap(1, 2);
    p.swap(0, 1); // modifies p
    assertEquals(0.0, v1.mul(p).sub(v2).norm1());
  }

  public void testRandPerm() {
    Permutation p = Permutation.rand(5, RNG);
    p = Permutation.create(p.array()); // check validity

    Matrix M = Matrix.rand(5, 3, RNG);
    assertEquals(0.0, p.mul(p.inv().mul(M)).sub(M).norm1());
  }
}
