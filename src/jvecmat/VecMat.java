package jvecmat;

import java.util.Random;

/**
 * Shared interface of vectors and matrices.
 */
public interface VecMat {

  /**
   * Returns <code>true</code> if the vector/matrix has a NaN element.
   *
   * @return <code>true</code> if the vector/matrix has a NaN element
   */
  boolean hasNaN();

  /**
   * Returns <code>true</code> if the vector/matrix has an infinite element.
   *
   * @return <code>true</code> if the vector/matrix has an infinite element
   */
  boolean hasInf();

  /**
   * Replaces the NaN and infinite elements of a vector/matrix.
   *
   * @param nan replacement value for NaN elements
   * @param negInf replacement value for negative infinity elements
   * @param posInf replacement value for positive infinity elements
   */
  void replaceNaNandInf(double nan, double negInf, double posInf);

  //----------------------------------------------------------------------------

  /**
   * Returns a new copy of the vector/matrix.
   *
   * @return copy of the vector/matrix
   */
  VecMat copy();

  //----------------------------------------------------------------------------

  /**
   * Set all vector/matrix elements to <code>c</code>.
   *
   * @param c the new value for all vector/matrix elements
   * @return <code>this</code> vector/matrix
   */
  VecMat setToConstant(double c);

  /**
   * Set all vector/matrix elements to zero.
   *
   * @return <code>this</code> vector/matrix
   */
  VecMat setToZero();

  /**
   * Set all vector/matrix elements to one.
   *
   * @return <code>this</code> vector/matrix
   */
  VecMat setToOne();

  /**
   * Set all vector/matrix elements randomly
   * drawing the new values from the uniform distribution on [0,1].
   *
   * @param rng random number generator
   * @return <code>this</code> vector/matrix
   */
  VecMat setToRand(Random rng);

  /**
   * Set all vector elements randomly
   * drawing the new values from the standard normal distribution.
   *
   * @param rng random number generator
   * @return <code>this</code> vector
   */
  VecMat setToRandN(Random rng);
}
