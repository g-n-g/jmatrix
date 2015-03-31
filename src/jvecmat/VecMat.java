package jvecmat;

import java.util.Random;

/**
 * Shared interface of vectors and matrices.
 */
public interface VecMat {

  /**
   * Returns <code>true</code> if there is a NaN element.
   *
   * @return <code>true</code> if there is a NaN element
   */
  boolean hasNaN();

  /**
   * Returns <code>true</code> if there is an infinite element.
   *
   * @return <code>true</code> if there is an infinite element
   */
  boolean hasInf();

  /**
   * Replaces the NaN and infinite elements by the specified values.
   *
   * @param nan replacement value for NaN elements
   * @param negInf replacement value for negative infinity elements
   * @param posInf replacement value for positive infinity elements
   */
  void replaceNaNandInf(double nan, double negInf, double posInf);

  //----------------------------------------------------------------------------

  /**
   * Returns a new copy.
   *
   * @return a copy of the object
   */
  VecMat copy();

  //----------------------------------------------------------------------------

  /**
   * Sets all elements to <code>c</code>.
   *
   * @param c the new value for all elements
   * @return <code>this</code> object
   */
  VecMat setToConstant(double c);

  /**
   * Sets all elements to zero.
   *
   * @return <code>this</code> object
   */
  VecMat setToZero();

  /**
   * Sets all elements to one.
   *
   * @return <code>this</code> object
   */
  VecMat setToOne();

  /**
   * Sets all elements randomly
   * drawing the new values from the uniform distribution on [0,1].
   *
   * @param rng random number generator
   * @return <code>this</code> object
   */
  VecMat setToRand(Random rng);

  /**
   * Sets all elements randomly
   * drawing the new values from the standard normal distribution.
   *
   * @param rng random number generator
   * @return <code>this</code> object
   */
  VecMat setToRandN(Random rng);

  //----------------------------------------------------------------------------

  /**
   * Entrywise absolute value operation.
   *
   * @return entrywise absolute value
   */
  VecMat abs();

  /**
   * Sets <code>this</code> to its entrywise absolute value.
   *
   * @return <code>this</code> object
   */
  VecMat absL();

  //----------------------------------------------------------------------------

  /**
   * Entrywise sign operation with zero replacement.
   *
   * @param zeroReplacement value replacing 0.0 values
   * @return entrywise sign value
   */
  VecMat sign(double zeroReplacement);

  /**
   * Entrywise sign operation.
   *
   * @return entrywise sign value
   */
  VecMat sign();

  /**
   * Sets <code>this</code> to its entrywise sign value with zero replacement.
   *
   * @param zeroReplacement value replacing 0.0 values
   * @return <code>this</code> object
   */
  VecMat signL(double zeroReplacement);

  /**
   * Sets <code>this</code> to its entrywise sign value.
   *
   * @return <code>this</code> object
   */
  VecMat signL();
}
