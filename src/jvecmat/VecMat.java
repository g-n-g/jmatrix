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
   * Entrywise absolute value (in new object).
   *
   * @return entrywise absolute value
   */
  VecMat abs();

  /**
   * Entrywise absolute value (in place).
   *
   * @return <code>this</code> object
   */
  VecMat absL();

  //----------------------------------------------------------------------------

  /**
   * Entrywise sign operation with zero replacement (in new object).
   *
   * @param zeroReplacement value replacing 0.0 values
   * @return entrywise sign value
   */
  VecMat sign(double zeroReplacement);

  /**
   * Entrywise sign operation (in new object).
   *
   * @return entrywise sign values
   */
  VecMat sign();

  /**
   * Entrywise sign operation with zero replacement (in place).
   *
   * @param zeroReplacement value replacing 0.0 values
   * @return <code>this</code> object
   */
  VecMat signL(double zeroReplacement);

  /**
   * Entrywise sign operation (in place).
   *
   * @return <code>this</code> object
   */
  VecMat signL();

  //----------------------------------------------------------------------------

  /**
   * Entrywise negation (in new object).
   *
   * @return object with the negated elements
   */
  VecMat neg();

  /**
   * Entrywise negation (in place).
   *
   * @return <code>this</code> object with the negated elements
   */
  VecMat negL();

  //----------------------------------------------------------------------------

  /**
   * Constant division (in new object). Divides all elements of
   * <code>this</code> by constant <code>c</code>.
   *
   * @param c constant to divide with
   * @return new object with values of
   *         <code>this</code> divided by <code>c</code>
   */
  VecMat div(double c);

  /**
   * Constant division (in place). Divides all elements of
   * <code>this</code> by constant <code>c</code>.
   *
   * @param c constant to divide with
   * @return <code>this</code>
   *         by its values having divided by <code>c</code>
   */
  VecMat divL(double c);

  //----------------------------------------------------------------------------

  /**
   * Constant multipication (in new object). Multiplies all elements of
   * <code>this</code> by constant <code>c</code>.
   *
   * @param c constant to divide with
   * @return new object with values of
   *         <code>this</code> multiplied by <code>c</code>
   */
  VecMat mul(double c);

  /**
   * Constant multipication (in place). Multiplies all elements of
   * <code>this</code> by constant <code>c</code>.
   *
   * @param c constant to divide with
   * @return <code>this</code>
   *         by its values having multiplied by <code>c</code>
   */
  VecMat mulL(double c);

  //----------------------------------------------------------------------------

  /**
   * Entrywise (signed) remainder with respect to modulus <code>m</code>
   * (in new object).
   *
   * @param m modulus
   * @return new object with the remainders of the values of
   *         <code>this</code> with respect to modulus <code>m</code>
   */
  VecMat mod(double m);

  /**
   * Entrywise (signed) remainder with respect to modulus <code>m</code>
   * (in place).
   *
   * @param m modulus
   * @return <code>this</code> by taking the remainders of its values
   *         with respect to modulus <code>m</code>
   */
  VecMat modL(double m);
}
