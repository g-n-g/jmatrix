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
   * Entrywise negation (in new object).
   *
   * @return object with the negated elements
   */
  VecMat neg();

  /**
   * Constant addition (in new object). Adds constant <code>c</code> to all
   * elements of <code>this</code> object.
   *
   * @param c constant to add
   * @return new object with values of
   *         <code>this</code> shifted by <code>c</code>
   */
  VecMat add(double c);

  /**
   * Constant subtraction (in new object). Subtracts constant <code>c</code>
   * from all elements of <code>this</code> object.
   *
   * @param c constant to subtract
   * @return new object with values of
   *         <code>this</code> shifted by <code>-c</code>
   */
  VecMat sub(double c);

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
   * Constant multipication (in new object). Multiplies all elements of
   * <code>this</code> by constant <code>c</code>.
   *
   * @param c constant to divide with
   * @return new object with values of
   *         <code>this</code> multiplied by <code>c</code>
   */
  VecMat mul(double c);

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
   * Takes the reciproc of all elements (in new object).
   *
   * @return elementwise reciproc in a new object
   */
  VecMat reciproc();
}
