package ai.gandg.jmatrix;


/**
 * Interface of unary operations.
 */
public interface UnaryOperation {

  /**
   * Applies the unary operation to value <code>v</code>.
   *
   * @param v value
   * @return result of the operation
   */
  double apply(double v);
}
