package ai.gandg.jmatrix;


/**
 * Interface of binary operations.
 */
public interface BinaryOperation {

  /**
   * Applies the binary operation to values <code>x</code> and <code>y</code>.
   *
   * @param x first argument
   * @param y second argument
   * @return the result of the operation
   */
  double apply(double x, double y);
}
