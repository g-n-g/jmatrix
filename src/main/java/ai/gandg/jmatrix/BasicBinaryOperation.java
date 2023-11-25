package ai.gandg.jmatrix;

/**
 * Implementation of basic binary operations.
 */
public enum BasicBinaryOperation implements BinaryOperation {

  /** Addition of two scalars. */
  ADD {
    @Override
    public double apply(double x, double y) { return x + y; }
  },

  /** Subtracting a scalar from another. */
  SUB {
    @Override
    public double apply(double x, double y) { return x - y; }
  },

  /** Multiplying two scalars. */
  MUL {
    @Override
    public double apply(double x, double y) { return x * y; }
  },

  /** Dividing a scalar by another. */
  DIV {
    @Override
    public double apply(double x, double y) { return x / y; }
  },

  /** Taking the modulo of a scalar by another. */
  MOD {
    @Override
    public double apply(double x, double y) { return x % y; }
  },

  /** Taking the minimum of two scalars. */
  MIN {
    @Override
    public double apply(double x, double y) { return Math.min(x,y); }
  },

  /** Taking the maximum of two scalars. */
  MAX {
    @Override
    public double apply(double x, double y) { return Math.max(x,y); }
  },

  /** Raising a scalar to the power of another. */
  POW {
    @Override
    public double apply(double x, double y) { return Math.pow(x,y); }
  };
}
