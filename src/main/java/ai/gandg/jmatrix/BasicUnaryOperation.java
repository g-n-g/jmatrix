package ai.gandg.jmatrix;


/**
 * Implementations of basic unary operations.
 */
public enum BasicUnaryOperation implements UnaryOperation {

  /** Taking the absolute value. */
  ABS {
    @Override
    public double apply(double v) { return Math.abs(v); }
  },

  /** Taking the natural logarithm. */
  LOG {
    @Override
    public double apply(double v) { return Math.log(v); }
  },

  /** Taking the power of Euler's number. */
  EXP {
    @Override
    public double apply(double v) { return Math.exp(v); }
  },

  /** Taking the square root. */
  SQRT {
    @Override
    public double apply(double v) { return Math.sqrt(v); }
  },

  /** Taking the trigonometric sine. */
  SIN {
    @Override
    public double apply(double v) { return Math.sin(v); }
  },

  /** Taking the trigonometric cosine. */
  COS {
    @Override
    public double apply(double v) { return Math.cos(v); }
  },

  /** Taking the trigonometric tangent. */
  TAN {
    @Override
    public double apply(double v) { return Math.tan(v); }
  },

  /** Taking the arc sine. */
  ASIN {
    @Override
    public double apply(double v) { return Math.asin(v); }
  },

  /** Taking the arc cosine. */
  ACOS {
    @Override
    public double apply(double v) { return Math.acos(v); }
  },

  /** Taking the arc tangent. */
  ATAN {
    @Override
    public double apply(double v) { return Math.atan(v); }
  },

  /** Rounding down (floor). */
  FLOOR {
    @Override
    public double apply(double v) { return Math.floor(v); }
  },

  /** Rounding up (ceil). */
  CEIL {
    @Override
    public double apply(double v) { return Math.ceil(v); }
  },

  /** Rounding. */
  ROUND {
    @Override
    public double apply(double v) { return Math.round(v); }
  },

  /** Negating. */
  NEG {
    @Override
    public double apply(double v) { return -v; }
  },

  /** Taking the reciproc. */
  RECIPROC {
    @Override
    public double apply(double v) { return 1.0 / v; }
  },

  /** Taking the sign. */
  SIGN {
    @Override
    public double apply(double v) { return Math.signum(v); }
  };
}
