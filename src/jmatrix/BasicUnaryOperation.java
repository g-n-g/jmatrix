package jmatrix;

/**
 * Implementations of basic unary operations.
 */
public enum BasicUnaryOperation implements UnaryOperation {
  ABS {
    @Override
    public double apply(double v) { return Math.abs(v); }
  },
  LOG {
    @Override
    public double apply(double v) { return Math.log(v); }
  },
  EXP {
    @Override
    public double apply(double v) { return Math.exp(v); }
  },
  SQRT {
    @Override
    public double apply(double v) { return Math.sqrt(v); }
  },
  SIN {
    @Override
    public double apply(double v) { return Math.sin(v); }
  },
  COS {
    @Override
    public double apply(double v) { return Math.cos(v); }
  },
  TAN {
    @Override
    public double apply(double v) { return Math.tan(v); }
  },
  ASIN {
    @Override
    public double apply(double v) { return Math.asin(v); }
  },
  ACOS {
    @Override
    public double apply(double v) { return Math.acos(v); }
  },
  ATAN {
    @Override
    public double apply(double v) { return Math.atan(v); }
  },
  FLOOR {
    @Override
    public double apply(double v) { return Math.floor(v); }
  },
  CEIL {
    @Override
    public double apply(double v) { return Math.ceil(v); }
  },
  ROUND {
    @Override
    public double apply(double v) { return Math.round(v); }
  },
  NEG {
    @Override
    public double apply(double v) { return -v; }
  },
  SIGN {
    @Override
    public double apply(double v) { return Math.signum(v); }
  };
}
