package jmatrix;

/**
 * Implementation of basic binary operations.
 */
public enum BasicBinaryOperation implements BinaryOperation {
  ADD {
    @Override
    public double apply(double x, double y) { return x + y; }
  },
  SUB {
    @Override
    public double apply(double x, double y) { return x - y; }
  },
  MUL {
    @Override
    public double apply(double x, double y) { return x * y; }
  },
  DIV {
    @Override
    public double apply(double x, double y) { return x / y; }
  },
  MOD {
    @Override
    public double apply(double x, double y) { return x % y; }
  },
  MIN {
    @Override
    public double apply(double x, double y) { return Math.min(x,y); }
  },
  MAX {
    @Override
    public double apply(double x, double y) { return Math.max(x,y); }
  },
  POW {
    @Override
    public double apply(double x, double y) { return Math.pow(x,y); }
  };
}
