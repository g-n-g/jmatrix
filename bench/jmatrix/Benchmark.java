package jmatrix;

import java.util.Random;

import static jmatrix.BasicUnaryOperation.SIGN;
import static jmatrix.BasicBinaryOperation.MUL;

/** General framework for matrix computation benchmarks. */
public abstract class Benchmark
{
  private boolean debug = false;
  private double tol;

  private double time_min, time_avg, time_std, time_max;
  private double delta_max;

  Benchmark() {
    reset(0.0);
  }

  /** Resets result data. */
  public void reset(double tol) {
    this.tol = tol;

    time_min = Double.MAX_VALUE;
    time_avg = 0.0;
    time_std = 0.0;
    time_max = 0.0;
    delta_max = Double.NEGATIVE_INFINITY;
  }

  //---------------------------------------------------------------------------

  /** Returns the name of the benchmark. */
  final public String name() {
    String name = getClass().getSimpleName();
    if (name.endsWith("Benchmark")) {
      name = name.substring(0, name.length()-9);
    }
    return name;
  }

  /** Matrix input types (RG: regular, PD: positive-definite). */
  protected enum BenchmarkType {
    A_RG, A_PD, Ab_RG, Ab_PD, AB_RG, Ab_FCR
  }

  /** Chooses input type for the benchmark. */
  protected BenchmarkType type() {
    return null; // undefined
  }

  /** Benchmark specific scaling. */
  protected double customScaling(Matrix A, Matrix bB) {
    return 1.0; // no scaling
  }

  /** Performs the benchmarked computation. */
  protected abstract void compute(Matrix A, Matrix bB) throws BenchmarkException;

  /** Measures the numerical error. */
  protected double check(Matrix A, Matrix bB) throws BenchmarkException {
    return Double.NEGATIVE_INFINITY; // ignored check
  }

  //---------------------------------------------------------------------------

  public static int randInt(Random rng, int min, int max) {
    assert (max >= min);
    if (max == min) { return min; }
    return min + rng.nextInt(max-min);
  }

  /** Generates a random matrix. */
  public static Matrix randMatrix(Random rng,
                                  int minrows, int maxrows,
                                  int mincols, int maxcols,
                                  double nnzratio, double scale, double tol) {
    int Arows = randInt(rng, minrows, maxrows);
    int Acols = randInt(rng, mincols, maxcols);
    Matrix A = Matrix.randN(Arows, Acols, rng).mul(rng.nextDouble()*scale);

    if (nnzratio > 0) {
      int nnz = (int)Math.round(nnzratio*(Arows*Acols));
      for (int k = 0; k < nnz; ++k) {
        int i = rng.nextInt(Arows);
        int j = rng.nextInt(Acols);
        A.set(i, j, (2.0*rng.nextDouble()-1.0)*tol);
      }
    }
    return A;
  }

  public static Matrix randPdMatrix(Random rng,
                                    int minsize, int maxsize,
                                    double mineig, double maxeig,
                                    double tol) {
    int size = minsize + rng.nextInt(maxsize-minsize);

    Matrix M = Matrix.randN(size, size, rng);
    Matrix Q = null;
    while (true) {
      Matrix[] QR = M.QR();
      Matrix R = QR[1];
      double minabs = Double.MAX_VALUE;
      for (int i = 0; i < size; ++i) {
        minabs = Math.min(minabs, Math.abs(R.get(i,i)));
      }
      if (minabs > tol) {
        Q = QR[0].ewb(MUL, R.getDiag().ewu(SIGN).T());
        break;
      }
    }

    Matrix S = Matrix.rand(size, 1, rng).mul(maxeig-mineig).add(mineig);
    return Q.ewb(MUL, S.T()).mul(Q.T());
  }

  public static Matrix randOrthoMatrix(Random rng, int size, double tol) {
    Matrix Q = null;
    for (;;) {
      Matrix QR[] = Matrix.randN(size, size, rng).QR();
      Q = QR[0];
      boolean isFound = false;
      for (int i = 0; i < size; ++i) {
        if (Math.abs(Q.get(i,i)) < tol) {
          isFound = true;
          break;
        }
      }
      if (!isFound) { break; }
    }
    return Q;
  }

  public static Matrix randDiagMatrix(Random rng, int nrows, int ncols,
                                      double scale, double offset) {
    Matrix D = Matrix.zeros(nrows, ncols);
    for (int i = 0; i < Math.min(nrows, ncols); ++i) {
      double v = rng.nextDouble()*scale;
      double s = Math.signum(v);
      if (s == 0.0) { s = 1.0; }
      D.set(i, i, v + s*offset);
    }
    return D;
  }

  /** Runs the benchmark and saves report data to a file. */
  public void run(String[] args) throws Exception {
    int narg = 1; // first arg (classname) is ignored here

    debug = Boolean.valueOf(args[narg++]);
    if (debug) { System.err.println(args[0] + "\n"); }

    final long seed = Long.valueOf(args[narg++]);
    final int nwarmups = Integer.valueOf(args[narg++]);
    if (nwarmups <= 0) {
      throw new BenchmarkException("Invalid nwarmups: " + nwarmups + "!");
    }
    final int nrepeats = Integer.valueOf(args[narg++]);
    if (nrepeats <= 0) {
      throw new BenchmarkException("Invalid nrepeats: " + nrepeats + "!");
    }

    final int minrows = Integer.valueOf(args[narg++]);
    final int maxrows = Integer.valueOf(args[narg++]);
    if (minrows < 0 || maxrows < 0 || minrows > maxrows) {
      throw new BenchmarkException("Invalid row spec: "
                                   + minrows + ", " + maxrows + "!");
    }
    final int mincols = Integer.valueOf(args[narg++]);
    final int maxcols = Integer.valueOf(args[narg++]);
    if (mincols < 0 || maxcols < 0 || mincols > maxcols) {
      throw new BenchmarkException("Invalid col spec: "
                                   + mincols + ", " + maxcols + "!");
    }

    final double nnzratio = Double.valueOf(args[narg++]);
    if (nnzratio < 0.0 || nnzratio > 1.0) {
      throw new BenchmarkException("Invalid nnzratio: " + nnzratio + "!");
    }
    final double scale = Double.valueOf(args[narg++]);
    if (scale < 0.0) {
      throw new BenchmarkException("Invalid scale: " + nnzratio + "!");
    }
    final double tol = Double.valueOf(args[narg++]);
    if (tol < 0.0 || tol > 1.0) {
      throw new BenchmarkException("Invalid tol: " + nnzratio + "!");
    }
    final double mineig = Double.valueOf(args[narg++]);
    final double maxeig = Double.valueOf(args[narg++]);
    if (mineig < 0 || maxeig < 0.0 || mineig > maxeig) {
      throw new BenchmarkException("Invalid eig spec: "
                                   + mineig + ", " + maxeig + "!");
    }

    Random rngA, rngbB;
    {
      Random rng = new Random(seed);
      rngA = new Random(rng.nextLong());
      rngbB = new Random(rng.nextLong());
    }
    reset(tol);

    int minsize = (int)Math.ceil(Math.sqrt(minrows*mincols));
    int maxsize = (int)Math.ceil(Math.sqrt(maxrows*maxcols));

    Matrix A = null, Acopy = null;
    Matrix bB = null, bBcopy = null;
    for (int r = 0; r < nwarmups+nrepeats; ++r) {
      switch (type()) {
      case A_RG:
        A = randMatrix(rngA,
                       minrows, maxrows, mincols, maxcols,
                       nnzratio, scale, tol);
        break;
      case A_PD:
        A = randPdMatrix(rngA,
                         minsize, maxsize,
                         mineig, maxeig, tol);
        break;
      case Ab_RG:
        A = randMatrix(rngA,
                       minrows, maxrows, mincols, maxcols,
                       nnzratio, scale, tol);
        bB = randMatrix(rngbB,
                        A.rows(), A.rows(), 1, 1,
                        0.0, scale, tol);
        break;
      case Ab_PD:
        A = randPdMatrix(rngA,
                         minsize, maxsize,
                         mineig, maxeig, tol);
        bB = randMatrix(rngbB,
                        A.rows(), A.rows(), 1, 1,
                        0.0, scale, tol);
        break;
      case AB_RG:
        A = randMatrix(rngA,
                       minrows, maxrows, mincols, maxcols,
                       nnzratio, scale, tol);
        bB = randMatrix(rngbB,
                        A.rows(), A.rows(), A.cols(), A.cols(),
                        nnzratio, scale, tol);
        break;
      case Ab_FCR: {
        int nrows = randInt(rngA, minrows, maxrows);
        int ncols = randInt(rngA, 1, nrows);
        A = randOrthoMatrix(rngA, nrows, tol)
          .mul(randDiagMatrix(rngA, nrows, ncols, scale, 2*tol))
          .mul(randOrthoMatrix(rngA, ncols, tol));
        bB = randMatrix(rngbB,
                        A.rows(), A.rows(), 1, 1,
                        0.0, scale, tol);
        break;
      }
      default:
        throw new BenchmarkException("Unhandled input type: " + type() + "!");
      }
      double cs = customScaling(A, bB);
      if (cs != 1.0) {
        if (A != null) { A.mul(cs, A); }
        if (bB != null) { bB.mul(cs, bB); }
      }
      if (A != null) { Acopy = A.copy(); }
      if (bB != null) { bBcopy = bB.copy(); }

      try { Thread.sleep(1); } catch (Exception e) {}

      if (debug) { System.err.println("COMPUTE " + r); }
      double ts = System.nanoTime();
      compute(A, bB);
      ts = (System.nanoTime() - ts) / 1e6; // ms
      if (debug) { System.err.println("DONE"); }

      if (debug) { System.err.println("CHECK"); }
      double delta = check(Acopy, bBcopy);
      if (delta > delta_max) { delta_max = delta; }
      if (debug) { System.err.println("DONE\n"); }

      if (r >= nwarmups) { // ignore time measurements for warmup rounds
        if (ts < time_min) { time_min = ts; }
        if (ts > time_max) { time_max = ts; }
        time_avg += ts / nrepeats;
        time_std += ts*ts / nrepeats;
      }
    }
    time_std = Math.sqrt(time_std - time_avg*time_avg);

    System.out.println(new BenchmarkData(name(),
                                         (int)Math.ceil(time_min),
                                         (int)Math.ceil(time_max),
                                         (int)Math.ceil(time_avg),
                                         (int)Math.ceil(time_std),
                                         delta_max));
    System.out.flush();
  }

  //---------------------------------------------------------------------------

  protected void checkDelta(double delta)
    throws BenchmarkException {
    if (delta > tol) {
      throw new BenchmarkException("Too large delta: "
                                   + delta + " > " + tol + "!");
    }
  }

  /** Checks that two matrices have the same size. */
  protected void checkMatrixSizeEqual(Matrix A, Matrix B)
    throws BenchmarkException {
    if (A.rows() != B.rows()) {
      throw new BenchmarkException("Rows are not equal: "
                                   + A.rows() + " != " + B.rows() + "!");
    }
    if (A.cols() != B.cols()) {
      throw new BenchmarkException("Columns are not equal: "
                                   + A.cols() + " != " + B.cols() + "!");
    }
  }

  /** Checks that matrix is square. */
  protected void checkMatrixSquareSize(Matrix A)
    throws BenchmarkException {
    if (A.rows() != A.cols()) {
      throw new BenchmarkException("Matrix is not square: "
                                   + A.rows() + " != " + A.cols());
    }
  }

  /** Checks that two matrices are equal. */
  protected double checkMatrixEquals(Matrix A, Matrix B)
    throws BenchmarkException {
    checkMatrixSizeEqual(A, B);
    double delta = 0.0;
    for (int i = 0; i < A.rows(); ++i) {
      for (int j = 0; j < A.cols(); ++j) {
        delta = Math.max(delta, Math.abs(A.get(i,j)-B.get(i,j)));
        checkDelta(delta);
      }
    }
    return delta;
  }

  /** Checks that matrix is lower triangular. */
  protected double checkMatrixLT(Matrix A)
    throws BenchmarkException {
    double delta = 0.0;
    for (int i = 0; i < A.rows(); ++i) {
      for (int j = i+1; j < A.cols(); ++j) {
        delta = Math.max(delta, Math.abs(A.get(i,j)));
        checkDelta(delta);
      }
    }
    return delta;
  }

  /** Checks that matrix is unit lower triangular. */
  protected double checkMatrixUnitLT(Matrix A)
  throws BenchmarkException {
    double delta = checkMatrixLT(A);
    for (int i = 0; i < Math.min(A.rows(), A.cols()); ++i) {
      delta = Math.max(delta, Math.abs(1.0 - A.get(i,i)));
      checkDelta(delta);
    }
    return delta;
  }

  /** Checks that matrix is upper triangular. */
  protected double checkMatrixUT(Matrix A)
    throws BenchmarkException {
    return checkMatrixLT(A.T());
  }

  /** Checks that matrix is unit upper triangular. */
  protected double checkMatrixUnitUT(Matrix A)
    throws BenchmarkException {
    return checkMatrixUnitLT(A.T());
  }

  /** Checks that matrix is diagonal. */
  protected double checkMatrixDiag(Matrix A)
    throws BenchmarkException {
    return Math.max(checkMatrixLT(A),
                    checkMatrixUT(A));
  }

  /** Checks that matrix is zero. */
  protected double checkMatrixZero(Matrix A)
    throws BenchmarkException {
    double delta = 0.0;
    for (int i = 0; i < A.rows(); ++i) {
      for (int j = 0; j < A.cols(); ++j) {
        delta = Math.max(delta, Math.abs(A.get(i,j)));
        checkDelta(delta);
      }
    }
    return delta;
  }

  /** Checks that matrix is identity. */
  protected double checkMatrixEye(Matrix A)
    throws BenchmarkException {
    checkMatrixSquareSize(A);
    return Math.max(checkMatrixUnitLT(A),
                    checkMatrixUT(A));
  }

  /** Checks that matrix is square diagonal with zero or one elements only. */
  protected double checkMatrixEyeWithZeros(Matrix A)
    throws BenchmarkException {
    checkMatrixSquareSize(A);
    double delta = checkMatrixDiag(A);
    for (int i = 0; i < A.rows(); ++i) {
      double Aii = A.get(i,i);
      delta = Math.max(delta, Math.min(Math.abs(Aii), Math.abs(1.0-Aii)));
      checkDelta(delta);
    }
    return delta;
  }

  /** Checks that matrix has orthogonal rows. */
  protected double checkMatrixOrthoRows(Matrix A)
    throws BenchmarkException {
    return checkMatrixEye(A.mul(A.T()));
  }

  /** Checks that matrix has orthogonal columns. */
  protected double checkMatrixOrthoCols(Matrix A)
    throws BenchmarkException {
    return checkMatrixEye(A.T().mul(A));
  }

  /** Checks that matrix is orthogonal. */
  protected double checkMatrixOrtho(Matrix A)
    throws BenchmarkException {
    return Math.max(checkMatrixOrthoRows(A),
                    checkMatrixOrthoCols(A));
  }

  /** Checks that matrix has orthogonal or zero rows. */
  protected double checkMatrixOrthoOrZeroRows(Matrix A)
    throws BenchmarkException {
    return checkMatrixEyeWithZeros(A.mul(A.T()));
  }

  /** Checks that matrix has orthogonal or zero columns. */
  protected double checkMatrixOrthoOrZeroCols(Matrix A)
    throws BenchmarkException {
    return checkMatrixEyeWithZeros(A.T().mul(A));
  }

  /** Checks that matrix has orthogonal or zero rows and columns. */
  protected double checkMatrixOrthoOrZero(Matrix A)
    throws BenchmarkException {
    return Math.max(checkMatrixOrthoOrZeroRows(A),
                    checkMatrixOrthoOrZeroRows(A));
  }

  /** Checks that matrix is symmetric. */
  protected double checkMatrixSymmetric(Matrix A)
    throws BenchmarkException {
    checkMatrixSquareSize(A);
    double delta = 0.0;
    for (int i = 0; i < A.rows(); ++i) {
      for (int j = 0; j < i; ++j) {
        delta = Math.max(delta, Math.abs(A.get(i,j)-A.get(j,i)));
        checkDelta(delta);
      }
    }
    return delta;
  }

  //---------------------------------------------------------------------------

  public static void main(String[] args) {
    try {
      ((Benchmark)Class.forName(args[0]).newInstance()).run(args);
    }
    catch (Exception e) {
      System.out.println(new BenchmarkData(args[0], 0, 0, 0, 0, Double.NaN));
      System.out.flush();
      e.printStackTrace();
    }
  }
}
