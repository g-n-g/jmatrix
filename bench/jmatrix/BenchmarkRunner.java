package jmatrix;

import java.io.File;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.lang.Process;
import java.lang.ProcessBuilder;
import java.util.Random;

/** Running benchmarks. */
public class BenchmarkRunner
{
  private static final String[] BENCHMARKS = {
    "MatMulAtA",
    "MatMulAtB",
    "Norm2",
    "Determinant",
    "LU",
    "QR",
    "ReducedQR",
    "QRnoQ",
    "OrthoNorm",
    "SVD",
    "ReducedSVD",
    "CholeskyL",
    "CholeskyLD",
    "SolvePDEqnLL",
    "SolvePDEqnLU",
    "SolveLSFCREqnLL",
    "SolveLSFCREqnQRnoQ",
    "SolveLSFCREqnReducedQR",
    "MatInv",
    "MatInvPsd"
  };

  private static boolean existsBenchmark(String benchmark) {
    for (String b : BENCHMARKS) {
      if (benchmark.equals(b)) { return true; }
    }
    return false;
  }

  // Default random seed (automatically generated if set to 0).
  private static final long SEED = 0;

  // Number of warmup and time measured runs of a benchmark.
  private static final int NWARMUPS = 20;
  private static final int NREPEATS = 200;

  // Benchmark input parameters.
  private static final int MINROWS = 100, MAXROWS = 500;
  private static final int MINCOLS =  50, MAXCOLS = 200;
  private static final double NNZRATIO = 0.25;
  private static final double SCALE = 1000.0;

  // Eigenvalue range for positive definite (pd) matrices.
  private static final double MINEIG = 0.0001;
  private static final double MAXEIG = 1000.0;

  // Error tolerance triggering an abort by an exception.
  private static final double TOL = 5e-6;

  // Turns matrix statistics report on/off.
  private static final boolean MATRIX_STAT = false;

  // Turns benchmark debugging on/off.
  private static final boolean DEBUG = true;
  // Anything printed on standard error within the benchmarks
  // will be written into the debug files named with this prefix.
  private static final String DEBUG_FILE = "JMATRIX_BENCHMARK_ERROR_";

  // Java command and classpath.
  private static final String JAVA_CMD = "java";
  private static final String CP = "bin:bin" + File.separator + "bench";

  //---------------------------------------------------------------------------

  private static int countDec(long n) {
    return (int)Math.ceil(Math.log10((double)n+1));
  }

  private static int maximum(int... values) {
    int maxi = Integer.MIN_VALUE;
    for (int i : values) {
      if (i > maxi) { maxi = i; }
    }
    return maxi;
  }

  public void run(long seed, String[] benchmarks) throws Exception {
    // determine random seed
    if (0 == seed) {
      seed = new Random().nextLong();
    }
    System.out.println("Random seed: " + seed);
    System.out.println();

    // gather and report matrix statistics
    if (MATRIX_STAT) { reportMatrixStatistics(seed); }

    // run the benchmarks
    double ts = System.nanoTime();
    System.out.println("Running benchmarks...");
    System.out.format("  %25s : %9s\n", "", "start (s)");
    System.out.flush();
    BenchmarkData[] data = new BenchmarkData[BENCHMARKS.length];
    for (int b = 0; b < benchmarks.length; ++b) {
      data[b] = null;
      
      String benchmark = benchmarks[b];
      if (!existsBenchmark(benchmark)) {
        System.out.println("  " + benchmark + " does not exist!");
        continue;
      }
      System.out.format("  %25s : %9.1f\n", benchmark, (System.nanoTime()-ts) / 1e9);
      System.out.flush();
      
      String classname = benchmark + "Benchmark";
      ProcessBuilder pb = new ProcessBuilder(JAVA_CMD,
                                             "-d64", "-Xms512m", "-Xmx512m",
                                             "-cp", CP,
                                             "jmatrix.Benchmark",
                                             "jmatrix." + classname,
                                             "" + DEBUG,
                                             "" + seed,
                                             "" + NWARMUPS,
                                             "" + NREPEATS,
                                             "" + MINROWS, "" + MAXROWS,
                                             "" + MINCOLS, "" + MAXCOLS,
                                             "" + NNZRATIO,
                                             "" + SCALE,
                                             "" + TOL,
                                             "" + MINEIG,
                                             "" + MAXEIG);
      if (DEBUG) { pb.redirectError(new File(DEBUG_FILE + b)); }

      Process proc = pb.start();
      BufferedReader in =
        new BufferedReader(new InputStreamReader(proc.getInputStream()));

      String line = in.readLine();
      if (line == null) {
        throw new BenchmarkException("Cannot read benchmark data from " +
                                     classname + "!");
      }
      data[b] = BenchmarkData.valueOf(line);
      proc.waitFor();
    }
    System.out.println();
    System.out.flush();

    // report the results
    reportBenchmarkResults(data);
  }

  private void reportMatrixStatistics(long seed) {
    Random rng = new Random(seed);

    // regular matrices

    int rowmin = Integer.MAX_VALUE, rowmax = 0; double rowavg = 0, rowstd = 0;
    int colmin = Integer.MAX_VALUE, colmax = 0; double colavg = 0, colstd = 0;
    int sizmin = Integer.MAX_VALUE, sizmax = 0; double sizavg = 0, sizstd = 0;
    int nnzmin = Integer.MAX_VALUE, nnzmax = 0; double nnzavg = 0, nnzstd = 0;

    System.out.println("Matrix statistics:");
    System.out.println();

    final int nrepeats = NWARMUPS+NREPEATS;
    for (int r = 0; r < nrepeats; ++r) {
      Matrix A = Benchmark.randMatrix(rng,
                                      MINROWS, MAXROWS,
                                      MINCOLS, MAXCOLS,
                                      NNZRATIO, SCALE, TOL);

      int rows = A.rows(), cols = A.cols();
      int size = rows * cols;
      int nnz = 0;
      for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
          double v = Math.abs(A.get(i,j));
          if (v < TOL) { ++nnz; }
        }
      }

      rowmin = Math.min(rowmin, rows);
      rowmax = Math.max(rowmax, rows);
      rowavg += rows / nrepeats;
      rowstd += rows*rows / nrepeats;

      colmin = Math.min(colmin, cols);
      colmax = Math.max(colmax, cols);
      colavg += cols / nrepeats;
      colstd += cols*cols / nrepeats;

      sizmin = Math.min(sizmin, size);
      sizmax = Math.max(sizmax, size);
      sizavg += size / nrepeats;
      sizstd += size*size / nrepeats;

      nnzmin = Math.min(nnzmin, nnz);
      nnzmax = Math.max(nnzmax, nnz);
      nnzavg += nnz / nrepeats;
      nnzstd += nnz*nnz / nrepeats;
    }
    rowstd = Math.ceil(Math.sqrt(rowstd - rowavg*rowavg));
    colstd = Math.ceil(Math.sqrt(colstd - colavg*colavg));
    sizstd = Math.ceil(Math.sqrt(sizstd - sizavg*sizavg));
    nnzstd = Math.ceil(Math.sqrt(nnzstd - nnzavg*nnzavg));
    rowavg = Math.ceil(rowavg);
    colavg = Math.ceil(colavg);
    sizavg = Math.ceil(sizavg);
    nnzavg = Math.ceil(nnzavg);

    int max1 = maximum(countDec(rowmin), countDec(colmin),
                       countDec(sizmin), countDec(nnzmin));
    int max2 = maximum(countDec((int)rowavg), countDec((int)colavg),
                       countDec((int)sizavg), countDec((int)nnzavg));
    int max3 = maximum(countDec((int)rowstd), countDec((int)colstd),
                       countDec((int)sizstd), countDec((int)nnzstd));
    int max4 = maximum(countDec(rowmax), countDec(colmax),
                       countDec(sizmax), countDec(nnzmax));

    System.out.format("  rg        : %" + max1 + "s" +
                      " | %" + max2 + "s" +
                      " +- %-" + max3 + "s" +
                      " | %" + max4 + "s\n",
                      "min", "avg", "std", "max");
    int ccount = 12 + max1 + 3 + max2 + 4 + max3 + 3 + max4;
    System.out.print("  ");
    for (int i = 0; i < ccount; ++i) {
      System.out.print("-");
    }
    System.out.println();

    String fmt =
      " %" + max1 + "d"
      + " | %" + max2 + "d"
      + " +- %" + max3 + "d"
      + " | %" + max4 + "d"
      + "\n";
    System.out.format("  rows      :" + fmt, rowmin, (int)rowavg, (int)rowstd, rowmax);
    System.out.format("  columns   :" + fmt, colmin, (int)colavg, (int)colstd, colmax);
    System.out.format("  entries   :" + fmt, sizmin, (int)sizavg, (int)sizstd, sizmax);
    System.out.format("  non-zeros :" + fmt, nnzmin, (int)nnzavg, (int)nnzstd, nnzmax);
    System.out.println();

    // pd matrices

    sizmin = Integer.MAX_VALUE; sizmax = 0; sizavg = 0.0; sizstd = 0.0;
    final int minsize = (int)Math.ceil(Math.sqrt(MINROWS*MINCOLS));
    final int maxsize = (int)Math.ceil(Math.sqrt(MAXROWS*MAXCOLS));
    for (int r = 0; r < nrepeats; ++r) {
      Matrix A = Benchmark.randPdMatrix(rng,
                                        minsize, maxsize,
                                        MINEIG, MAXEIG, TOL);
      int size = A.rows() * A.cols();
      sizmin = Math.min(sizmin, size);
      sizmax = Math.max(sizmax, size);
      sizavg += size / nrepeats;
      sizstd += size*size / nrepeats;
    }
    sizstd = Math.ceil(Math.sqrt(sizstd - sizavg*sizavg));
    sizavg = Math.ceil(sizavg);

    max1 = maximum(3, countDec(sizmin));
    max2 = maximum(3, countDec((int)sizavg));
    max3 = maximum(3, countDec((int)sizstd));
    max4 = maximum(3, countDec(sizmax));

    System.out.format("  pd        : %" + max1 + "s" +
                      " | %" + max2 + "s" +
                      " +- %-" + max3 + "s" +
                      " | %" + max4 + "s\n",
                      "min", "avg", "std", "max");
    ccount = 12 + max1 + 3 + max2 + 4 + max3 + 3 + max4;
    System.out.print("  ");
    for (int i = 0; i < ccount; ++i) {
      System.out.print("-");
    }
    System.out.println();

    fmt =
      " %" + max1 + "d"
      + " | %" + max2 + "d"
      + " +- %" + max3 + "d"
      + " | %" + max4 + "d"
      + "\n";
    System.out.format("  size      :" + fmt, sizmin, (int)sizavg, (int)sizstd, sizmax);

    System.out.println();
    System.out.flush();
  }

  private void reportBenchmarkResults(BenchmarkData[] data) {
    System.out.println("Running times (ms) and numerical errors:");
    System.out.println();

    int namecmax = 3, mincmax = 3, avgcmax = 3, stdcmax = 3, maxcmax = 3;
    for (BenchmarkData bd : data) {
      if (bd == null) { continue; }
      namecmax = Math.max(namecmax, bd.name.length());
      mincmax = Math.max(mincmax, countDec(bd.time_min));
      avgcmax = Math.max(avgcmax, countDec(bd.time_avg));
      stdcmax = Math.max(stdcmax, countDec(bd.time_std));
      maxcmax = Math.max(maxcmax, countDec(bd.time_max));
    }
    System.out.format("  %" + namecmax + "s" +
                      " : %" + mincmax + "s" +
                      " | %" + avgcmax + "s" +
                      " +- %-" + stdcmax + "s" +
                      " | %" + maxcmax + "s" +
                      " ,   max error \n",
                      " ", "min", "avg", "std", "max");
    int ccount = namecmax + 3 + mincmax + 3 + avgcmax + 4 + stdcmax + 3 + maxcmax + 16;
    System.out.print("  ");
    for (int i = 0; i < ccount; ++i) {
      System.out.print("-");
    }
    System.out.println();

    for (BenchmarkData bd : data) {
      if (bd == null) { continue; }
      System.out.format("  %-" + namecmax + "s" +
                        " : %" + mincmax + "d" +
                        " | %" + avgcmax + "d" +
                        " +- %" + stdcmax + "d" +
                        " | %" + maxcmax + "d",
                        bd.name,
                        bd.time_min,
                        bd.time_avg,
                        bd.time_std,
                        bd.time_max);

      if (Double.isNaN(bd.delta_max)) {
        System.out.print(" , error!");
      }
      else if (bd.delta_max >= 0.0) {
        System.out.format(" , %e", bd.delta_max);
      }
      System.out.println();
    }
  }

  //---------------------------------------------------------------------------

  public static void main(String[] args) {
    BenchmarkRunner runner = new BenchmarkRunner();
    try {
      long seed = SEED;
      if (args.length >= 1 && !args[0].isEmpty()) {
        seed = Long.valueOf(args[0]);
      }

      String[] benchmarks = BENCHMARKS;
      if (args.length >= 2 && !args[1].isEmpty()) {
        benchmarks = args[1].split(",");
      }

      runner.run(seed, benchmarks);
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }
}
