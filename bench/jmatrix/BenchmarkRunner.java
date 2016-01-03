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
  private final String JAVA_CMD = "java";
  private final String CP = "bin:bin" + File.separator + "bench";

  private final String[] BENCHMARKS = {
    "MatMulBenchmark",
    "LUBenchmark",
    "QRBenchmark",
    "ReducedQRBenchmark",
    "QRnoQBenchmark",
    "OrthoNormBenchmark",
    "SVDBenchmark",
    "ReducedSVDBenchmark",
    "CholeskyLBenchmark",
    "CholeskyLDBenchmark",
    "SolveEqnLUBenchmark",
    "SolveEqnLLBenchmark",
    "MatInvBenchmark",
    "MatInvPsdBenchmark"
  };

  private final long SEED = 1927311;
  private final int NWARMUPS = 10;
  private final int NREPEATS = 100;

  private final int MINROWS = 100, MAXROWS = 300;
  private final int MINCOLS =  50, MAXCOLS = 150;

  private final double NNZRATIO = 0.2;
  private final double SCALE = 1000.0;
  private final double TOL = 1e-6;

  private final double MINEIG = 0.0001;
  private final double MAXEIG = 1000.0;

  // Turns matrix statistics report on/off.
  private final boolean MATRIX_STAT = true;

  // Turns benchmark debugging on/off.
  private final boolean DEBUG = false;
  // Anything printed on standard error within the benchmarks
  // will be written into the debug file.
  private final String DEBUG_FILE = "JMATRIX_BENCHMARK_ERROR.txt";

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

  public void run() throws Exception {
    // gather and report matrix statistics
    if (MATRIX_STAT) { reportMatrixStatistics(); }

    // run the benchmarks
    System.out.println("Running benchmarks...");
    System.out.println();
    BenchmarkData[] data = new BenchmarkData[BENCHMARKS.length];
    for (int b = 0; b < BENCHMARKS.length; ++b) {
      String classname = BENCHMARKS[b];
      System.out.println("  " + classname);

      ProcessBuilder pb = new ProcessBuilder(JAVA_CMD,
                                             "-d64", "-Xms512m", "-Xmx512m",
                                             "-cp", CP,
                                             "jmatrix." + classname,
                                             "" + SEED,
                                             "" + NWARMUPS,
                                             "" + NREPEATS,
                                             "" + MINROWS, "" + MAXROWS,
                                             "" + MINCOLS, "" + MAXCOLS,
                                             "" + NNZRATIO,
                                             "" + SCALE,
                                             "" + TOL,
                                             "" + MINEIG,
                                             "" + MAXEIG);
      if (DEBUG) { pb.redirectError(new File(DEBUG_FILE)); }

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

  private void reportMatrixStatistics() {
    Random rng = new Random(SEED);

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
                        " ,    max error\n",
                        " ", "min", "avg", "std", "max");
    int ccount = namecmax + 3 + mincmax + 3 + avgcmax + 4 + stdcmax + 3 + maxcmax + 15;
    System.out.print("  ");
    for (int i = 0; i < ccount; ++i) {
      System.out.print("-");
    }
    System.out.println();

    for (BenchmarkData bd : data) {
      System.out.format("  %-" + namecmax + "s" +
                        " : %" + mincmax + "d" +
                        " | %" + avgcmax + "d" +
                        " +- %" + stdcmax + "d" +
                        " | %" + maxcmax + "d" +
                        " , %e" +
                        "\n",
                        bd.name,
                        bd.time_min,
                        bd.time_avg,
                        bd.time_std,
                        bd.time_max,
                        bd.delta_max);
    }
  }

  //---------------------------------------------------------------------------

  public static void main(String[] args) {
    BenchmarkRunner runner = new BenchmarkRunner();
    try {
      runner.run();
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }
}
