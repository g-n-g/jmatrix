package jmatrix.bench;

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
    "MatInv",
    "MatInvPsd",
    "SolvePDEqnLL",
    "SolvePDEqnLU",
    "SolveLSFCREqnLL",
    "SolveLSFCREqnQRnoQ",
    "SolveLSFCREqnReducedQR"
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

  // Turns benchmark debugging on/off.
  private static final boolean DEBUG = false;
  // Anything printed on standard error within the benchmarks
  // will be written into the debug files named with this prefix.
  private static final String DEBUG_FILE = "JMATRIX_BENCHMARK_ERROR_";

  // Java command and classpath.
  private static final String JAVA_CMD = "java";
  private static final String CP = "bin:bin" + File.separator + "bench";

  //---------------------------------------------------------------------------

  public void run(long seed, String[] benchmarks, boolean isDebug)
    throws Exception {

    // determine random seed
    if (0 == seed) {
      seed = new Random().nextLong();
    }
    System.out.println("Random seed: " + seed);
    System.out.println();

    // run the benchmarks
    double ts = System.nanoTime();
    System.out.println("Running benchmarks...");
    System.out.format("  %25s : %10s\n", "", "elapsed(s)");
    System.out.flush();
    BenchmarkData[] data = new BenchmarkData[BENCHMARKS.length];
    for (int b = 0; b < benchmarks.length; ++b) {
      data[b] = null;
      
      String benchmark = benchmarks[b];
      if (!existsBenchmark(benchmark)) {
        System.out.println("  " + benchmark + " does not exist!");
        continue;
      }
      System.out.format("  %25s : %10.0f\n", benchmark, (System.nanoTime()-ts) / 1e9);
      System.out.flush();
      
      String classname = benchmark + "Benchmark";
      ProcessBuilder pb = new ProcessBuilder(JAVA_CMD,
                                             "-d64", "-Xms512m", "-Xmx512m",
                                             "-cp", CP,
                                             "jmatrix.bench.Benchmark",
                                             "jmatrix.bench." + classname,
                                             "" + isDebug,
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
      if (isDebug) { pb.redirectError(new File(DEBUG_FILE + b)); }

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

    // reporting
    BenchmarkReport report = new BenchmarkReport();
    report.printMatrixStats(data);
    report.printResults(data);
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

      boolean debug = DEBUG;
      if (args.length >= 3 && !args[2].isEmpty()) {
        debug = Boolean.valueOf(args[2]);
      }

      runner.run(seed, benchmarks, debug);
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }
}
