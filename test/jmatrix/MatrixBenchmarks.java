package jmatrix;

import java.lang.System;
import java.util.Random;
import java.util.List;
import java.util.LinkedList;

import org.junit.Test;
import org.junit.AfterClass;
import static org.junit.Assert.assertTrue;

import static jmatrix.Matrix.NR;
import static jmatrix.Matrix.TOL;
import static jmatrix.BasicUnaryOperation.*;
import static jmatrix.BasicBinaryOperation.*;

/**
 * Benchmarks for the matrix class.
 */
public class MatrixBenchmarks {

  public static final double PREC = 1e-8;
  public static final int REPEAT_COUNT = 100;
  public static final int SEED = 19273;
  public static final int MINSIZE = 50;
  public static final int MAXSIZE = 100;

  //---------------------------------------------------------------------------

  private interface Bench
  {
    void run(Matrix A);
  }

  private class Result
  {
    Result(String name, long min, long max, long avg, long std) {
      this.name = name;
      this.min = min;
      this.max = max;
      this.avg = avg;
      this.std = std;
    }

    String name;
    long min, max, avg, std;
  }
  private static List<Result> results = new LinkedList<Result>();

  private void bench(String name, Bench test) {
    Random rng = new Random(SEED);

    long[] times = new long[REPEAT_COUNT];
    for (int repeat = 0; repeat < REPEAT_COUNT; ++repeat) {
      int rows = MINSIZE+rng.nextInt(MAXSIZE-MINSIZE);
      int cols = MINSIZE+rng.nextInt(MAXSIZE-MINSIZE);
      // System.out.println("rows = " + rows + ", cols = " + cols);
      Matrix A = Matrix.randN(rows, cols, rng);

      long ts = System.currentTimeMillis();
      test.run(A);
      times[repeat] = System.currentTimeMillis() - ts;
    }

    long min = Long.MAX_VALUE, max = 0, avg = 0;
    for (int repeat = 0; repeat < REPEAT_COUNT; ++repeat) {
      long time = times[repeat];
      if (min > time) { min = time; }
      if (max < time) { max = time; }
      avg += time;
    }
    avg /= REPEAT_COUNT;

    long std = 0;
    for (int repeat = 0; repeat < REPEAT_COUNT; ++repeat) {
      long diff = times[repeat] - avg;
      std += diff*diff;
    }
    std /= REPEAT_COUNT-1;
    std = (long)Math.sqrt((double)std);

    results.add(new Result(name, min, max, avg, std));
  }

  @AfterClass
  public static void logResults() {
    System.out.println("\n");
    System.out.println("Times (ms):");
    System.out.println();
    for (Result result : results) {
      System.out.format("  %-15s : %5d | %5d +- %5d | %5d\n", result.name,
                        result.min, result.avg, result.std, result.max);
    }
  }

  //---------------------------------------------------------------------------

  @Test public void QR() { bench("QR", new BenchQR()); }
  @Test public void reducedSVD() { bench("SVD (reduced)", new BenchReducedSVD()); }

  //---------------------------------------------------------------------------

  private class BenchQR implements Bench
  {
    public void run(Matrix A) {
      Matrix[] QR = A.QR();
      Matrix Q = QR[0]; Matrix R = QR[1];
      assertTrue(PREC > A.sub(Q.mul(R)).norm1());
      assertTrue(PREC > Matrix.eye(A.rows()).sub(Q.T().mul(Q)).norm1());
      assertTrue(PREC > Matrix.eye(A.rows()).sub(Q.mul(Q.T())).norm1());
    }
  }

  private class BenchReducedSVD implements Bench
  {
    public void run(Matrix A) {
      Matrix[] USV = A.reducedSVD();
      Matrix U = USV[0], S = USV[1], V = USV[2];
      assertTrue(PREC > U.ewb(MUL, S.T()).mul(V.T()).sub(A).normF());
      int rank = S.rows();
      assertTrue(PREC > U.T().mul(U).sub(Matrix.eye(rank)).normF());
      assertTrue(PREC > V.T().mul(V).sub(Matrix.eye(rank)).normF());

      Matrix A1 = A.div(TOL);
      USV = A1.reducedSVD();
      Matrix U1 = USV[0], S1 = USV[1], V1 = USV[2];
      assertTrue(PREC > U1.sub(U).normF());
      assertTrue(PREC > V1.sub(V).normF());
      assertTrue(PREC > S1.mul(TOL).sub(S).normF());

      Matrix A2 = A.mul(TOL);
      USV = A2.reducedSVD();
      Matrix U2 = USV[0], S2 = USV[1], V2 = USV[2];
      assertTrue(PREC > U2.sub(U).normF());
      assertTrue(PREC > V2.sub(V).normF());
      assertTrue(PREC > S2.div(TOL).sub(S).normF());
    }
  }
}
