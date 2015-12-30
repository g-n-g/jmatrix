package jmatrix;

import java.lang.System;
import java.util.Random;
import java.util.List;
import java.util.LinkedList;

import org.junit.Test;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import static jmatrix.MatrixAssert.assertMatrixEquals;
import static jmatrix.MatrixAssert.assertMatrixUnitLT;
import static jmatrix.MatrixAssert.assertMatrixLT;
import static jmatrix.MatrixAssert.assertMatrixUT;
import static jmatrix.MatrixAssert.assertMatrixDiag;
import static jmatrix.MatrixAssert.assertMatrixOrtho;
import static jmatrix.MatrixAssert.assertMatrixOrthoCols;

import static jmatrix.Matrix.NR;
import static jmatrix.Matrix.TOL;
import static jmatrix.BasicUnaryOperation.*;
import static jmatrix.BasicBinaryOperation.*;

/**
 * Benchmarks for the matrix class.
 */
public class MatrixBenchmarks {

  public static final int REPEAT_COUNT = 100;
  public static final int SEED = 19273;
  public static final int MINSIZE = 50;
  public static final int MAXSIZE = 100;

  //---------------------------------------------------------------------------

  private interface Bench
  {
    void compute(Matrix A);
    void check(Matrix A);
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

  private void bench(String name, Bench test, boolean isPsd) {
    Random rng = new Random(SEED);

    long[] times = new long[REPEAT_COUNT];
    for (int repeat = 0; repeat < REPEAT_COUNT; ++repeat) {
      int rows = MINSIZE+rng.nextInt(MAXSIZE-MINSIZE);
      int cols = MINSIZE+rng.nextInt(MAXSIZE-MINSIZE);
      // System.out.println("rows = " + rows + ", cols = " + cols);
      Matrix A = Matrix.randN(rows, cols, rng);
      if (isPsd) { A = A.mul(A.T()).add(Matrix.eye(rows).mul(2*TOL)); }

      long ts = System.currentTimeMillis();
      test.compute(A.copy());
      times[repeat] = System.currentTimeMillis() - ts;
      test.check(A);
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

  @Test public void LU() { bench("LU", new BenchLU(), false); }
  @Test public void QR() { bench("QR", new BenchQR(), false); }
  @Test public void CholeskyL() { bench("Cholesky L",
                                        new BenchCholeskyL(),
                                        true); }
  @Test public void CholeskyLD() { bench("Cholesky LD",
                                         new BenchCholeskyLD(),
                                         true); }
  @Test public void reducedSVD() { bench("SVD (reduced)",
                                         new BenchReducedSVD(),
                                         false); }

  //---------------------------------------------------------------------------

  private class BenchLU implements Bench
  {
    public void compute(Matrix A) {
      Matrix[] LUP = A.LU();
      L = LUP[0];
      U = LUP[1];
      P = LUP[2];
    }

    public void check(Matrix A) {
      assertMatrixEquals(A, P.T().mul(L).mul(U), TOL);
      assertMatrixUnitLT(L, TOL);
      assertMatrixUT(U, TOL);
    }

    private Matrix L, U, P;
  }

  private class BenchQR implements Bench
  {
    public void compute(Matrix A) {
      Matrix[] QR = A.QR();
      Q = QR[0];
      R = QR[1];
    }

    public void check(Matrix A) {
      assertMatrixEquals(A, Q.mul(R), TOL);
      assertMatrixOrtho(Q, TOL);
      assertMatrixUT(R, TOL);
    }

    private Matrix Q, R;
  }

  private class BenchCholeskyL implements Bench
  {
    public void compute(Matrix A) {
      L = A.choleskyL();
    }

    public void check(Matrix A) {
      assertEquals(L.rows(), L.cols());
      assertMatrixEquals(A, L.mul(L.T()), TOL);
      assertMatrixLT(L, TOL);
    }

    private Matrix L;
  }

  private class BenchCholeskyLD implements Bench
  {
    public void compute(Matrix A) {
      Matrix[] LD = A.choleskyLD();
      L = LD[0];
      D = LD[1];
    }

    public void check(Matrix A) {
      assertEquals(L.rows(), L.cols());
      assertMatrixEquals(A, L.mul(D).mul(L.T()), TOL);
      assertMatrixLT(L, TOL);
      assertMatrixDiag(D, TOL);
    }

    private Matrix L, D;
  }

  private class BenchReducedSVD implements Bench
  {
    public void compute(Matrix A) {
      Matrix[] USV = A.reducedSVD();
      U = USV[0];
      S = USV[1];
      V = USV[2];
    }

    public void check(Matrix A) {
      assertMatrixEquals(A, U.ewb(MUL, S.T()).mul(V.T()), TOL);
      assertMatrixOrthoCols(U, TOL);
      assertMatrixOrthoCols(V, TOL);
    }

    private Matrix U, S, V;
  }
}
