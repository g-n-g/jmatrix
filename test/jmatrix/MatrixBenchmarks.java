package jmatrix;

import java.lang.System;
import java.util.Random;
import java.util.List;
import java.util.LinkedList;

import org.junit.Test;
import org.junit.AfterClass;
import org.junit.FixMethodOrder;
import org.junit.runners.MethodSorters;
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
@FixMethodOrder(MethodSorters.NAME_ASCENDING)
public class MatrixBenchmarks {

  public static final int REPEAT_COUNT = 100;
  public static final int SEED = 19273;
  public static final int MINSIZE = 50;
  public static final int MAXSIZE = 150;
  public static final double SCALE = 100.0;
  public static final double NNZ_RATIO = 0.2;

  //---------------------------------------------------------------------------

  @Test public void LU() { bench("LU", new BenchLU(), false); }
  @Test public void QR() { bench("QR", new BenchQR(), false); }
  @Test public void reducedQR() { bench("QR (reduced)", new BenchReducedQR(), false); }
  @Test public void noQQR() { bench("QR (noQ)", new BenchNoQQR(), false); }
  @Test public void CholeskyL() { bench("Cholesky L", new BenchCholeskyL(), true); }
  @Test public void CholeskyLD() { bench("Cholesky LD", new BenchCholeskyLD(), true); }
  @Test public void reducedSVD() { bench("SVD (reduced)", new BenchReducedSVD(), false); }

  //---------------------------------------------------------------------------

  private interface Bench
  {
    void compute(Matrix A);
    void check(Matrix A);
  }

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

  private class BenchReducedQR implements Bench
  {
    public void compute(Matrix A) {
      Matrix[] QR = A.reducedQR();
      Q = QR[0];
      R = QR[1];
    }

    public void check(Matrix A) {
      assertMatrixEquals(A, Q.mul(R), TOL);
      if (A.rows() > A.cols()) {
        assertMatrixOrthoCols(Q, TOL);
      }
      else {
        assertMatrixOrtho(Q, TOL);
      }
      assertMatrixUT(R, TOL);
    }

    private Matrix Q, R;
  }

  private class BenchNoQQR implements Bench
  {
    public void compute(Matrix A) {
      R = Matrix.create(A.rows(), A.cols());
      final int t = Math.min(A.rows()-1, A.cols());
      A.QR(null, R, Matrix.create(t,1));
    }

    public void check(Matrix A) {
      assertMatrixUT(R, TOL);
      Matrix[] QR = A.QR();
      Matrix Q = QR[0];
      assertMatrixEquals(A, Q.mul(R), TOL);
      assertMatrixOrtho(Q, TOL);
    }

    private Matrix R;
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

  //---------------------------------------------------------------------------

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

  private static Matrix genA(Random rng) {
    int Arows = MINSIZE+rng.nextInt(MAXSIZE-MINSIZE);
    int Acols = MINSIZE+rng.nextInt(MAXSIZE-MINSIZE);
    // System.out.println("rows = " + rows + ", cols = " + cols);
    Matrix A = Matrix.randN(Arows, Acols, rng).mul(rng.nextDouble()*SCALE);
    if (NNZ_RATIO > 0) {
      int nnz = (int)Math.round(NNZ_RATIO*(Arows*Acols));
      for (int k = 0; k < nnz; ++k) {
        int i = rng.nextInt(Arows);
        int j = rng.nextInt(Acols);
        A.set(i, j, (2.0*rng.nextDouble()-1.0)*TOL);
      }
    }
    return A;
  }

  private void bench(String name, Bench test, boolean isPsd) {
    Random rng = new Random(SEED);

    final double nrepeats = REPEAT_COUNT;
    long times_min = Long.MAX_VALUE, times_max = 0;
    double times_avg = 0, times_std = 0;

    for (int repeat = 0; repeat < REPEAT_COUNT; ++repeat) {
      System.gc();

      Matrix A = genA(rng);
      if (isPsd) { A = A.mul(A.T()).add(Matrix.eye(A.rows()).mul(2*TOL)); }

      long ts = System.currentTimeMillis();
      test.compute(A.copy());
      ts = System.currentTimeMillis() - ts;
      test.check(A);

      if (ts < times_min) { times_min = ts; }
      if (ts > times_max) { times_max = ts; }
      times_avg += ts / nrepeats;
      times_std += ts*ts / nrepeats;
    }
    times_std = Math.ceil(Math.sqrt(times_std - times_avg*times_avg));

    results.add(new Result(name, times_min, times_max,
                           (long)times_avg, (long)times_std));
  }

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

  private static double cmputStd(double avgsum, double stdsum, int count) {
    double avg = avgsum / count;
    return Math.sqrt((stdsum - avg*avg) / (count-1));
  }

  @AfterClass
  public static void logResults() {
    System.out.println("\n");
    {
      System.out.println("Matrix statistics:");
      System.out.println();

      final double nrepeats = REPEAT_COUNT;
      int rowmin = Integer.MAX_VALUE, rowmax = 0; double rowavg = 0, rowstd = 0;
      int colmin = Integer.MAX_VALUE, colmax = 0; double colavg = 0, colstd = 0;
      int sizmin = Integer.MAX_VALUE, sizmax = 0; double sizavg = 0, sizstd = 0;
      int nnzmin = Integer.MAX_VALUE, nnzmax = 0; double nnzavg = 0, nnzstd = 0;
      Random rng = new Random(SEED);
      for (int repeat = 0; repeat < REPEAT_COUNT; ++repeat) {
        Matrix A = genA(rng);

        int rows = A.rows(), cols = A.cols();
        int size = rows * cols;
        int nnz = 0;
        for (int i = 0; i < rows; ++i) {
          for (int j = 0; j < cols; ++j) {
            if (Math.abs(A.get(i,j)) < TOL) { ++nnz; }
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
    }
    System.out.println();
    {
      System.out.println("Running times (ms):");
      System.out.println();

      int namecmax = 0, mincmax = 0, avgcmax = 0, stdcmax = 0, maxcmax = 0;
      for (Result result : results) {
        namecmax = Math.max(namecmax, result.name.length());
        mincmax = Math.max(mincmax, countDec(result.min));
        avgcmax = Math.max(avgcmax, countDec(result.avg));
        stdcmax = Math.max(stdcmax, countDec(result.std));
        maxcmax = Math.max(maxcmax, countDec(result.max));
      }

      for (Result result : results) {
        System.out.format("  %-" + namecmax + "s :"
                          + " %" + mincmax + "d" +
                          " | %" + avgcmax + "d" +
                          " +- %" + stdcmax + "d" +
                          " | %" + maxcmax + "d\n",
                          result.name,
                          result.min,
                          result.avg,
                          result.std,
                          result.max);
      }
    }
  }
}
