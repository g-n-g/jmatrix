package ai.gandg.jmatrix_bench;

import java.util.Map;
import java.util.HashMap;
import ai.gandg.jmatrix_bench.BenchmarkData.MatrixStat;
import ai.gandg.jmatrix_bench.Benchmark.BenchmarkType;


/** Collection of benchmark matrix statistics data. */
public class BenchmarkReport
{
  public void printResults(BenchmarkData[] data) {
    System.out.println("Input types, running times (ms), and numerical errors:");
    System.out.println();

    int namecmax = 3, inputcmax = 5;
    int mincmax = 3, avgcmax = 3, stdcmax = 3, maxcmax = 3;
    for (BenchmarkData bd : data) {
      if (bd == null) { continue; }
      namecmax = Math.max(namecmax, bd.name.length());
      inputcmax = Math.max(inputcmax, ("" + bd.type).length());
      mincmax = Math.max(mincmax, countDec(bd.time_min));
      avgcmax = Math.max(avgcmax, countDec(bd.time_avg));
      stdcmax = Math.max(stdcmax, countDec(bd.time_std));
      maxcmax = Math.max(maxcmax, countDec(bd.time_max));
    }
    System.out.format("  %" + namecmax + "s" +
                      " : %" + inputcmax + "s" +
                      " | %" + mincmax + "s" +
                      " | %" + avgcmax + "s" +
                      " +- %-" + stdcmax + "s" +
                      " | %" + maxcmax + "s" +
                      " |   max error \n",
                      " ", "input", "min", "avg", "std", "max");
    int ccount = namecmax + 3 + inputcmax + 3
      + mincmax + 3 + avgcmax + 4 + stdcmax + 3 + maxcmax + 16;
    System.out.print("  ");
    for (int i = 0; i < ccount; ++i) {
      System.out.print("-");
    }
    System.out.println();

    for (BenchmarkData entry : data) {
      if (entry == null) { continue; }
      System.out.format("  %-" + namecmax + "s" +
                        " : %" + inputcmax + "s" +
                        " | %" + mincmax + "d" +
                        " | %" + avgcmax + "d" +
                        " +- %" + stdcmax + "d" +
                        " | %" + maxcmax + "d",
                        entry.name,
                        "" + entry.type,
                        entry.time_min,
                        entry.time_avg,
                        entry.time_std,
                        entry.time_max);

      if (Double.isNaN(entry.delta_max)) {
        System.out.print(" | error!");
      }
      else if (entry.delta_max >= 0.0) {
        System.out.format(" | %e", entry.delta_max);
      }
      System.out.println();
    }
  }

  private static class ReportStat
  {
    static class Entry
    {
      int min, max;
      double avg, std;

      Entry(MatrixStat.Entry entry) {
        min = (int)entry.min;
        max = (int)entry.max;
        avg = entry.avg;
        std = Math.ceil(Math.sqrt(entry.var - avg*avg));
        avg = Math.ceil(avg);
      }
    }
    private BenchmarkType type;
    private Entry rows, cols, size, nnzs;

    ReportStat(BenchmarkData entry) {
      if (entry == null || entry.statA == null) { return; }

      this.type = entry.type;
      MatrixStat stat = entry.statA;
      this.rows = new Entry(stat.rows);
      this.cols = new Entry(stat.cols);
      this.size = new Entry(stat.size);
      this.nnzs = new Entry(stat.nnzs);
    }

    void print() {
      int max1 = maximum(countDec(rows.min), countDec(cols.min),
                         countDec(size.min), countDec(nnzs.min), 3);
      int max2 = maximum(countDec((int)rows.avg), countDec((int)cols.avg),
                         countDec((int)size.avg), countDec((int)nnzs.avg), 3);
      int max3 = maximum(countDec((int)rows.std), countDec((int)cols.std),
                         countDec((int)size.std), countDec((int)nnzs.std), 3);
      int max4 = maximum(countDec(rows.max), countDec(cols.max),
                         countDec(size.max), countDec(nnzs.max), 3);

      String fmt = "  %-8s  : %" + max1 + "s | %"
        + max2 + "s +- %-" + max3 + "s | %" + max4 + "s\n";
      System.out.format(fmt, type, "min", "avg", "std", "max");

      int ccount = 12 + max1 + 3 + max2 + 4 + max3 + 3 + max4;
      System.out.print("  ");
      for (int i = 0; i < ccount; ++i) { System.out.print("-"); }
      System.out.println();

      fmt = " %" + max1 + "d | %" + max2
        + "d +- %" + max3 + "d | %" + max4 + "d\n";
      System.out.format("  rows      :" + fmt,
                        rows.min, (int)rows.avg, (int)rows.std, rows.max);
      System.out.format("  columns   :" + fmt,
                        cols.min, (int)cols.avg, (int)cols.std, cols.max);
      System.out.format("  sizes     :" + fmt,
                        size.min, (int)size.avg, (int)size.std, size.max);
      System.out.format("  non-zeros :" + fmt,
                        nnzs.min, (int)nnzs.avg, (int)nnzs.std, nnzs.max);
      System.out.println();
    }
  }

  public void printMatrixStats(BenchmarkData[] data) {
    Map<BenchmarkType, ReportStat> stats =
      new HashMap<BenchmarkType, ReportStat>();

    System.out.format("Matrix statistics:\n\n");
    for (BenchmarkData entry : data) {
      if (entry == null || entry.statA == null) { continue; }

      BenchmarkType type = entry.type;
      ReportStat stat = new ReportStat(entry);
      if (stats.containsKey(type)) {
        assert (stat.equals(stats.get(type)));
      }
      else {
        stats.put(type, stat);
        stat.print();
      }
    }
    System.out.flush();
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
}
