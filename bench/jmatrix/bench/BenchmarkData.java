package jmatrix.bench;

import jmatrix.Matrix;
import jmatrix.bench.Benchmark.BenchmarkType;

/** Represents benchmark data. */
public final class BenchmarkData
{
  public final String name;
  BenchmarkType type;
  public final int time_min, time_avg, time_std, time_max;
  public final double delta_max;
  public final MatrixStat statA;

  public BenchmarkData(String name) {
    this.name = name;
    type = null;
    time_min = time_avg = time_std = time_max = 0;
    delta_max = Double.NaN;
    statA = null;
  }

  public BenchmarkData(String name, BenchmarkType type,
                       int time_min, int time_max,
                       int time_avg, int time_std,
                       double delta_max, MatrixStat statA) {
    this.name = name;
    this.type = type;
    this.time_min = time_min;
    this.time_max = time_max;
    this.time_avg = time_avg;
    this.time_std = time_std;
    this.delta_max = delta_max;
    this.statA = statA;
  }

  //---------------------------------------------------------------------------

  @Override
  public String toString() {
    return name + DELIM + type + DELIM + time_min + DELIM + time_max
      + DELIM + time_avg + DELIM + time_std + DELIM + delta_max + DELIM + statA;
  }

  public static BenchmarkData valueOf(String s) throws BenchmarkException {
    String[] tokens = s.split(DELIM);
    if (tokens.length != 8) {
      throw new BenchmarkException("Invalid BenchmarkData string: " + s);
    }
    return new BenchmarkData(tokens[0],
                             BenchmarkType.valueOf(tokens[1]),
                             Integer.valueOf(tokens[2]),
                             Integer.valueOf(tokens[3]),
                             Integer.valueOf(tokens[4]),
                             Integer.valueOf(tokens[5]),
                             Double.valueOf(tokens[6]),
                             MatrixStat.valueOf(tokens[7]));
  }

  //---------------------------------------------------------------------------

  public static class MatrixStat
  {
    int count;
    Entry rows, cols, size, nnzs;

    MatrixStat(int count, Entry rows, Entry cols, Entry size, Entry nnzs) {
      this.count = count;
      this.rows = rows;
      this.cols = cols;
      this.size = size;
      this.nnzs = nnzs;
    }

    MatrixStat() {
      count = 0;
      rows = new Entry();
      cols = new Entry();
      size = new Entry();
      nnzs = new Entry();
    }

    public void add(Matrix M) {
      ++count;
      double rows = M.rows();
      this.rows.add(rows, count);
      double cols = M.cols();
      this.cols.add(cols, count);
      double size = rows * cols;
      this.size.add(size, count);
      double nnzs = M.nnz();
      this.nnzs.add(nnzs, count);
    }

    @Override
    public String toString() {
      return "" + count
        + DELIM_MATRIX_STAT + rows + DELIM_MATRIX_STAT + cols
        + DELIM_MATRIX_STAT + size + DELIM_MATRIX_STAT + nnzs;
    }

    public static MatrixStat valueOf(String s) throws BenchmarkException {
      String[] tokens = s.split(DELIM_MATRIX_STAT);
      if (tokens.length != 5) {
        throw new BenchmarkException("Invalid MatrixStat string: " + s);
      }
      return new MatrixStat(Integer.valueOf(tokens[0]),
                            Entry.valueOf(tokens[1]),
                            Entry.valueOf(tokens[2]),
                            Entry.valueOf(tokens[3]),
                            Entry.valueOf(tokens[4]));
    }

    @Override
    public boolean equals(Object aThat) {
      if (this == aThat) { return true; }
      if (!(aThat instanceof MatrixStat)) { return false; }
      MatrixStat that = (MatrixStat)aThat;
      return
        this.count == that.count &&
        this.rows.equals(that.rows) &&
        this.cols.equals(that.cols) &&
        this.size.equals(that.size) &&
        this.nnzs.equals(that.nnzs);
    }


    public static class Entry
    {
      // minimum, maximum, average, "variance"
      double min, max, avg, var;

      Entry(double min, double max, double avg, double var) {
        this.min = min;
        this.max = max;
        this.avg = avg;
        this.var = var;
      }

      Entry() {
        min = Double.MAX_VALUE;
        max = Double.MIN_VALUE;
        avg = var = 0.0;
      }

      void add(double value, int count) {
        min = Math.min(min, value);
        max = Math.max(max, value);
        double ratio = (count-1) / (double)count;
        avg = avg*ratio + value/count;
        var = var*ratio + (value*value)/count;
      }

      @Override
      public String toString() {
        return "" + min + DELIM_MATRIX_STAT_ENTRY + max +
          DELIM_MATRIX_STAT_ENTRY + avg + DELIM_MATRIX_STAT_ENTRY + var;
      }

      public static Entry valueOf(String s) throws BenchmarkException {
        String[] tokens = s.split(DELIM_MATRIX_STAT_ENTRY);
        if (tokens.length != 4) {
          throw new BenchmarkException("Invalid Entry string: " + s);
        }
        return new Entry(Double.valueOf(tokens[0]),
                         Double.valueOf(tokens[1]),
                         Double.valueOf(tokens[2]),
                         Double.valueOf(tokens[3]));
      }

      @Override
      public boolean equals(Object aThat) {
        if (this == aThat) { return true; }
        if (!(aThat instanceof Entry)) { return false; }
        Entry that = (Entry)aThat;
        return
          this.min == that.min &&
          this.max == that.max &&
          this.avg == that.avg &&
          this.var == that.var;
      }
    }
  }

  //---------------------------------------------------------------------------

  private static final String DELIM = "!!!";
  private static final String DELIM_MATRIX_STAT = "!!";
  private static final String DELIM_MATRIX_STAT_ENTRY = "!";
}
