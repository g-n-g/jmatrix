package jmatrix;

/** Represents benchmark data. */
public final class BenchmarkData
{
  public final String name;
  public final int time_min, time_avg, time_std, time_max;
  public final double delta_max;

  public BenchmarkData(String name,
                       int time_min, int time_max,
                       int time_avg, int time_std,
                       double delta_max) {
    this.name = name;
    this.time_min = time_min;
    this.time_max = time_max;
    this.time_avg = time_avg;
    this.time_std = time_std;
    this.delta_max = delta_max;
  }

  //---------------------------------------------------------------------------

  private static final String DELIM = "!";

  @Override
  public String toString() {
    return name + DELIM + time_min + DELIM + time_max
      + DELIM + time_avg + DELIM + time_std + DELIM + delta_max;
  }

  public static BenchmarkData valueOf(String s) throws BenchmarkException {
    String[] tokens = s.split(DELIM);
    if (tokens.length != 6) {
      throw new BenchmarkException("Invalid BenchmarkData string: " + s + "!");
    }
    return new BenchmarkData(tokens[0],
                             Integer.valueOf(tokens[1]),
                             Integer.valueOf(tokens[2]),
                             Integer.valueOf(tokens[3]),
                             Integer.valueOf(tokens[4]),
                             Double.valueOf(tokens[5]));
  }
}
