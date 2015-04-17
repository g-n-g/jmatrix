package jvecmat;

import junit.framework.TestCase;

/**
 * Enable the implementation assertions.
 */
public class AssertionBaseTest extends TestCase {

  public AssertionBaseTest(String name) {
    super(name);
  }

  static {
    ClassLoader.getSystemClassLoader().setDefaultAssertionStatus(true);
  }

  //----------------------------------------------------------------------------

  protected void log(Vector v) {
    System.out.println(v);
  }

  protected void log(Matrix m) {
    System.out.println(m);
  }

  protected void log(Permutation p) {
    System.out.println(p);
  }
}
