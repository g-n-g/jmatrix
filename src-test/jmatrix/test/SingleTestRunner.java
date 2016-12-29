package jmatrix.test;

import org.junit.runner.JUnitCore;
import org.junit.runner.Request;
import org.junit.runner.Result;
import org.junit.runner.notification.Failure;

public class SingleTestRunner {

  public static void main(String... args) throws ClassNotFoundException {
    if (args.length == 0 || args[0].equals("jmatrix.")) {
      System.out.println("Test case is not specified, aborting...");
      System.exit(1);
    }
    String[] classAndMethod = args[0].split("#");
    Request request = Request.method(Class.forName(classAndMethod[0]),
                                     classAndMethod[1]);

    System.out.println("Test case: " + args[0]);
    Result result = new JUnitCore().run(request);
    for (Failure fail : result.getFailures()) {
      System.out.println("Fail reason: " + fail.getMessage());
      System.out.println(fail.getTrace());
    }
    System.out.println("TEST " + (result.getFailureCount() == 1 ?
                                  "FAILED!" : "PASSED."));
    System.exit(0);
  }
}
