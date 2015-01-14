/**
 * Created by krodgers on 1/14/15.
 */

object Main {

  def main(args:Array[String]): Unit ={

    // Create Figaro Model

    // Write Model to UAI File

    // Run Inference on UAI Model

    //bp_interface bpi = new bp_interface();
    //bpi.doBP("string!");
    //bpi.getSolution();


    val libMap = System.mapLibraryName("BP");
    println(libMap);
    //System.loadLibrary(libMap);

    val bpi = new bp_interface()
    bpi.doBP("../data/test.uai")
    val res = bpi.getSolution
    println("solution can be found in " + res)

    // Get the solution back

    // Print out solution


  }
}
