/*
  Interface to call bp methods from figaro

  Kathryn Rodgers 2015
 */

class bp_interface{

    /**
     Does inference using BP on the model described in the UAI file
     @params uaiFileName the name of the uai file that describes the model
     */
    //TODO: be able to choose the task
    public native boolean doBP(String uaiFileName);

    /**
     Returns the solution
     */
    public native String getSolution();
    
    static{
	System.out.println("loading library");
        System.loadLibrary("BP");
    }

}
