/*
  Interface to call bp methods from figaro

  Kathryn Rodgers 2015
 */

class bp_interface{
    private long interfacePtr; // pointer to BpInterface object

    public native long initBpInterface();
    public native void destroyBpInterface();

    public bp_interface(){
	interfacePtr = this.initBpInterface(); // create the object
    }
    public void destroy(){
	destroyBpInterface();// delete the object
	interfacePtr = 0L;
    }
    protected void finalize() throws Throwable{
	destroyBpInterface();
	interfacePtr = 0L;
    }

    /**
     Does inference using BP on the model described in the UAI file
     @params uaiFileName the name of the uai file that describes the model
     */
    //TODO: be able to choose the task
    public native boolean doBP( String uaiFileName);

    /**
     Returns the solution
     */
    public native String getSolution();
    
    static{
	System.out.println("loading library");
	System.loadLibrary("BP");
    }

}
