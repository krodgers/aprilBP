#include "iostream"
#include "appril_java_interface.h"



/*
 * Class:     appril_java_interface
 * Method:    MapErrorCode2String
 * Signature: (I)[C
 */
JNIEXPORT jcharArray JNICALL Java_appril_1java_1interface_MapErrorCode2String(JNIEnv *, jobject, jint){

  printf( "map2string\n");

  return NULL;

}

/*
 * Class:     appril_java_interface
 * Method:    Initialize
 * Signature: (Z)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_Initialize__Z(JNIEnv *, jobject, jboolean)
[
  printf( "map2string\n");

  return 0;
  }

/*
 * Class:     appril_java_interface
 * Method:    Initialize
 * Signature: ()I
 */
  JNIEXPORT jint JNICALL Java_appril_1java_1interface_Initialize__(JNIEnv *, jobject){
  printf( "map2string\n");
  return 0;
}

/*
 * Class:     appril_java_interface
 * Method:    Terminate
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_Terminate(JNIEnv *, jobject){
  printf( "map2string\n");

  return 0;
}

/*
 * Class:     appril_java_interface
 * Method:    GetCurrentState
 * Signature: ([CI)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_GetCurrentState(JNIEnv *, jobject, jcharArray, jint){

  printf( "map2string\n");
  return 0;
}

/*
 * Class:     appril_java_interface
 * Method:    DefineProblem
 * Signature: ([C[C[C)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_DefineProblem(JNIEnv *, jobject, jcharArray, jcharArray, jcharArray){

  printf( "map2string\n");
  return 0;
  }


/*
 * Class:     appril_java_interface
 * Method:    GetProblemComplexity
 * Signature: ([CI)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_GetProblemComplexity(JNIEnv *, jobject, jcharArray, jint)]
  printf( "map2string\n");
return 0;
}

/*
 * Class:     appril_java_interface
 * Method:    GetSolutionComplexity
 * Signature: ([CI)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_GetSolutionComplexity(JNIEnv *, jobject, jcharArray, jint){
  printf( "map2string\n");
  return 0;
}

/*
 * Class:     appril_java_interface
 * Method:    StartQueryComputation
 * Signature: ([C)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_StartQueryComputation(JNIEnv *, jobject, jcharArray){

  printf( "map2string\n");
  return 0;
}

/*
 * Class:     appril_java_interface
 * Method:    GetSolution
 * Signature: ([CI)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_GetSolution(JNIEnv *, jobject, jcharArray, jint){
  printf( "map2string\n");
  return 0;
}

/*
 * Class:     appril_java_interface
 * Method:    StopComputation
 * Signature: ([CI)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_StopComputation(JNIEnv *, jobject, jcharArray, jint){
  printf( "map2string\n");
  return 0;
}
