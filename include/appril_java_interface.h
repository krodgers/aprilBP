/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class appril_java_interface */

#ifndef _Included_appril_java_interface
#define _Included_appril_java_interface
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     appril_java_interface
 * Method:    MapErrorCode2String
 * Signature: (I)[C
 */
JNIEXPORT jcharArray JNICALL Java_appril_1java_1interface_MapErrorCode2String
  (JNIEnv *, jobject, jint);

/*
 * Class:     appril_java_interface
 * Method:    Initialize
 * Signature: (Z)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_Initialize__Z
  (JNIEnv *, jobject, jboolean);

/*
 * Class:     appril_java_interface
 * Method:    Initialize
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_Initialize__
  (JNIEnv *, jobject);

/*
 * Class:     appril_java_interface
 * Method:    Terminate
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_Terminate
  (JNIEnv *, jobject);

/*
 * Class:     appril_java_interface
 * Method:    GetCurrentState
 * Signature: ([CI)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_GetCurrentState
  (JNIEnv *, jobject, jcharArray, jint);

/*
 * Class:     appril_java_interface
 * Method:    DefineProblem
 * Signature: ([C[C[C)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_DefineProblem
  (JNIEnv *, jobject, jcharArray, jcharArray, jcharArray);

/*
 * Class:     appril_java_interface
 * Method:    GetProblemComplexity
 * Signature: ([CI)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_GetProblemComplexity
  (JNIEnv *, jobject, jcharArray, jint);

/*
 * Class:     appril_java_interface
 * Method:    GetSolutionComplexity
 * Signature: ([CI)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_GetSolutionComplexity
  (JNIEnv *, jobject, jcharArray, jint);

/*
 * Class:     appril_java_interface
 * Method:    StartQueryComputation
 * Signature: ([C)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_StartQueryComputation
  (JNIEnv *, jobject, jcharArray);

/*
 * Class:     appril_java_interface
 * Method:    GetSolution
 * Signature: ([CI)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_GetSolution
  (JNIEnv *, jobject, jcharArray, jint);

/*
 * Class:     appril_java_interface
 * Method:    StopComputation
 * Signature: ([CI)I
 */
JNIEXPORT jint JNICALL Java_appril_1java_1interface_StopComputation
  (JNIEnv *, jobject, jcharArray, jint);

#ifdef __cplusplus
}
#endif
#endif
