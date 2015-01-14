#include "bp_interface.h"
#include "jni.h"
#include "stdio.h"

/*
 * Class:     bp_interface
 * Method:    doBP
 * Signature: (Ljava/lang/String;)Z
 */
JNIEXPORT jboolean JNICALL Java_bp_1interface_doBP (JNIEnv * env, jobject job, jstring fileName){
  printf("doBP\n");

  return true;
}

/*
 * Class:     bp_interface
 * Method:    getSolution
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_bp_1interface_getSolution (JNIEnv * env, jobject jobj){
  printf("getSolution\n");
  return env->NewStringUTF(":D");
  
}
