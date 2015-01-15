#include "bp_interface.h"
#include "jni.h"
#include "stdio.h"
#include "BPinterface.h"



/*
 * Class:     bp_interface
 * Method:    doBP
 * Signature: (Ljava/lang/String;)Z
 */
JNIEXPORT jboolean JNICALL Java_bp_1interface_doBP (JNIEnv * env, jobject job, jstring fileName){
printf("doBP\n");
const char *str= env->GetStringUTFChars(fileName, NULL);
int len = env->GetStringLength(fileName);

char file[len];
for(int i=0; i < len; i++){
file[i] = str[i];
}

  printf("The string is %s", file);

lgbp::BpInterface* alg = new lgbp::BpInterface();
 
char task[] = {'P','R'};
  alg->initialize(60, task, file, NULL, NULL, true);
  alg->runInference();
 

  //need to release this string when done with it in
  //order to avoid memory leak
  env->ReleaseStringUTFChars(fileName, str);
 
  return true;
}

/*
 * Class:     bp_interface
 * Method:    getSolution
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_bp_1interface_getSolution (JNIEnv * env, jobject jobj){
  printf("getSolution\n");
  double res;
 //  alg->getSolution(res);
// FILE* output = fopen("bp_results.txt", "w");
// fprintf(output, "%g\n", res);
// fclose(output);
  return env->NewStringUTF(":O");
  
}
