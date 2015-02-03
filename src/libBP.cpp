#include "bp_interface.h"
#include "stdio.h"
#include "BPinterface.h"

/*
 * Class:     bp_interface
 * Method:    initBpInterface
 * Signature: ()J
 */
JNIEXPORT jlong JNICALL Java_bp_1interface_initBpInterface (JNIEnv *env, jobject jobj){
  jfieldID fid;
  jclass cls;
printf("initBpInterface");
  cls = env->GetObjectClass(jobj);
  fid = env->GetFieldID(cls, "interfacePtr", "J"); // get the pointer to the bp object



  lgbp::BpInterface* alg = new lgbp::BpInterface();
    
  env->SetLongField(jobj, fid, (jlong)alg);
  return 0;
   
}

/*
 * Class:     bp_interface
 * Method:    destroyBpInterface
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_bp_1interface_destroyBpInterface (JNIEnv * env , jobject jobj){
  jfieldID fid;
  jclass cls;
  lgbp::BpInterface* alg;
  printf("destroyBpInterfact()\n");
  cls = env->GetObjectClass(jobj);
  fid = env->GetFieldID(cls, "interfacePtr", "J"); // get the pointer to the bp object
  alg = (lgbp::BpInterface*)env->GetLongField(jobj,fid);
    
  delete alg;
}

/*
 * Class:     bp_interface
 * Method:    doBP
 * Signature: (Ljava/lang/String;)Z
 */
JNIEXPORT jboolean JNICALL Java_bp_1interface_doBP (JNIEnv * env, jobject jobj, jstring fileName){

  // Get the bpObject
  jclass cls = env->GetObjectClass(jobj);
  jfieldID fid = env->GetFieldID(cls, "interfacePtr", "L");
  lgbp::BpInterface* alg = (lgbp::BpInterface*)env->GetLongField(jobj,fid);
  


  printf("\n\ndoBP\n");
  const char *str= env->GetStringUTFChars(fileName, NULL);
  int len = env->GetStringLength(fileName);

  char file[len];
  printf("\n\n");
  for(int i=0; i < len; i++){
    file[i] = str[i];
    printf("%c  ", (char)str[i]);
  }
  //printf("\n\n");

  printf("The string is %s\n\n", file);
  printf("The string is %s\n\n", str);

 
  char task[2];
  task[0] = 'P';
  task[1] = 'R';
  printf("Initializing\n");
  alg->initialize(60, task, file, NULL, NULL, true);
  printf("Running inference \n");
  alg->runInference();
 

  //need to release this string when done with it in
  //order to avoid memory leak
  env->ReleaseStringUTFChars(fileName, str);

 printf("Ending doBP()\n");
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
  //   alg->getSolution(res);
  // FILE* output = fopen("bp_results.txt", "w");
  // fprintf(output, "%g\n", res);
  // fclose(output);
  return env->NewStringUTF(":O");
  
}

