/**
 * Created by Kathryn on 1/25/2015.
 */

import com.cra.figaro.language._
import com.cra.figaro.algorithm.factored.VariableElimination

object Main {

  def main(args:Array[String]): Unit ={
    val A = Flip(.3)
    val B =  Flip(.4)
    val C = Select(.2->A, .3 -> B, .5 -> true)

    val alg = VariableElimination(C)
    alg.start()
    val prob = alg.probability(C, true)
    alg.kill

    println("Probability of C being true: " + prob)






  }

}
