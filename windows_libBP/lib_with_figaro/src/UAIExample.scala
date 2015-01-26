import com.cra.figaro.algorithm.factored.MPEVariableElimination
import com.cra.figaro.language._
import com.cra.figaro.library.compound._
import com.cra.figaro.algorithm._
import com.cra.figaro.algorithm.factored.factors.Variable

object UAIExample {

  /*
   * This is for the output of flip. 
   * In the internal list of values, 'false' happens to be at index 1 and 'true' is at index 0.
   */
  def convertBoolToIndex(b: Boolean) = {
    if (b == false) {
      1
    } else {
      0
    }
  }

  def TestProblem1 = {
    
    Universe.createNew
    val e1 = Flip(0.85)
    val e2 = Flip(0.25)
    val e3 = Flip(0.50)
    val e4 = If(e3, e1, e2)
    
    val e5 = Select(0.25 -> 0, 0.25 -> 1, .50 -> 2)
    val e6 = Select(0.15 -> 0, 0.75 -> 1, .10 -> 2)
    val e7 = If(e4, e5, e6)
    
    e4.observe(false)
    e7.observe(0)
    val algorithm = MPEVariableElimination()
    algorithm.start

    println(algorithm.mostLikelyValue(e1))
    println(algorithm.mostLikelyValue(e2))
    println(algorithm.mostLikelyValue(e3))
    println(algorithm.mostLikelyValue(e4))
    println(algorithm.mostLikelyValue(e5))
    println(algorithm.mostLikelyValue(e6))
    println(algorithm.mostLikelyValue(e7))

    val UAIConverter = UAIFormatConverter(e1, e2, e3, e4, e5, e6, e7)
    val result = UAIConverter.ConvertToUAIString
    val problem = result._1
    val evidence = result._2
    val variableMap = result._3
    println()
    println(problem)
    println(evidence)
    println()
    println(result._3(Variable(e1).id) + " " + convertBoolToIndex(algorithm.mostLikelyValue(e1)))
    println(result._3(Variable(e2).id) + " " + convertBoolToIndex(algorithm.mostLikelyValue(e2)))
    println(result._3(Variable(e3).id) + " " + convertBoolToIndex(algorithm.mostLikelyValue(e3)))
    println(result._3(Variable(e4).id) + " " + convertBoolToIndex(algorithm.mostLikelyValue(e4)))
    println(result._3(Variable(e5).id) + " " + algorithm.mostLikelyValue(e5))
    println(result._3(Variable(e6).id) + " " + algorithm.mostLikelyValue(e6))
    println(result._3(Variable(e7).id) + " " + algorithm.mostLikelyValue(e7))
    algorithm.stop
    algorithm.kill

    //                 0 1 2 3 4 5 6        
    //				   0 2 1 0 1 0 1
    //                 0 2 1 0 1 0 1


  }

  def TestProblem2 = {
    Universe.createNew
    val e1 = Flip(0.50)
    val e2 = Flip(1.0)
    val e3 = Flip(0.20)
    val e4 = If(e1, e2, e3)
    e1.observe(false)
    e4.observe(false)
    val algorithm = MPEVariableElimination()
    algorithm.start

    algorithm.stop
    algorithm.kill
    println(algorithm.mostLikelyValue(e1))
    println(algorithm.mostLikelyValue(e2))
    println(algorithm.mostLikelyValue(e3))
    println(algorithm.mostLikelyValue(e4))

    val UAIConverter = UAIFormatConverter(e1, e2, e3, e4)
    val result = UAIConverter.ConvertToUAIString
    val problem = result._1
    val evidence = result._2
    val variableMap = result._3
    println(problem)
    println(evidence)
    println(result._3(Variable(e1).id) + " " + convertBoolToIndex(algorithm.mostLikelyValue(e1)))
    println(result._3(Variable(e2).id) + " " + convertBoolToIndex(algorithm.mostLikelyValue(e2)))
    println(result._3(Variable(e3).id) + " " + convertBoolToIndex(algorithm.mostLikelyValue(e3)))
    println(result._3(Variable(e4).id) + " " + convertBoolToIndex(algorithm.mostLikelyValue(e4)))


    


  }

  def TestProblem3 = {
    Universe.createNew
    val e1 = Flip(0.436)
    val e2 = Flip(0.128)
    val e3 = Flip(0.920)
    val e4 = If(e1, e2, e3)
    val e5 = Select(0.210 -> 0, 0.333 -> 1, 0.457 -> 2)
    val e6 = Select(0.811 -> 0, 0.000 -> 1, 0.189 -> 2)
    val e7 = If(e4, e5, e6)

    e4.observe(true)
    e7.observe(1)
    println("Number of conditions: " + e5.allConditions.size)
    println("Number of constraints: " + e5.allConstraints.size)

    //Compare with Figaro algorithm:
    val algorithm = MPEVariableElimination()
    algorithm.start

    val UAIConverter = UAIFormatConverter(e1, e2, e3, e4, e5, e6, e7)
    val result = UAIConverter.ConvertToUAIString
    val problem = result._1
    val evidence = result._2
    val variableMap = result._3
    println(problem)
    println(evidence)

    println(result._3(Variable(e1).id) + " " + algorithm.mostLikelyValue(e1))
    println(result._3(Variable(e4).id) + " " + algorithm.mostLikelyValue(e4))
    println(result._3(Variable(e7).id) + " " + algorithm.mostLikelyValue(e7))
    println(result._3(Variable(e2).id) + " " + algorithm.mostLikelyValue(e2))
    println(result._3(Variable(e3).id) + " " + algorithm.mostLikelyValue(e3))
    println(result._3(Variable(e5).id) + " " + algorithm.mostLikelyValue(e5))
    println(result._3(Variable(e6).id) + " " + algorithm.mostLikelyValue(e6))

    println()
    println(result._3(Variable(e1).id) + " " + convertBoolToIndex(algorithm.mostLikelyValue(e1)))
    println(result._3(Variable(e4).id) + " " + convertBoolToIndex(algorithm.mostLikelyValue(e4)))
    println(result._3(Variable(e7).id) + " " + algorithm.mostLikelyValue(e7))

    println(result._3(Variable(e2).id) + " " + convertBoolToIndex(algorithm.mostLikelyValue(e2)))
    println(result._3(Variable(e3).id) + " " + convertBoolToIndex(algorithm.mostLikelyValue(e3)))
    println(result._3(Variable(e5).id) + " " + algorithm.mostLikelyValue(e5))
    println(result._3(Variable(e6).id) + " " + algorithm.mostLikelyValue(e6))

	//0 1 2 3 4 5 6
	//Figaro output: 1 1 1 0 0 1 0
    //   UCI output: 1 0 1 1 0 1 0
    //With 0.01 probability for e6:
    //   UCI output: 1 1 1 0 0 1 0
	
    algorithm.stop
    algorithm.kill
  }

  def TestProblem4 = {
    
    Universe.createNew()
    val e1 = Flip(0.5)
    e1.setConstraint((b: Boolean) => if (b) 3.0; else 1.0)
    val e2TrueOutcome = Flip(0.4)
    val e2FalseOutcome = Flip(0.9)
    val e2 = If(e1, Flip(0.4), Flip(0.9))

    val e3TrueOutcome = Flip(0.52)
    val e3FalseOutcome = Flip(0.4)
    val e3 = If(e1, e3TrueOutcome, e3FalseOutcome)

    val e4 = e2 === e3
    e4.observe(true)
    val alg = MPEVariableElimination()
    alg.start()
    println(alg.mostLikelyValue(e1))
    println(alg.mostLikelyValue(e2))
    println(alg.mostLikelyValue(e3))
    println(alg.mostLikelyValue(e4))

    val UAIConverter = UAIFormatConverter(e1, e2TrueOutcome, e2FalseOutcome, e2, e3TrueOutcome, e3FalseOutcome, e3, e4)
    val result = UAIConverter.ConvertToUAIString
    println(result._1) //Problem
    println(result._2) //Evidence
    println(result._3(Variable(e1).id) + " " + convertBoolToIndex(alg.mostLikelyValue(e1)))
    println(result._3(Variable(e2TrueOutcome).id) + " " + convertBoolToIndex(alg.mostLikelyValue(e2TrueOutcome)))
    println(result._3(Variable(e2FalseOutcome).id) + " " + convertBoolToIndex(alg.mostLikelyValue(e2FalseOutcome)))
    println(result._3(Variable(e2).id) + " " + convertBoolToIndex(alg.mostLikelyValue(e2)))
    println(result._3(Variable(e3TrueOutcome).id) + " " + convertBoolToIndex(alg.mostLikelyValue(e3TrueOutcome)))
    println(result._3(Variable(e3FalseOutcome).id) + " " + convertBoolToIndex(alg.mostLikelyValue(e3FalseOutcome)))
    println(result._3(Variable(e3).id) + " " + convertBoolToIndex(alg.mostLikelyValue(e3)))
    println(result._3(Variable(e4).id) + " " + convertBoolToIndex(alg.mostLikelyValue(e4)))

    //The === implicitly creates two more elements which we can't see here.
    //      0 1 2 3 4 5 6 7 8 9
    //      1 0 x 1 1 0 0 1 x 1
    //      1 0 0 1 1 0 0 1 1 1


    alg.kill();
  }

  def main(args: Array[String]): Unit = {

    //TestProblem1
    //TestProblem2
    TestProblem3
    //TestProblem4

  }

}
