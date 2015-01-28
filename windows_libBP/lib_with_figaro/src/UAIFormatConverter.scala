import com.cra.figaro.algorithm.Values
import com.cra.figaro.algorithm.factored.factors.{Factory, Variable, Factor}
import com.cra.figaro.language._

import scala.collection.immutable.List
import scala.collection.immutable.Map
import scala.collection.mutable
import scala.collection.mutable.ListBuffer
import scala.util.Properties


class UAIFormatConverter(val universe: Universe, val targets: Element[_]*)(
  val showTiming: Boolean,
  val dependentUniverses: List[(Universe, List[NamedEvidence[_]])],
  val debug: Boolean) {
  val targetElements = targets.toList

//  val header = "BAYES"
val header = "MARKOV"
  protected var targetFactors: Map[Element[_], Factor[Double]] = Map()

  private def expand(): Unit = universe.activeElements

  private def getVariables(elements: Seq[Element[_]]): Map[Int, Variable[_]] = {
    val allElements = universe.activeElements
    val variables = allElements map (Variable(_))
    val result = mutable.Map.empty[Int, Variable[_]]
    var nextID = 0
    for (v <- variables) {
      result += nextID -> v
      nextID += 1
    }
    result.toMap
  }

  //The size and domain of the factor.
  def getFactorInformation(f: Factor[Double], variableOrdering: Map[Int, Int]): String = {
    val result = new StringBuilder

    result.append(f.variables.size + " ")
    val thisFactorOrdering = mutable.Map.empty[Int, Int]
    val listOfAbsoluteIDs = f.variables.map(v => variableOrdering(v.id))
    val orderedList = listOfAbsoluteIDs.sorted
    for ((v, i) <- f.variables.zipWithIndex) {
      thisFactorOrdering += i -> orderedList.indexOf(variableOrdering(v.id))
    }
    orderedList.foreach(x => result.append(x + " "))
    result.append(Properties.lineSeparator)
    result.toString
  }

  def compareLists(list1: List[Int], list2: List[Int], orderMap: Map[Int, Int]): Int = {
    //These lists must be the same length, since they are indices to the same factor.
    var counter = 0
    while (counter < list1.size) {
      val c = list1(orderMap(counter)) compare list2(orderMap(counter))
      if (c != 0) return c
      counter += 1
    }
    0 //Equal
  }

  //Variables must be consistently ordered throughout all factors
  def getUAIFactorString(f: Factor[Double], variableOrdering: Map[Int, Int]): String = {
    val result = new StringBuilder
    var newLine = 1

    for (v <- f.variables) {
      newLine *= v.range.size
    }

      newLine = newLine / (f.variables(f.variables.size - 1).range.size)
    val indices = f.allIndices

    //There are faster/more elegant ways of doing this in Scala
    val indexListBuffer = ListBuffer.empty[Int]
    for (x <- 0 to f.variables.size - 1) {
      indexListBuffer.append(x)
    }

    val indexList = indexListBuffer.toList

    val thisFactorOrdering = mutable.Map.empty[Int, Int]
    val listOfAbsoluteIDs = f.variables.map(v => variableOrdering(v.id))
    val orderedList = listOfAbsoluteIDs.sorted
    for ((v, i) <- f.variables.zipWithIndex) {
      thisFactorOrdering += i -> orderedList.indexOf(variableOrdering(v.id))
    }

    val indicesOrdering = new Ordering[List[Int]] { def compare(list1: List[Int], list2: List[Int]) = compareLists(list1, list2, thisFactorOrdering.toMap) }
    val allIndicesOrdered = indices.sorted(indicesOrdering)

    var counter = 0

    for (i <- allIndicesOrdered) {
      //i(0) corresponds to f.variables(0), i(1) corresponds to f.variables(1), etc.
      //So i.size == f.variables.size
      result.append(f.get(i) + " ")
      counter += 1
      if (f.variables.size > 1 && counter % newLine == 0) {
        result.append("\n")
      }

    }

    result.toString
  }

  def convertEvidenceToUCIFormat(elements: List[Element[_]], variableOrdering: Map[Int, Int]): String = {

    val result = StringBuilder.newBuilder
    val numberOfEvidenceLines = universe.conditionedElements.size + universe.constrainedElements.size
    result.append(numberOfEvidenceLines)
    result.append(Properties.lineSeparator)
    for (e <- elements) {
      val hasConditionsOrConstraints = !(e.allConditions.isEmpty && e.allConstraints.isEmpty)

      if (hasConditionsOrConstraints == true) {
        val v = Variable(e)
        
        result.append(variableOrdering(v.id) + " " + (v.range.indexOf(e.value)))
        result.append(Properties.lineSeparator)
      }
    }
    
    result.toString

  }

  def getFactorsWithID: Map[Int, Factor[Double]] = {
    val result = mutable.Map.empty[Int, Factor[Double]]
    val allElements = universe.activeElements
    var factorID = 0

    for (e <- allElements) {
      Values()(e)
      val factors = Factory.make(e)

      for (f <- factors) {
        result += factorID -> f
        factorID += 1
      }
    }
    result.toMap
  }

  def ConvertToUAIString(): (String, String,Map[Int, Int] ) = {
    expand()
    val problemString = StringBuilder.newBuilder
    problemString.append(header)
    problemString.append(Properties.lineSeparator)
    val targetVariables = getVariables(targetElements) //0 based index to variable
    val targetVariables2 = targetElements map (Variable(_))

    //Decide order
    val variablesToAbsoluteOrderPlacement = mutable.Map.empty[Int, Int] //Variable id to zero based index
    //Everything must use the same ordering
    targetVariables.map(v => variablesToAbsoluteOrderPlacement += v._2.id -> v._1)
    problemString.append(targetVariables.keys.size)
    problemString.append(Properties.lineSeparator)
    for (k <- 0 to targetVariables.keys.size - 1) {
      problemString.append(targetVariables(k).range.size + " ")
    }
    problemString.append(Properties.lineSeparator)
    val allFactors = getFactorsWithID
    problemString.append(allFactors.keys.size)
    problemString.append(Properties.lineSeparator)
    allFactors.keys.map(k => problemString.append(getFactorInformation(allFactors(k), variablesToAbsoluteOrderPlacement.toMap)))
    allFactors.keys.map(k => problemString.append(allFactors(k).allIndices.size + Properties.lineSeparator + getUAIFactorString(allFactors(k), variablesToAbsoluteOrderPlacement.toMap) + Properties.lineSeparator))

    val evidenceString = convertEvidenceToUCIFormat(targetElements,variablesToAbsoluteOrderPlacement.toMap)

    (problemString.toString, evidenceString.toString(),variablesToAbsoluteOrderPlacement.toMap)

  }

}

object UAIFormatConverter {

  def apply(targets: Element[_]*)(implicit universe: Universe) =
    new UAIFormatConverter(universe, targets: _*)(
      false,
      List(),
      false)

}
