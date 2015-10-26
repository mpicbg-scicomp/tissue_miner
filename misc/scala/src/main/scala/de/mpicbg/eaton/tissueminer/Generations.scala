package de.mpicbg.eaton.tissueminer

import de.mpicbg.eaton.tissueminer.Tables._
import de.mpicbg.eaton.tissueminer.TraversalUtils._
import de.mpicbg.eaton.tissueminer._
import slick.driver.SQLiteDriver.api._
import slick.lifted.TableQuery

import scala.collection.parallel.immutable.ParMap
import scala.concurrent.Await
import scala.concurrent.duration.Duration._


/**
  * Reimplementation of generation calculation
  *
  * @author Holger Brandl
  */
object Generations extends App {

  // define a db context to work with
  val db = Database.forURL("jdbc:sqlite:/Users/brandl/Desktop/demo_ForTissueMiner.sqlite", driver = "org.sqlite.JDBC")
  implicit val cells: Seq[Tables.CellHistoriesRow] = Await.result(db.run(TableQuery[CellHistories].result), Inf)


  case class GenCounter(cell: Tables.CellHistoriesRow, generation: Int)

  // define a recursive function to calculate the cell generation
  def calcGeneration(genCounter: GenCounter): GenCounter = {
    genCounter.cell.mother match {
      case Some(motherCell) => calcGeneration(GenCounter(motherCell, genCounter.generation + 1))
      case None => genCounter
    }
  }


  val generations = cells.par.map(cell => cell -> calcGeneration(GenCounter(cell, 0)).generation).toMap

  // calculate a histogram (http://langref.org/scala/maps/algorithms/histogram)
  generations.groupBy(_._2).mapValues(_.size).mkString("\n")
}
