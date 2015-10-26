package de.mpicbg.eaton.tissueminer

import Generations
import de.mpicbg.eaton.tissueminer.Tables
import de.mpicbg.eaton.tissueminer.Tables._
import slick.driver.SQLiteDriver.api._
import slick.lifted.TableQuery

import scala.concurrent.Await
import scala.concurrent.duration.Duration._
import de.mpicbg.eaton.tissueminer.TraversalUtils._


/**
  * Implicit reference getters for cells and bonds
  *
  * @author Holger Brandl
  */
object TraversalUtils {


  implicit class Traversal(val cell: Tables.CellHistoriesRow) {

    def leftDaughter(implicit tissueHist: Seq[Tables.CellHistoriesRow]) = {
      cell.leftDaughterCellId match {
        case Some(id) => tissueHist.find(_.cellId == id)
        case None => None
      }
    }


    def rightDaughter(implicit tissueHist: Seq[Tables.CellHistoriesRow]): Option[Tables.CellHistoriesRow] = {
      cell.rightDaughterCellId match {
        case Some(id) => tissueHist.find(_.cellId == id)
        case None => None
      }
    }


    def mother(implicit tissueHist: Seq[Tables.CellHistoriesRow]): Option[Tables.CellHistoriesRow] = {
      tissueHist.find(someCell => {
        //        val children = List(someCell.leftDaughterCellId.getOrElse(-1), someCell.leftDaughterCellId.getOrElse(-1))
        val children = List(someCell.leftDaughterCellId.getOrElse(-1), someCell.leftDaughterCellId.getOrElse(-1))
        children.contains(cell.cellId)
      })
    }
  }
}


package object TraversalTest {


  val db = Database.forURL("jdbc:sqlite:/Users/brandl/Desktop/demo_ForTissueMiner.sqlite", driver = "org.sqlite.JDBC")

  //  http://slick.typesafe.com/doc/3.1.0-RC3/gettingstarted.html#database-configuration
  //  db.createSession()

  val cellHistories = TableQuery[CellHistories]

  implicit val cells: Seq[Tables.CellHistoriesRow] = Await.result(db.run(cellHistories.result), Inf)


  val id = 10012
  val query = CellHistories.filter(_.cellId === id)
  val cellInfo = Await.result(db.run(query.result), Inf).head

  val daughter = cellInfo.leftDaughter
  val mother = cellInfo.mother
}
