package de.mpicbg.eaton.tissueminer.misc

import de.mpicbg.eaton.tissueminer.Tables
import de.mpicbg.eaton.tissueminer.Tables._
import de.mpicbg.eaton.tissueminer.TraversalUtils._
import slick.driver.SQLiteDriver.api._
import slick.lifted.TableQuery

import scala.concurrent.Await
import scala.concurrent.duration.Duration._


/**
  * Document me!
  *
  * @author Holger Brandl
  */
object TraversalTest {


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
