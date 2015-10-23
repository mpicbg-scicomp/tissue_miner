
import de.mpicbg.eaton.tissueminer.Tables
import de.mpicbg.eaton.tissueminer.Tables._
import slick.lifted.TableQuery

import scala.concurrent.Await
import slick.driver.SQLiteDriver.api._
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration.Duration._


/**
  * Document me!
  *
  * @author Holger Brandl
  */
object SqliteTest extends App {


  // do the same with slick
  //  SourceCodeGenerator.main(
  //    Array("slick.driver.SQLiteDriver", "org.sqlite.JDBC", "jdbc:sqlite:/Users/brandl/Desktop/demo_ForTissueMiner.sqlite", "/Users/brandl/projects/learning/flyscala/src/main/scala", "de.mpicbg.eaton.tissueminer")
  //  )

  // https://www.playframework.com/documentation/2.4.x/ScalaDatabase


  //  val h2mem1 = {
  //    url = "jdbc:sqlite:/Users/brandl/Desktop/demo_ForTissueMiner.sqlite"
  //    driver = org.sqlite.JDBC
  //    connectionPool = disabled
  //    keepAliveConnection = true
  //  }


  val db = Database.forURL("jdbc:sqlite:/Users/brandl/Desktop/demo_ForTissueMiner.sqlite", driver = "org.sqlite.JDBC")

  //  http://slick.typesafe.com/doc/3.1.0-RC3/gettingstarted.html#database-configuration
  //  db.createSession()

  val cellHistories = TableQuery[CellHistories]
  //  db.run(ch.result).map(println(_))

  //  val setup = DBIO.seq(
  //    (CellHistories.schema ++ Cells.schema).create
  //  )
  //
  //  val setupFuture = db.run(setup)
  //  Await.ready(setupFuture, 15.seconds)

  db.run(cellHistories.result).map(_.foreach {
    case (cellHistory: Tables.CellHistoriesRow) =>
      println(cellHistory.cellId)
  })

  // join a table (http://slick.typesafe.com/doc/3.1.0-RC3/gettingstarted.html#database-configuration)

  // This fails without type constraint because of
  // http://stackoverflow.com/questions/22204868/slick-2-0-define-generic-find-by-field-method

  //  val q2 = for {
  //    ch <- CellHistories if ch.generation < 9.0
  //    c <- Cells if ch.cellId === c.cellId
  //  } yield (c.cellId, c.frame, c.centerX, c.centerY)


  val q3 = for {
    c <- Cells
  } yield (c.cellId, c.frame, c.centerX, c.centerY)


  db.run(q3.result).map(_.foreach {
    case (cell) => println(cell)
  })

  import scala.concurrent.duration._


  // define custom query methods (ie DAO)
  // see http://stackoverflow.com/questions/31057153/slick-3-0-dao-fails-to-compile
  def getCellInfo(id: Int) = {
    val query = CellHistories.filter(_.cellId === id)
    Await.result(db.run(query.result), 10.hours).headOption
    db.run(query.result)
  }


  private val generation = getCellInfo(10629)


  //  def getCellInfo(id: Int) = CellHistories.map(_.cellId)


  //  db.run(getCellInfo(10313).result).result(10.hours) ==>  error: Don't call `Awaitable` methods directly, use the `Await` object.
  //  result.foreach(println(_))

  // see http://slick.typesafe.com/doc/3.0.0/sql-to-slick.html

  // plain queries are also possible (see http://slick.typesafe.com/doc/3.0.0/sql.html)
  val action = sql"SELECT cell_id, generation FROM cell_histories".as[(Int, Int)]
  db.run(action).foreach(println(_))
  //  ==> works!!
  // even better
  val action2 = sql"SELECT cell_id, generation FROM cell_histories".as[CellHistoriesRow]
  db.run(action2).foreach(println(_))


  //  intellij allows to define sql dialect per project, which considerable eases sql typing
  // also it allows to define a datasource which gives autocompletion on the db schema

  //  db.run(Tables.Vertices.result).map(_.foreach({
  //    case (frame, vertexId, xPos, yPos) =>
  //      println(frame + vertexId, xPos, yPos)

}




package object traversalTesting {
  val db = Database.forURL("jdbc:sqlite:/Users/brandl/Desktop/demo_ForTissueMiner.sqlite", driver = "org.sqlite.JDBC")

  //  http://slick.typesafe.com/doc/3.1.0-RC3/gettingstarted.html#database-configuration
  //  db.createSession()

  val cellHistories = TableQuery[CellHistories]

  implicit val cells: Seq[Tables.CellHistoriesRow] = Await.result(db.run(cellHistories.result), Inf)

  import de.mpicbg.eaton.tissueminer.TraversalUtils._

  val id = 10012
  val query = CellHistories.filter(_.cellId === id)
  val cellInfo = Await.result(db.run(query.result), Inf).head

  val daughter: Option[Tables.CellHistoriesRow] = cellInfo.leftDaughter
  val mother = cellInfo.mother

  //  cells.map(_.mother).flatten


}