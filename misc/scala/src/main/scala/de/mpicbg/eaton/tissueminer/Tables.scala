package de.mpicbg.eaton.tissueminer

// AUTO-GENERATED Slick data model
/** Stand-alone Slick data model for immediate use */
object Tables extends {

  val profile = slick.driver.SQLiteDriver
} with Tables


/** Slick data model trait for extension, choice of backend or usage in the cake pattern. (Make sure to initialize this late.) */
trait Tables {

  val profile: slick.driver.JdbcProfile

  import profile.api._
  import slick.model.ForeignKeyAction

  // NOTE: GetResult mappers for plain SQL are only generated for tables where Slick knows how to map the types of all columns.

  import slick.jdbc.{GetResult => GR}

  /** DDL for all tables. Call .create to execute. */
  lazy val schema: profile.SchemaDescription = Array(Bonds.schema, CellHistories.schema, Cells.schema, DirectedBonds.schema, Frames.schema, Vertices.schema).reduceLeft(_ ++ _)

  @deprecated("Use .schema instead of .ddl", "3.0")
  def ddl = schema


  /** Entity class storing rows of table Bonds
    * @param frame Database column frame SqlType(INTEGER)
    * @param bondId Database column bond_id SqlType(INTEGER), PrimaryKey
    * @param bondLength Database column bond_length SqlType(REAL) */
  case class BondsRow(frame: Int, bondId: Int, bondLength: Double)


  /** GetResult implicit for fetching BondsRow objects using plain SQL queries */
  implicit def GetResultBondsRow(implicit e0: GR[Int], e1: GR[Double]): GR[BondsRow] = GR {
    prs => import prs._
      BondsRow.tupled((<<[Int], <<[Int], <<[Double]))
  }


  /** Table description of table bonds. Objects of this class serve as prototypes for rows in queries. */
  class Bonds(_tableTag: Tag) extends Table[BondsRow](_tableTag, "bonds") {

    def * = (frame, bondId, bondLength) <>(BondsRow.tupled, BondsRow.unapply)


    /** Maps whole row to an option. Useful for outer joins. */
    def ? = (Rep.Some(frame), Rep.Some(bondId), Rep.Some(bondLength)).shaped.<>({ r => import r._; _1.map(_ => BondsRow.tupled((_1.get, _2.get, _3.get))) }, (_: Any) => throw new Exception("Inserting into ? projection not supported."))


    /** Database column frame SqlType(INTEGER) */
    val frame: Rep[Int] = column[Int]("frame")
    /** Database column bond_id SqlType(INTEGER), PrimaryKey */
    val bondId: Rep[Int] = column[Int]("bond_id", O.PrimaryKey)
    /** Database column bond_length SqlType(REAL) */
    val bondLength: Rep[Double] = column[Double]("bond_length")

    /** Foreign key referencing Frames (database name frames_FK_1) */
    lazy val framesFk = foreignKey("frames_FK_1", frame, Frames)(r => r.frame, onUpdate = ForeignKeyAction.NoAction, onDelete = ForeignKeyAction.NoAction)
  }


  /** Collection-like TableQuery object for table Bonds */
  lazy val Bonds = new TableQuery(tag => new Bonds(tag))


  /** Entity class storing rows of table CellHistories
    * @param cellId Database column cell_id SqlType(INTEGER), PrimaryKey
    * @param tissueAnalyzerGroupId Database column tissue_analyzer_group_id SqlType(INTEGER)
    * @param firstOcc Database column first_occ SqlType(INTEGER)
    * @param lastOcc Database column last_occ SqlType(INTEGER)
    * @param leftDaughterCellId Database column left_daughter_cell_id SqlType(INTEGER)
    * @param rightDaughterCellId Database column right_daughter_cell_id SqlType(INTEGER)
    * @param appearsBy Database column appears_by SqlType(TEXT)
    * @param disappearsBy Database column disappears_by SqlType(TEXT)
    * @param lineageGroup Database column lineage_group SqlType(TEXT)
    * @param generation Database column generation SqlType(INTEGER) */
  case class CellHistoriesRow(cellId: Int, tissueAnalyzerGroupId: Int, firstOcc: Int, lastOcc: Int, leftDaughterCellId: Option[Int], rightDaughterCellId: Option[Int], appearsBy: String, disappearsBy: String, lineageGroup: String, generation: Option[Int])


  /** GetResult implicit for fetching CellHistoriesRow objects using plain SQL queries */
  implicit def GetResultCellHistoriesRow(implicit e0: GR[Int], e1: GR[Option[Int]], e2: GR[String]): GR[CellHistoriesRow] = GR {
    prs => import prs._
      CellHistoriesRow.tupled((<<[Int], <<[Int], <<[Int], <<[Int], <<?[Int], <<?[Int], <<[String], <<[String], <<[String], <<?[Int]))
  }


  /** Table description of table cell_histories. Objects of this class serve as prototypes for rows in queries. */
  class CellHistories(_tableTag: Tag) extends Table[CellHistoriesRow](_tableTag, "cell_histories") {

    def * = (cellId, tissueAnalyzerGroupId, firstOcc, lastOcc, leftDaughterCellId, rightDaughterCellId, appearsBy, disappearsBy, lineageGroup, generation) <>(CellHistoriesRow.tupled, CellHistoriesRow.unapply)


    /** Maps whole row to an option. Useful for outer joins. */
    def ? = (Rep.Some(cellId), Rep.Some(tissueAnalyzerGroupId), Rep.Some(firstOcc), Rep.Some(lastOcc), leftDaughterCellId, rightDaughterCellId, Rep.Some(appearsBy), Rep.Some(disappearsBy), Rep.Some(lineageGroup), generation).shaped.<>({ r => import r._; _1.map(_ => CellHistoriesRow.tupled((_1.get, _2.get, _3.get, _4.get, _5, _6, _7.get, _8.get, _9.get, _10))) }, (_: Any) => throw new Exception("Inserting into ? projection not supported."))


    /** Database column cell_id SqlType(INTEGER), PrimaryKey */
    val cellId: Rep[Int] = column[Int]("cell_id", O.PrimaryKey)
    /** Database column tissue_analyzer_group_id SqlType(INTEGER) */
    val tissueAnalyzerGroupId: Rep[Int] = column[Int]("tissue_analyzer_group_id")
    /** Database column first_occ SqlType(INTEGER) */
    val firstOcc: Rep[Int] = column[Int]("first_occ")
    /** Database column last_occ SqlType(INTEGER) */
    val lastOcc: Rep[Int] = column[Int]("last_occ")
    /** Database column left_daughter_cell_id SqlType(INTEGER) */
    val leftDaughterCellId: Rep[Option[Int]] = column[Option[Int]]("left_daughter_cell_id")
    /** Database column right_daughter_cell_id SqlType(INTEGER) */
    val rightDaughterCellId: Rep[Option[Int]] = column[Option[Int]]("right_daughter_cell_id")
    /** Database column appears_by SqlType(TEXT) */
    val appearsBy: Rep[String] = column[String]("appears_by")
    /** Database column disappears_by SqlType(TEXT) */
    val disappearsBy: Rep[String] = column[String]("disappears_by")
    /** Database column lineage_group SqlType(TEXT) */
    val lineageGroup: Rep[String] = column[String]("lineage_group")
    /** Database column generation SqlType(INTEGER) */
    val generation: Rep[Option[Int]] = column[Option[Int]]("generation")

    /** Foreign key referencing CellHistories (database name cell_histories_FK_1) */
    lazy val cellHistoriesFk = foreignKey("cell_histories_FK_1", (rightDaughterCellId, leftDaughterCellId), CellHistories)(r => (Rep.Some(r.cellId), Rep.Some(r.cellId)), onUpdate = ForeignKeyAction.NoAction, onDelete = ForeignKeyAction.NoAction)
    /** Foreign key referencing Frames (database name frames_FK_2) */
    lazy val framesFk = foreignKey("frames_FK_2", (lastOcc, firstOcc), Frames)(r => (r.frame, r.frame), onUpdate = ForeignKeyAction.NoAction, onDelete = ForeignKeyAction.NoAction)
  }


  /** Collection-like TableQuery object for table CellHistories */
  lazy val CellHistories = new TableQuery(tag => new CellHistories(tag))


  /** Entity class storing rows of table Cells
    * @param frame Database column frame SqlType(INTEGER)
    * @param cellId Database column cell_id SqlType(INTEGER)
    * @param centerX Database column center_x SqlType(REAL)
    * @param centerY Database column center_y SqlType(REAL)
    * @param area Database column area SqlType(REAL)
    * @param elongXx Database column elong_xx SqlType(REAL)
    * @param elongXy Database column elong_xy SqlType(REAL) */
  case class CellsRow(frame: Int, cellId: Int, centerX: Double, centerY: Double, area: Double, elongXx: Double, elongXy: Double)


  /** GetResult implicit for fetching CellsRow objects using plain SQL queries */
  implicit def GetResultCellsRow(implicit e0: GR[Int], e1: GR[Double]): GR[CellsRow] = GR {
    prs => import prs._
      CellsRow.tupled((<<[Int], <<[Int], <<[Double], <<[Double], <<[Double], <<[Double], <<[Double]))
  }


  /** Table description of table cells. Objects of this class serve as prototypes for rows in queries. */
  class Cells(_tableTag: Tag) extends Table[CellsRow](_tableTag, "cells") {

    def * = (frame, cellId, centerX, centerY, area, elongXx, elongXy) <>(CellsRow.tupled, CellsRow.unapply)


    /** Maps whole row to an option. Useful for outer joins. */
    def ? = (Rep.Some(frame), Rep.Some(cellId), Rep.Some(centerX), Rep.Some(centerY), Rep.Some(area), Rep.Some(elongXx), Rep.Some(elongXy)).shaped.<>({ r => import r._; _1.map(_ => CellsRow.tupled((_1.get, _2.get, _3.get, _4.get, _5.get, _6.get, _7.get))) }, (_: Any) => throw new Exception("Inserting into ? projection not supported."))


    /** Database column frame SqlType(INTEGER) */
    val frame: Rep[Int] = column[Int]("frame")
    /** Database column cell_id SqlType(INTEGER) */
    val cellId: Rep[Int] = column[Int]("cell_id")
    /** Database column center_x SqlType(REAL) */
    val centerX: Rep[Double] = column[Double]("center_x")
    /** Database column center_y SqlType(REAL) */
    val centerY: Rep[Double] = column[Double]("center_y")
    /** Database column area SqlType(REAL) */
    val area: Rep[Double] = column[Double]("area")
    /** Database column elong_xx SqlType(REAL) */
    val elongXx: Rep[Double] = column[Double]("elong_xx")
    /** Database column elong_xy SqlType(REAL) */
    val elongXy: Rep[Double] = column[Double]("elong_xy")

    /** Primary key of Cells (database name pk_cells) */
    val pk = primaryKey("pk_cells", (cellId, frame))

    /** Foreign key referencing CellHistories (database name cell_histories_FK_1) */
    lazy val cellHistoriesFk = foreignKey("cell_histories_FK_1", cellId, CellHistories)(r => r.cellId, onUpdate = ForeignKeyAction.NoAction, onDelete = ForeignKeyAction.NoAction)
    /** Foreign key referencing Frames (database name frames_FK_2) */
    lazy val framesFk = foreignKey("frames_FK_2", frame, Frames)(r => r.frame, onUpdate = ForeignKeyAction.NoAction, onDelete = ForeignKeyAction.NoAction)
  }


  /** Collection-like TableQuery object for table Cells */
  lazy val Cells = new TableQuery(tag => new Cells(tag))


  /** Entity class storing rows of table DirectedBonds
    * @param frame Database column frame SqlType(INTEGER)
    * @param cellId Database column cell_id SqlType(INTEGER)
    * @param dbondId Database column dbond_id SqlType(INTEGER), PrimaryKey
    * @param conjDbondId Database column conj_dbond_id SqlType(INTEGER)
    * @param bondId Database column bond_id SqlType(INTEGER)
    * @param vertexId Database column vertex_id SqlType(INTEGER)
    * @param leftDbondId Database column left_dbond_id SqlType(INTEGER) */
  case class DirectedBondsRow(frame: Int, cellId: Int, dbondId: Int, conjDbondId: Int, bondId: Int, vertexId: Int, leftDbondId: Int)


  /** GetResult implicit for fetching DirectedBondsRow objects using plain SQL queries */
  implicit def GetResultDirectedBondsRow(implicit e0: GR[Int]): GR[DirectedBondsRow] = GR {
    prs => import prs._
      DirectedBondsRow.tupled((<<[Int], <<[Int], <<[Int], <<[Int], <<[Int], <<[Int], <<[Int]))
  }


  /** Table description of table directed_bonds. Objects of this class serve as prototypes for rows in queries. */
  class DirectedBonds(_tableTag: Tag) extends Table[DirectedBondsRow](_tableTag, "directed_bonds") {

    def * = (frame, cellId, dbondId, conjDbondId, bondId, vertexId, leftDbondId) <>(DirectedBondsRow.tupled, DirectedBondsRow.unapply)


    /** Maps whole row to an option. Useful for outer joins. */
    def ? = (Rep.Some(frame), Rep.Some(cellId), Rep.Some(dbondId), Rep.Some(conjDbondId), Rep.Some(bondId), Rep.Some(vertexId), Rep.Some(leftDbondId)).shaped.<>({ r => import r._; _1.map(_ => DirectedBondsRow.tupled((_1.get, _2.get, _3.get, _4.get, _5.get, _6.get, _7.get))) }, (_: Any) => throw new Exception("Inserting into ? projection not supported."))


    /** Database column frame SqlType(INTEGER) */
    val frame: Rep[Int] = column[Int]("frame")
    /** Database column cell_id SqlType(INTEGER) */
    val cellId: Rep[Int] = column[Int]("cell_id")
    /** Database column dbond_id SqlType(INTEGER), PrimaryKey */
    val dbondId: Rep[Int] = column[Int]("dbond_id", O.PrimaryKey)
    /** Database column conj_dbond_id SqlType(INTEGER) */
    val conjDbondId: Rep[Int] = column[Int]("conj_dbond_id")
    /** Database column bond_id SqlType(INTEGER) */
    val bondId: Rep[Int] = column[Int]("bond_id")
    /** Database column vertex_id SqlType(INTEGER) */
    val vertexId: Rep[Int] = column[Int]("vertex_id")
    /** Database column left_dbond_id SqlType(INTEGER) */
    val leftDbondId: Rep[Int] = column[Int]("left_dbond_id")

    /** Foreign key referencing Bonds (database name bonds_FK_1) */
    lazy val bondsFk = foreignKey("bonds_FK_1", bondId, Bonds)(r => r.bondId, onUpdate = ForeignKeyAction.NoAction, onDelete = ForeignKeyAction.NoAction)
    /** Foreign key referencing CellHistories (database name cell_histories_FK_2) */
    lazy val cellHistoriesFk = foreignKey("cell_histories_FK_2", cellId, CellHistories)(r => r.cellId, onUpdate = ForeignKeyAction.NoAction, onDelete = ForeignKeyAction.NoAction)
    /** Foreign key referencing DirectedBonds (database name directed_bonds_FK_3) */
    lazy val directedBondsFk = foreignKey("directed_bonds_FK_3", (leftDbondId, conjDbondId), DirectedBonds)(r => (r.dbondId, r.dbondId), onUpdate = ForeignKeyAction.NoAction, onDelete = ForeignKeyAction.NoAction)
    /** Foreign key referencing Frames (database name frames_FK_4) */
    lazy val framesFk = foreignKey("frames_FK_4", frame, Frames)(r => r.frame, onUpdate = ForeignKeyAction.NoAction, onDelete = ForeignKeyAction.NoAction)
    /** Foreign key referencing Vertices (database name vertices_FK_5) */
    lazy val verticesFk = foreignKey("vertices_FK_5", vertexId, Vertices)(r => r.vertexId, onUpdate = ForeignKeyAction.NoAction, onDelete = ForeignKeyAction.NoAction)

    /** Uniqueness Index over (conjDbondId) (database name sqlite_autoindex_directed_bonds_1) */
    val index1 = index("sqlite_autoindex_directed_bonds_1", conjDbondId, unique = true)
    /** Uniqueness Index over (leftDbondId) (database name sqlite_autoindex_directed_bonds_2) */
    val index2 = index("sqlite_autoindex_directed_bonds_2", leftDbondId, unique = true)
  }


  /** Collection-like TableQuery object for table DirectedBonds */
  lazy val DirectedBonds = new TableQuery(tag => new DirectedBonds(tag))


  /** Entity class storing rows of table Frames
    * @param frame Database column frame SqlType(INTEGER), PrimaryKey
    * @param timeSec Database column time_sec SqlType(INTEGER) */
  case class FramesRow(frame: Int, timeSec: Int)


  /** GetResult implicit for fetching FramesRow objects using plain SQL queries */
  implicit def GetResultFramesRow(implicit e0: GR[Int]): GR[FramesRow] = GR {
    prs => import prs._
      FramesRow.tupled((<<[Int], <<[Int]))
  }


  /** Table description of table frames. Objects of this class serve as prototypes for rows in queries. */
  class Frames(_tableTag: Tag) extends Table[FramesRow](_tableTag, "frames") {

    def * = (frame, timeSec) <>(FramesRow.tupled, FramesRow.unapply)


    /** Maps whole row to an option. Useful for outer joins. */
    def ? = (Rep.Some(frame), Rep.Some(timeSec)).shaped.<>({ r => import r._; _1.map(_ => FramesRow.tupled((_1.get, _2.get))) }, (_: Any) => throw new Exception("Inserting into ? projection not supported."))


    /** Database column frame SqlType(INTEGER), PrimaryKey */
    val frame: Rep[Int] = column[Int]("frame", O.PrimaryKey)
    /** Database column time_sec SqlType(INTEGER) */
    val timeSec: Rep[Int] = column[Int]("time_sec")
  }


  /** Collection-like TableQuery object for table Frames */
  lazy val Frames = new TableQuery(tag => new Frames(tag))


  /** Entity class storing rows of table Vertices
    * @param frame Database column frame SqlType(INTEGER)
    * @param vertexId Database column vertex_id SqlType(INTEGER), PrimaryKey
    * @param xPos Database column x_pos SqlType(REAL)
    * @param yPos Database column y_pos SqlType(REAL) */
  case class VerticesRow(frame: Int, vertexId: Int, xPos: Double, yPos: Double)


  /** GetResult implicit for fetching VerticesRow objects using plain SQL queries */
  implicit def GetResultVerticesRow(implicit e0: GR[Int], e1: GR[Double]): GR[VerticesRow] = GR {
    prs => import prs._
      VerticesRow.tupled((<<[Int], <<[Int], <<[Double], <<[Double]))
  }


  /** Table description of table vertices. Objects of this class serve as prototypes for rows in queries. */
  class Vertices(_tableTag: Tag) extends Table[VerticesRow](_tableTag, "vertices") {

    def * = (frame, vertexId, xPos, yPos) <>(VerticesRow.tupled, VerticesRow.unapply)


    /** Maps whole row to an option. Useful for outer joins. */
    def ? = (Rep.Some(frame), Rep.Some(vertexId), Rep.Some(xPos), Rep.Some(yPos)).shaped.<>({ r => import r._; _1.map(_ => VerticesRow.tupled((_1.get, _2.get, _3.get, _4.get))) }, (_: Any) => throw new Exception("Inserting into ? projection not supported."))


    /** Database column frame SqlType(INTEGER) */
    val frame: Rep[Int] = column[Int]("frame")
    /** Database column vertex_id SqlType(INTEGER), PrimaryKey */
    val vertexId: Rep[Int] = column[Int]("vertex_id", O.PrimaryKey)
    /** Database column x_pos SqlType(REAL) */
    val xPos: Rep[Double] = column[Double]("x_pos")
    /** Database column y_pos SqlType(REAL) */
    val yPos: Rep[Double] = column[Double]("y_pos")

    /** Foreign key referencing Frames (database name frames_FK_1) */
    lazy val framesFk = foreignKey("frames_FK_1", frame, Frames)(r => r.frame, onUpdate = ForeignKeyAction.NoAction, onDelete = ForeignKeyAction.NoAction)
  }


  /** Collection-like TableQuery object for table Vertices */
  lazy val Vertices = new TableQuery(tag => new Vertices(tag))
}
