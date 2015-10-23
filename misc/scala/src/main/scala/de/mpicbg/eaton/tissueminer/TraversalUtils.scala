package de.mpicbg.eaton.tissueminer

/**
  * Implicit reference getters for cells and bonds
  *
  * @author Holger Brandl
  */
object TraversalUtils {


  implicit class Traversa(val cell: Tables.CellHistoriesRow) {

    def leftDaughter(implicit tissueHist: Seq[Tables.CellHistoriesRow]): Option[Tables.CellHistoriesRow] = {
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
