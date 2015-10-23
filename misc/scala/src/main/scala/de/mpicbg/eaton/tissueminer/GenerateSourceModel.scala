package de.mpicbg.eaton.tissueminer

import slick.codegen.SourceCodeGenerator


object GenerateSourceModel extends {

  SourceCodeGenerator.main(
      Array("slick.driver.SQLiteDriver", "org.sqlite.JDBC", "jdbc:sqlite:/Users/brandl/Desktop/demo_ForTissueMiner.sqlite", "/Users/brandl/projects/learning/flyscala/src/main/scala", "de.mpicbg.eaton.tissueminer")
    )

}