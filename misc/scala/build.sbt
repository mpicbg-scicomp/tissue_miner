name := "tm_scala"

version := "0.1"

scalaVersion := "2.11.7"

//resolvers += "Typesafe repository" at "http://repo.typesafe.com/typesafe/releases/"

// http://www.scala-sbt.org/0.13/tutorial/Library-Dependencies.html
resolvers += Resolver.mavenLocal

// http://repo1.maven.org/maven2/org/squeryl/squeryl_2.11/
//libraryDependencies += "org.squeryl" %% "squeryl" % "0.9.6-SNAPSHOT"

// http://mvnrepository.com/artifact/com.h2database/h2/1.4.189

//libraryDependencies += "com.h2database" % "h2" % "1.4.189"


// recently added 0.9.6 contains it
libraryDependencies += "org.xerial" % "sqlite-jdbc" % "3.8.11.2"


// http://mvnrepository.com/artifact/com.typesafe.slick/slick_2.11/3.1.0-M2
libraryDependencies += "com.typesafe.slick" %% "slick" % "3.1.0-RC3"
libraryDependencies += "com.typesafe.slick" %% "slick-codegen" % "3.1.0-RC3"

libraryDependencies += "org.slf4j" % "slf4j-nop" % "1.6.4"


libraryDependencies += "com.lihaoyi" % "ammonite-repl" % "0.4.9-SNAPSHOT" % "test" cross CrossVersion.full
initialCommands in(Test, console) := """ammonite.repl.Repl.run("")"""


//libraryDependencies  ++=  Seq(
//  "org.squeryl" %% "squeryl" % "0.9.5-6",
////  "com.typesafe.slick" %% "slick" % "3.0.3",
//"org.xerial" % "sqlite-jdbc" % "3.7.2"
//
//)
