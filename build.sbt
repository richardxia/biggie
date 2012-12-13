name := "Biggie"

scalaSource in Compile := file("./src")

mainClass in (Compile, run) := Some("Biggie")

artifactPath in Compile in packageBin <<= baseDirectory { base => base / "build" / "biggie.jar" }

//libraryDependencies += "net.sf" % "picard" % "1.8.1" from "http://downloads.sourceforge.net/project/picard/sam-jdk/1.81/sam-1.81.jar"
