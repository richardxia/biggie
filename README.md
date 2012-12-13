cs262a
======

A pipeline for distributed SNP calling

Requirements
------------
* Scala
* sbt

Dependencies
------------
* Picard library (I have a separate script which fetches it automatically (`download_picard.sh`), but if not, place sam-<version>.jar and picard-<version>.jar in lib/)

Compiling
---------
In the root directory:

  sbt compile

Running
-------
./bin/biggie cmd args,...

    sbt "run-main biggie.SimpleClassifier alignment.bam number_of_bases"
    sbt "run-main biggie.SnpCaller alignment.bam ref.fa regions.txt"

running directly:

    JAVA_OPTS=-Xmx16g -cp lib/sam-1.8.1.jar:target/scala-2.9.2/classes biggie.SnpCaller alignment.sam number_of_bases

regions.txt:
refSeq	start	end

e.g.
chr21	14800232	14804369
