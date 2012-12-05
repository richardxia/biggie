cs262a
======

A pipeline for distributed SNP calling

Requirements
------------
* Scala
* sbt

Compiling
---------
sbt compile

Running
-------
./bin/biggie cmd args,...

sbt "run-main biggie.SnpCaller alignment.bam ref.fa regions.txt"

regions.txt:
refSeq	start	end

e.g.
chr21	14800232	14804369
