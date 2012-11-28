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

sbt "run-main biggie.SnpCaller samFile rangeStart rangeEnd"
