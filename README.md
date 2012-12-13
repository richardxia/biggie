cs262a
======

A pipeline for distributed SNP calling

Requirements
------------
* Scala
* sbt

Dependencies
------------
* Picard library (I have a separate script which fetches it automatically (`download_picard.sh`), but if not, place `sam-<version>.jar` and `picard-<version>.jar` in lib/)

Compiling
---------
In the root directory:

    sbt package

This compiles and puts a biggie.jar file in build/. If you just want to compile without packaging, run:

    sbt compile

Running
-------
    ./bin/biggie cmd args ...

Commands

* all alignment.bam reference.fa
* classify sorted-BAM-file
* snp alignment.bam reference.fa regions.txt

regions.txt:
refSeq	start	end

e.g.
chr21	14800232	14804369
