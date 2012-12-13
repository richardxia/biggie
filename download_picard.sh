#!/bin/sh

VERSION=1.81
FOLDER=picard-tools-$VERSION
wget http://downloads.sourceforge.net/project/picard/picard-tools/$VERSION/$FOLDER.zip
unzip -j $FOLDER.zip $FOLDER/sam-$VERSION.jar $FOLDER/picard-$VERSION.jar -d lib/
rm $FOLDER.zip
