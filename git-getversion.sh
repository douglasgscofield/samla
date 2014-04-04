#!/bin/bash
revisioncount=`git log --oneline | wc -l | sed -e 's/^\s\+//g' -e 's/\s\+$//g'`
projectversion=`git describe --tags --long`
cleanversion=${projectversion%%-*}

echo "#define VERSION \"$projectversion-$revisioncount\""
#echo "#define VERSION \"$cleanversion.$revisioncount\""
