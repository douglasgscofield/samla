#!/bin/bash
revisioncount=`git log --oneline | wc -l | sed -e 's/^\s\+//g' -e 's/\s\+$//g'`
projectversion=`git describe --tags --long`
cleanversion=${projectversion%%-*}
reference=`git config --get remote.origin.url`

echo "#define VERSION \"$projectversion-$revisioncount\""
#echo "#define VERSION \"$cleanversion.$revisioncount\""

echo "#define REFERENCE \"$reference\""

