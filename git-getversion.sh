#!/bin/bash
revisioncount=`git log --oneline | wc -l | sed -e 's/^\s\+//g' -e 's/\s\+$//g'`
projectversion=`git describe --tags --long`
cleanversion=${projectversion%%-*}
reference=`git config --get remote.origin.url`
vcfrevisioncount=`(cd vcflib && git log --oneline | wc -l | sed -e 's/^\s\+//g' -e 's/\s\+$//g')`
vcfversion=`(cd vcflib && git describe --all --long)`
vcfreference=`(cd vcflib && git config --get remote.origin.url)`

echo "#define VERSION \"$projectversion-$revisioncount\""
#echo "#define VERSION \"$cleanversion.$revisioncount\""

echo "#define REFERENCE \"$reference\""

echo "#define VCFLIB_VERSION \"$vcfversion-$vcfrevisioncount\""

echo "#define VCFLIB_REFERENCE \"$vcfreference\""

# for the below to work, the makefile line should define CXX and CXXFLAGS:
# CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)" ./git-getversion.sh > version.h
#
# cxxversion=`${CXX} --version | head -n 1`
# cxxflags="${CXXFLAGS}"
# echo "#define CXX_VERSION \"$cxxversion\""
# echo "#define CXXFLAGS \"$cxxflags\""
