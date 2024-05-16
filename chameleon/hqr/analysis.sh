#!/usr/bin/env bash
###
#
#  @file analysis.sh
#  @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.0.0
#  @author Mathieu Faverge
#  @date 2021-03-25
#
###

# Performs an analysis of HQR source code:
# - we consider to be in HQR's source code root
# - we consider having the coverage file hqr.lcov in the root directory
# - we consider having cppcheck, rats, sonar-scanner programs available in the environment

if [ $# -gt 0 ]
then
    BUILDDIR=$1
fi
BUILDDIR=${BUILDDIR:-build}

# List source files:
rm -f filelist.txt
git ls-files | grep "\.[ch]$" > filelist.txt

# Generate coverage analysis report
python3 /usr/local/lib/python3.8/dist-packages/lcov_cobertura.py hqr.lcov --output hqr-coverage.xml

# to get it displayed and captured by gitlab to expose the badge on the main page
cat ./hqr-gcov.log

# Undefine this because not relevant in our configuration
export UNDEFINITIONS="-UWIN32 -UWIN64 -U_MSC_EXTENSIONS -U_MSC_VER -U__SUNPRO_C -U__SUNPRO_CC -U__sun -Usun -U__cplusplus"

# run cppcheck analysis
cppcheck -v -f --project=build/compile_commands.json --language=c --platform=unix64 --enable=all --xml --xml-version=2 --suppress=missingInclude ${UNDEFINITIONS} ${DEFINITIONS} --file-list=./filelist.txt 2> hqr-cppcheck.xml

# run rats analysis
rats -w 3 --xml  `cat filelist.txt` > hqr-rats.xml

# Set the default for the project key
SONARQUBE_PROJECTKEY=${SONARQUBE_PROJECTKEY:-hiepacs:hqr:gitlab:master}

# create the sonarqube config file
cat > sonar-project.properties << EOF
sonar.host.url=https://sonarqube.inria.fr/sonarqube
sonar.login=$SONARQUBE_LOGIN

sonar.links.homepage=$CI_PROJECT_URL
sonar.links.scm=$CI_REPOSITORY_URL
sonar.links.ci=$CI_PROJECT_URL/pipelines
sonar.links.issue=$CI_PROJECT_URL/issues

sonar.projectKey=$SONARQUBE_PROJECTKEY
sonar.projectDescription=Library for hierarchical QR/LQ reduction trees
sonar.projectVersion=master

sonar.sources=include, src, testings
sonar.inclusions=`cat filelist.txt | xargs echo | sed 's/ /, /g'`
sonar.sourceEncoding=UTF-8
sonar.c.errorRecoveryEnabled=true
sonar.c.compiler.charset=UTF-8
sonar.c.compiler.parser=GCC
sonar.c.compiler.regex=^(.*):(\\d+):\\d+: warning: (.*)\\[(.*)\\]$
sonar.c.compiler.reportPath=hqr-build.log
sonar.c.coverage.reportPath=hqr-coverage.xml
sonar.c.cppcheck.reportPath=hqr-cppcheck.xml
sonar.c.rats.reportPath=hqr-rats.xml
sonar.c.jsonCompilationDatabase=${BUILDDIR}/compile_commands.json
sonar.lang.patterns.c++: **/*.cxx,**/*.cpp,**/*.cc,**/*.hxx,**/*.hpp,**/*.hh
sonar.lang.patterns.c: **/*.c,**/*.h
sonar.lang.patterns.python: **/*.py
EOF

# run sonar analysis + publish on sonarqube-dev
sonar-scanner -X > sonar.log
