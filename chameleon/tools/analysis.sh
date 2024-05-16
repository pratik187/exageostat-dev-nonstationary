#!/usr/bin/env bash
###
#
#  @file analysis.sh
#  @copyright 2013-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 1.1.0
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 2021-04-20
#
###

# Performs an analysis of Chameleon source code:
# - we consider to be in Chameleon's source code root
# - we consider having the coverage file chameleon_coverage.xml in the root directory
# - we consider having cppcheck, rats, sonar-scanner programs available in the environment

# filter sources:
# - consider generated files in ${BUILDDIR}
# - exclude base *z* files to avoid duplication
# - exclude cblas.h and lapacke-.h because not really part of chameleon and make cppcheck analysis too long

if [ $# -gt 0 ]
then
    BUILDDIR=$1
fi
BUILDDIR=${BUILDDIR:-build}

TOOLSDIR=$(dirname $0)
$TOOLSDIR/find_sources.sh

# Generate coverage xml output
INPUT_FILES=""
for name in $( ls -1 chameleon_*.lcov | grep -v simgrid)
do
    INPUT_FILES="$INPUT_FILES -a $name";
done
lcov $INPUT_FILES -o chameleon.lcov
lcov --summary chameleon.lcov

python3 /usr/local/lib/python3.8/dist-packages/lcov_cobertura.py chameleon.lcov --output chameleon_coverage.xml

# Undefine this because not relevant in our configuration
export UNDEFINITIONS="-UCHAMELEON_USE_OPENCL -UWIN32 -UWIN64 -U_MSC_EXTENSIONS -U_MSC_VER -U__SUNPRO_C -U__SUNPRO_CC -U__sun -Usun -U__cplusplus"

# run cppcheck analysis
CPPCHECK_OPT=" -v -f --language=c --platform=unix64 --enable=all --xml --xml-version=2 --suppress=missingInclude ${UNDEFINITIONS}"
cppcheck $CPPCHECK_OPT --file-list=./filelist_none.txt 2> chameleon_cppcheck.xml
cppcheck $CPPCHECK_OPT -DPRECISION_s -UPRECISION_d -UPRECISION_c -UPRECISION_z --file-list=./filelist_s.txt 2>> chameleon_cppcheck.xml
cppcheck $CPPCHECK_OPT -DPRECISION_d -UPRECISION_s -UPRECISION_c -UPRECISION_z --file-list=./filelist_d.txt 2>> chameleon_cppcheck.xml
cppcheck $CPPCHECK_OPT -DPRECISION_c -UPRECISION_s -UPRECISION_d -UPRECISION_z --file-list=./filelist_c.txt 2>> chameleon_cppcheck.xml
cppcheck $CPPCHECK_OPT -DPRECISION_z -UPRECISION_s -UPRECISION_d -UPRECISION_c --file-list=./filelist_z.txt 2>> chameleon_cppcheck.xml

# Set the default for the project key
SONARQUBE_PROJECTKEY=${SONARQUBE_PROJECTKEY:-hiepacs:chameleon:gitlab:$CI_PROJECT_NAMESPACE:$CI_COMMIT_REF_NAME}

# create the sonarqube config file
cat > sonar-project.properties << EOF
sonar.host.url=https://sonarqube.inria.fr/sonarqube
sonar.login=$SONARQUBE_LOGIN

sonar.links.homepage=$CI_PROJECT_URL
sonar.links.scm=$CI_REPOSITORY_URL
sonar.links.ci=$CI_PROJECT_URL/pipelines
sonar.links.issue=$CI_PROJECT_URL/issues

sonar.projectKey=$SONARQUBE_PROJECTKEY
sonar.projectDescription=Dense linear algebra subroutines for heterogeneous and distributed architectures
sonar.projectVersion=master

sonar.sources=build-openmp/runtime/openmp, build-parsec/runtime/parsec, build-quark/runtime/quark, build-starpu, compute, control, coreblas, example, include, runtime, testing
sonar.inclusions=`cat filelist.txt | sed ':a;N;$!ba;s/\n/, /g'`
sonar.c.includeDirectories=$(echo | gcc -E -Wp,-v - 2>&1 | grep "^ " | tr '\n' ',').,$(find . -type f -name '*.h' | sed -r 's|/[^/]+$||' |sort |uniq | xargs echo | sed -e 's/ /,/g'),$PARSEC_DIR/include,$QUARK_DIR/include,$STARPU_DIR/include/starpu/1.2,$SIMGRID_DIR/include
sonar.sourceEncoding=UTF-8
sonar.c.errorRecoveryEnabled=true
sonar.c.gcc.charset=UTF-8
sonar.c.gcc.regex=(?<file>.*):(?<line>[0-9]+):[0-9]+:\\\x20warning:\\\x20(?<message>.*)\\\x20\\\[(?<id>.*)\\\]
sonar.c.gcc.reportPath=chameleon_build.log
sonar.c.coverage.reportPath=chameleon_coverage.xml
sonar.c.cppcheck.reportPath=chameleon_cppcheck.xml
sonar.c.clangsa.reportPath=build-openmp/analyzer_reports/*/*.plist, build-parsec/analyzer_reports/*/*.plist, build-quark/analyzer_reports/*/*.plist, build-starpu/analyzer_reports/*/*.plist, build-starpu_simgrid/analyzer_reports/*/*.plist
sonar.c.jsonCompilationDatabase=build-openmp/compile_commands.json, build-parsec/compile_commands.json, build-quark/compile_commands.json, build-starpu/compile_commands.json, build-starpu_simgrid/compile_commands.json
sonar.lang.patterns.c++: **/*.cxx,**/*.cpp,**/*.cc,**/*.hxx,**/*.hpp,**/*.hh
sonar.lang.patterns.c: **/*.c,**/*.h
sonar.lang.patterns.python: **/*.py
EOF

# run sonar analysis + publish on sonarqube-dev
sonar-scanner -X > sonar.log
