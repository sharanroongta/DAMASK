#!/bin/bash

cat  README
echo

if [ "$OSTYPE" == "linux-gnu" ] || [ "$OSTYPE" == 'linux' ]; then
  DAMASK_ROOT=$(readlink -f "`dirname $BASH_SOURCE`")
else
  [[ "${BASH_SOURCE::1}" == "/" ]] && BASE="" || BASE="`pwd`/"
  STAT=$(stat "`dirname $BASE$BASH_SOURCE`")
  DAMASK_ROOT=${STAT##* }
fi

BUILDROOT=$DAMASK_ROOT/build
BUILDDIR=spectral

# prepare building directory
# structure:
#   BUILD_DIR
#   |-BUILD_SPECTRAL
#   |-BUILD_FEM
#   |-BUILD_MARC
if [ ! -d $BUILDROOT ]; then
    mkdir $BUILDROOT
fi
cd $BUILDROOT
if [ -d $BUILDDIR ] ; then
    rm -rf $BUILDDIR
fi
mkdir $BUILDDIR
cd $BUILDDIR

##
# CMake call
# PETSC_DIR                |  PETSC directory
# CMAKE_BUILD_TYPE         |  Default set to release (no debugging output)
# CMAKE_VERBOSE_MAKEFILE   |  [ON/OFF] toggle makefile verbose output
# OPENMP                   |  [ON/OFF]
# OPTIMIZATION             |  [OFF,DEFENSIVE,AGGRESSIVE,ULTRA]
# DAMASK_DRIVER            |  [SPECTRAL, FEM]
# DAMASK_INSTALL           |  Directory to install binary output
# BUILDCMD_PRE             |  Compiler prefix
# BUILDCMD_POST            |  Compiler suffix,
#                          |    e.g. DAMASK_SUFFIX="-Wl,--verbose"
cmake -D PETSC_DIR=${PETSC_DIR}          \
      -D CMAKE_BUILD_TYPE=RELEASE        \
      -D CMAKE_VERBOSE_MAKEFILE=OFF      \
      -D OPENMP=ON                       \
      -D OPTIMIZATION=DEFENSIVE          \
      -D DAMASK_DRIVER=SPECTRAL          \
      -D DAMASK_INSTALL=${HOME}/bin      \
      -D BUILDCMD_PRE=""                 \
      -D BUILDCMD_POST="-Wl,--verbose"   \
      ../..

echo "Start compiling DAMASK_spectral"
make
make install