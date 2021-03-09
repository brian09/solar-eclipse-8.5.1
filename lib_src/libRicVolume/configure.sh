#!/bin/bash
SCRIPT=$0
SCRIPT_PATH=$(dirname $SCRIPT)
cd $SCRIPT_PATH
SCRIPT_PATH=$(pwd)
args=("$@")
INSTALL_PREFIX=/usr/local
CXX_COMPILER=$(which g++)
C_COMPILER=$(which gcc)
SOURCE_PATH=$SCRIPT_PATH/src
shared_ext="so"
OS="Linux"
if [[ "$OSTYPE" == "linux-gnu"* ]]
then
        OS="Linux"
        shared_ext="so"
elif [[ "$OSTYPE" == "darwin"* ]]
then
        OS="Mac"
        shared_ext="dylib"

else
        echo "Not setup for $OSTYPE build"
        echo "Only linux-gnu and darwin are accepted operating system"
        exit 1
fi

function help_display {
 	echo "Usage: configure.sh [Options]"
 	echo "Script Options:"
 	echo " "
 	echo "--install_prefix <Prefix of path to installation directories>"
 	echo "Default: /usr/local"
 	echo " "
 	echo "--cxx_compiler <C++ compiler> Default: $(which g++)"
 	echo " " 
 	echo "-c_compiler <C compiler> Default: $(which gcc)"
 	echo ""

}


for ((i = 0; i < $#; i++ )); do
	current_arg=${args[$i]}

 if [ "$current_arg" = "--help" ] || [ "$current_arg" = "--help" ]
 then
	help_display
	exit
 elif [ "$current_arg" = "--install_prefix" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	INSTALL_PREFIX=${args[$i]}
 elif [ "$current_arg" = "--cxx_compiler" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	CXX_COMPILER=${args[$i]} 
 
 elif [ "$current_arg" = "--c_compiler" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	C_COMPILER=${args[$i]}	
 
 else
 
 	echo "Invalid argument was entered: $current_arg"
 	echo "See --help for list of valid arguments"
 	exit 1
 fi

done
if [ ! command -v $CXX_COMPILER &>/dev/null ]
then
   echo "C++ compiler not found at $CXX_COMPILER"
   exit 1
fi
if [ ! command -v $C_COMPILER &>/dev/null ]
then
   echo "C compiler not found at $C_COMPILER"
   exit 1
fi


rm -f Makefile
rm -f sources.mk
echo "SOURCE_PATH=$SCRIPT_PATH/src" >> sources.mk
echo "CXX_OBJECTS= \\" >> sources.mk
echo "\$(SOURCE_PATH)/nemardr.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/NiftiOrientConvert.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/RicUtil.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/RicVolume.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/RicVolumeSet.o" >> sources.mk
echo "C_OBJECTS= \\" >> sources.mk
echo "\$(SOURCE_PATH)/adler32.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/compress.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/crc32.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/deflate.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/gzio.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/infback.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/inffast.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/inflate.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/inftrees.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/nifti1_io.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/trees.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/uncompr.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/znzlib.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/zutil.o" >> sources.mk
echo "HEADERS= \\" >> sources.mk
echo "\$(SOURCE_PATH)/crc32.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/deflate.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/dtypes.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/inffast.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/inffixed.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/inflate.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/inftrees.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/nemardr.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/nifti1.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/nifti1_io.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/RicMatrix.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/RicPoint.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/RicUtil.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/RicVolume.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/RicVolumeSet.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/trees.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/zconf.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/zconf.in.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/zlib.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/znzlib.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/zutil.h " >> sources.mk
echo "include sources.mk" >> Makefile
echo ".SUFFIXES: .c .cpp .o .h" >> Makefile
echo ".cpp.o:" >> Makefile
echo "	$CXX_COMPILER -c -O3 -std=c++11 -I\$(SOURCE_PATH) -w -o \$@ \$<" >> Makefile
echo ".c.o:" >> Makefile
echo "	$C_COMPILER -c -O3 -I\$(SOURCE_PATH) -w -o \$@ \$<" >> Makefile

echo ".PHONY: all" >> Makefile
echo "all: libRicVolume.a" >> Makefile
echo "libRicVolume.a:  \$(C_OBJECTS) \$(CXX_OBJECTS) \$(HEADERS)" >> Makefile
if [ $OS = "Linux" ]
then
echo "	ar rcs \$@ \$(C_OBJECTS) \$(CXX_OBJECTS)" >> Makefile
else
echo "	libtool -static -o \$@ \$(C_OBJECTS) \$(CXX_OBJECTS)" >> Makefile
fi
echo "install:" >> Makefile
echo "	mkdir -p $INSTALL_PREFIX/include"  >> Makefile
echo "	mkdir -p $INSTALL_PREFIX/lib"  >> Makefile
echo "	cp \$(HEADERS) $INSTALL_PREFIX/include" >> Makefile
echo "	cp libRicVolume.a $INSTALL_PREFIX/lib" >> Makefile
echo "clean:" >> Makefile
echo "	rm -f \$(C_OBJECTS) \$(CXX_OBJECTS) libRicVolume.a" >> Makefile
echo "libRicVolume makefile successfully made"
