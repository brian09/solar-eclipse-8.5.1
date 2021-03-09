#!/bin/bash

SCRIPT=$0
SCRIPT_PATH=$(dirname $SCRIPT)
cd $SCRIPT_PATH
SCRIPT_PATH=$(pwd)
args=("$@")
INSTALL_PREFIX=/usr/local
C_COMPILER=$(which gcc)
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
 	echo ""
 	echo "--install_prefix <Prefix of path to installation directories> Default: /usr/local" 
 	echo "libplinkio.a will be installed at PREFIX/lib and header files installed at PREFIX/include"
 	echo " "
 	echo "--c_compiler <C compiler> Default: $(which gcc)"
 	echo ""
}


for ((i = 0; i < $#; i++ ))
 do
	current_arg=${args[$i]}

 if [ "$current_arg" = "--help" ] || [ "$current_arg" = "--help" ]
 then
	help_display
	exit
  elif [ "$current_arg" = "--c_compiler" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	C_COMPILER=${args[$i]}
 			
 elif [ "$current_arg" = "--install_prefix" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	INSTALL_PREFIX=${args[$i]} 
 else
 	
 	echo "Invalid argument was entered: $current_arg"
 	echo "See --help for list of valid arguments"
 	exit 1
 fi

done
if [ ! command -v $C_COMPILER &>/dev/null ]
then
   echo "C compiler not found at $C_COMPILER"
   exit 1
fi
mkdir -p $INSTALL_PREFIX

rm -f Makefile
rm -f sources.mk
echo "SOURCE_PATH=$SCRIPT_PATH/src" >> sources.mk
echo "INCLUDE_PATH=$SCRIPT_PATH/include" >> sources.mk
echo "C_OBJECTS= \\" >> sources.mk
echo "\$(SOURCE_PATH)/plinkio.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/libcsv.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/file.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/fam_parse.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/fam.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/bim_parse.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/bim.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/bed_header.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/bed.o" >> sources.mk
echo "HEADERS= \\" >> sources.mk
echo "\$(INCLUDE_PATH)/csv.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/plinkio.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/plinkio/bed.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/plinkio/bed_header.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/plinkio/bim.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/plinkio/bim_parse.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/plinkio/fam.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/plinkio/fam_parse.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/plinkio/file.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/plinkio/snp_lookup.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/plinkio/snp_lookup_big.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/plinkio/snp_lookup_little.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/plinkio/status.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/plinkio/utarray.h" >> sources.mk


echo "include sources.mk" >> Makefile
echo ".SUFFIXES: .c .o .h" >> Makefile
echo ".c.o:" >> Makefile
echo "	$C_COMPILER -c -O3 -I\$(INCLUDE_PATH) -o \$@ \$<" >> Makefile

echo ".PHONY: all" >> Makefile
echo "all: libplinkio.a" >> Makefile
echo "libplinkio.a :  \$(C_OBJECTS) \$(HEADERS)" >> Makefile

if [ $OS = "Linux" ]
then

echo "	ar rcs \$@ \$(C_OBJECTS) " >> Makefile

else

echo "	libtool -static -o \$@ \$(C_OBJECTS) " >> Makefile

fi

echo "install:" >> Makefile	
echo "	mkdir -p $INSTALL_PREFIX/include"  >> Makefile
echo "	mkdir -p $INSTALL_PREFIX/lib"  >> Makefile
echo "	cp -R $SCRIPT_PATH/include/* $INSTALL_PREFIX/include" >> Makefile
echo "	cp libplinkio.a $INSTALL_PREFIX/lib" >> Makefile
echo "clean:" >> Makefile
echo "	rm -f \$(C_OBJECTS) libplinkio.a" >> Makefile
echo "libplinkio makefile successfully made"
