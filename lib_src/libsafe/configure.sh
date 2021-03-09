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
 	echo " "
 	echo "--install_path  <Prefix of path to installation directories>"
 	echo "Default: /usr/local"
 	echo " "
 	echo "--c_compiler <C compiler> Default: $(which gcc)"
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
if [ ! command -v $C_COMPILER &>/dev/null ]
then
   echo "C compiler not found at $C_COMPILER"
   exit 1
fi 
rm -f sources.mk Makefile
echo "SOURCE_PATH=$SCRIPT_PATH/src" >> sources.mk
echo "INCLUDE_PATH=$SCRIPT_PATH/include" >> sources.mk
echo "C_OBJECTS= \\" >> sources.mk
echo "\$(SOURCE_PATH)/pipeback.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/safelib.o" >> sources.mk
#echo "$(SOURCE_PATH)/w_mehd.c \\" >> sources.mk
echo "HEADERS= \\" >> sources.mk
echo "\$(INCLUDE_PATH)/mvn.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/nrutil.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/pipeback.h \\" >> sources.mk
echo "\$(INCLUDE_PATH)/safelib.h \\" >> sources.mk

echo "include sources.mk" >> Makefile
echo ".SUFFIXES: .c .o .h" >> Makefile

echo ".c.o:" >> Makefile
echo "	$C_COMPILER -c -O3 -w -I\$(INCLUDE_PATH) -o \$@ \$<" >> Makefile
echo ".PHONY: all" >> Makefile
echo "all: libsafe.a" >> Makefile
echo "libsafe.a :  \$(C_OBJECTS) \$(HEADERS)" >> Makefile
if [ $OS = "Linux" ]
then
echo "	ar rcs \$@ \$(C_OBJECTS) " >> Makefile
else
echo "	libtool -static -o \$@ \$(C_OBJECTS) " >> Makefile
fi
echo "install:" >> Makefile
echo "	mkdir -p $INSTALL_PREFIX/include"  >> Makefile
echo "	mkdir -p $INSTALL_PREFIX/lib"  >> Makefile
echo "	cp \$(HEADERS) $INSTALL_PREFIX/include" >> Makefile
echo "	cp libsafe.a $INSTALL_PREFIX/lib" >> Makefile
echo "clean:" >> Makefile
echo "	rm -f \$(C_OBJECTS) libsafe.a" >> Makefile
echo "libsafe makefile successfully made"

