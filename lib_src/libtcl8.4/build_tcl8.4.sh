#!/bin/bash
SCRIPT=$0
SCRIPT_PATH=$(dirname $SCRIPT)
cd $SCRIPT_PATH
SCRIPT_PATH=$(pwd)
INSTALL_PATH=/usr/local
args=("$@")
MAKE_CLEAN=0
C_COMPILER=$(which gcc)
function help_display {
 	echo "Usage: build_tcl8.4.sh [Options]"
 	echo "" 
 	echo "--install_prefix <Prefix of path to lib folder where libsafe.a and include folder where header files are installed>"
 	echo "Default: /usr/local"
 	echo " "
 	echo "--build_fresh Default: Disabled"
 	echo "Builds a completely fresh version libtcl8.4.a by removing previous compiled and linked code"
 	echo ""
 	echo "--c_compiler <Path to c compiler> Default: $(which gcc)"
 	echo " "
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
 elif [ "$current_arg" = "--build_fresh" ] 
 then
	MAKE_CLEAN=1
	
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

cd  unix
if [ $MAKE_CLEAN -eq 1 ]
then
   make clean
   
fi

CC=$C_COMPILER
export CC
chmod +x $SCRIPT_PATH/unix/configure 
$SCRIPT_PATH/unix/configure --enable-shared=no --enable-threads --enable-64bit --prefix=$INSTALL_PREFIX
if [ $? -ne 0 ]
then
	echo "Failed to configure tcl8.4"
	exit 1
fi
if [ $MAKE_CLEAN -eq 1 ]
then
   make clean
fi

make -C $SCRIPT_PATH/unix all 
if [  $? -ne 0 ]
then
	echo "Failed to build tcl8.4"
	exit 1
fi
make -C $SCRIPT_PATH/unix install
if [ $? -ne 0 ]
then
	echo "Failed to install tcl8.4"
	exit 1
fi
