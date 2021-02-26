SCRIPT=$(realpath $0)
SCRIPT_PATH=$(dirname $SCRIPT)
args=("$@")
INSTALL_PATH=/usr/local/bin
C_COMPILER=$(which gcc)
FC_COMPILER=$(which gfortran)
MAKE_CLEAN=0
SAFELIB_PATH=/usr/local
shared_ext="so"
OS="Linux"
if [[ "$OSTYPE" = "linux-gnu" ]]
then
        OS="Linux"
        shared_ext="so"
elif [[ "$OSTYPE" = "darwin" ]]
then
        OS="Mac"
        shared_ext="dylib"

else
        echo "Not setup for $OSTYPE build"
        echo "Only linux-gnu and darwin are accepted operating system"
        exit 1
fi
function help_display {
 	echo "Usage: build_thread_party_bins.sh [Options]"
 	echo "Script Options:"
 	echo ""
 	echo "--install_path <Path to solar-eclipse bin folder> Default: /usr/local/bin"
 	echo ""
 	echo "--c_compiler <C compiler> Default: $(which gcc)"
 	echo " "
 	echo "--fc_compiler <Fortran compiler> Default: $(which gfortran)"
 	echo ""
 	echo "--safelib_path <Path to libsafe.a and it's include files> Default: /usr/local "
 	echo "Path requires both include folder and lib folder" 
 	echo "Default locations are /usr/local/lib and /usr/local/include"
 	echo ""
 	echo "--build_fresh Default: Disabled"
 	echo "Builds completely fresh executable binaries by removing previous compiled and linked code" 
 	echo ""
}
for ((i = 0; i < $#; i++ )); do
	 current_arg=${args[$i]}
 if [ "$current_arg" = "--help" ] || [ "$current_arg" = "--help" ]
 then
	help_display
	exit
 elif [ "$current_arg" = "--install_path" ]  && [ $(($i + 1)) -lt $# ]
 then
 	i=$(($i + 1))
 	INSTALL_PATH=${args[$i]}
 elif [ "$current_arg" = "--c_compiler" ]  && [ $(($i + 1)) -lt $# ]
 then
 	i=$(($i + 1))
 	C_COMPILER=${args[$i]}
 elif [ "$current_arg" = "--fc_compiler" ]  && [ $(($i + 1)) -lt $# ]
 then
 	i=$(($i + 1))
 	FC_COMPILER=${args[$i]}	
 elif [ "$current_arg" = "--safelib_path" ]  && [ $(($i + 1)) -lt $# ]
 then
 	i=$(($i + 1))
 	SAFELIB_PATH=${args[$i]}	 			 
 elif [ "$current_arg" = "--build_fresh" ]  && [ $(($i + 1)) -lt $# ]
 then
 	MAKE_CLEAN = 1	
 else
 
 	echo "Invalid argument was entered: $current_arg"
 	echo "See --help for list of valid arguments"
 	exit 1
 fi
 done 
if [ ! command -v $FC_COMPILER &>/dev/null ]
then
   echo "Fortran compiler not found at $FC_COMPILER"
   exit 1
fi
if [ ! command -v $C_COMPILER &>/dev/null ]
then
   echo "C compiler not found at $C_COMPILER"
   exit 1
fi
if [ ! -f $SAFELIB_PATH/lib/libsafe.a ]; then
	echo "libsafe.a not found at $SAFELIB_PATH/lib"
	exit 1
fi
if [ ! -f $SAFELIB_PATH/include/safelib.h ]; then
	echo "safelib.h not found at $SAFELIB_PATH/include"
	exit 1
fi
if [ ! -f $SAFELIB_PATH/include/pipeback.h ]; then
	echo "pipeback.h not found at $SAFELIB_PATH/include"
	exit 1
fi
if [ ! -f $SAFELIB_PATH/include/nrutil.h ]; then
	echo "nrutil.h not found at $SAFELIB_PATH/include"
	exit 1
fi
if [ ! -f $SAFELIB_PATH/include/mvn.h ]; then
	echo "mvn.h not found at $SAFELIB_PATH/include"
	exit 1
fi
mkdir -p $INSTALL_PATH
CC=$C_COMPILER 
export CC  
cd $SCRIPT_PATH/fastlink-src
if [ $MAKE_CLEAN -eq 1 ]
then 
	make clean
	if [  $? -ne 0 ]
	then
		echo "Failed to run make clean on third party binary mlink"
		exit 1
	fi
fi
make all 
if [  $? -ne 0 ] 
then
  echo "Failed to build third party binary mlink"
  exit 1
fi

cp mlink $INSTALL_PATH
echo "Successfully installed mlink"
cd $SCRIPT_PATH/homo
if [ $MAKE_CLEAN -eq 1 ]
then 
	make clean
	if [  $? -ne 0 ]
	then
		echo "Failed to run make clean on third party binary homo"
		exit 1
	fi
fi
make 
if [  $? -ne 0 ]
then
  echo "Failed to build third party binary homo"
  exit 1
fi
cp homo $INSTALL_PATH
echo "Successfully installed homo"
FC=$FC_COMPILER
export FC
export SAFELIB_PATH
cd $SCRIPT_PATH/ibd-src
if [ $MAKE_CLEAN -eq 1 ]
then 
	make clean
	if [  $? -ne 0 ]
	then	           
		echo "Failed to run make clean on third party binaries allfreq, genfreq, ibdmat, ibdprep, multipnt, relate, unknown, getmeans, ibdmc, makeped, and mrgibd"
		exit 1
	fi	
fi
make all
if [  $? -ne 0 ]
then
  echo "Failed to build third party binaries allfreq, genfreq, ibdmat, ibdprep, multipnt, relate, unknown, getmeans, ibdmc, makeped, and mrgibd"
  exit 1
fi
cp allfreq $INSTALL_PATH
cp genfreq $INSTALL_PATH
cp ibdmat $INSTALL_PATH
cp ibdprep $INSTALL_PATH
cp multipnt $INSTALL_PATH
cp relate $INSTALL_PATH
cp unknown $INSTALL_PATH
cp getmeans $INSTALL_PATH 
cp ibdmc $INSTALL_PATH
cp makeped $INSTALL_PATH
cp mrgibd $INSTALL_PATH
cd $SCRIPT_PATH/third-party-scripts
cp domcibd $INSTALL_PATH
cp open-x11a $INSTALL_PATH
cp snpplotk $INSTALL_PATH
cp stringplotk $INSTALL_PATH

echo "Successfully install allfreq, domcibd, genfreq, ibdmat, ibdprep, multipnt, relate, unknown, getmeans, ibdmc, makeped, and mrgibd"


cd $SCRIPT_PATH/rlwrap-master
chmod +x ./configure
./configure
if [  $? -ne 0 ]
then 
 echo "Configure script for rlwrap failed"
 echo "Check for needed libraries"
 exit 1
fi
if [ $MAKE_CLEAN -eq 1 ]
then 
	make clean
	if [  $? -ne 0 ]
	then	           
		echo "Failed to run make clean on rlwrap"
		exit 1
	fi	
fi
make all 
if [  $? -ne 0 ]
then 
 echo "Failed to build rlwrap"
 exit 1
fi
cp ./src/rlwrap $INSTALL_PATH
echo "Third party binaries have been successfully installed"
