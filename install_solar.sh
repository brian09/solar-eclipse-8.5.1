#!/bin/bash
SOLAR_RELEASE=8.5.1
CXX_COMPILER=$(which g++)
C_COMPILER=$(which gcc)
FC_COMPILER=$(which gfortran)
SCRIPT=$0
SCRIPT_PATH=$(dirname $SCRIPT)
cd $SCRIPT_PATH
SCRIPT_PATH=$(pwd)
BUILD_DYNAMIC=1
CUDA_PREFIX=0
USE_GPU_TOOLS=0
LIBCUDART_PATH=0
LIBCUBLAS_PATH=0
MAKE_CLEAN=0
args=("$@")
BUILD_THIRD_PARTY_BINARIES=1
SOLAR_SCRIPT_NAME=solar
SOLAR_SCRIPT_PATH=/usr/local/bin
SOLAR_RELEASE_PATH=/opt/appl/solar/$SOLAR_RELEASE
USE_MKL=0
MKL_PATH=/opt/intel/mkl
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
 	echo "Usage: install_solar.sh [Options]"
 	echo "Script Options:"
 	echo "--script_name <Name of script used to startup solar> Default: solar"
 	echo " "
 	echo "--script_path <Folder containg solar startup script> Default: /usr/local/bin"
 	echo "" 
 	echo "--release_path <Folder containing binary and library folders used by solar> Defaut: /opt/appl/solar/$SOLAR_RELEASE"
 	echo "" 
 	echo "--disable_build_third_party_binaries When disabled binaries precompiled for Linux are used Default: Off"
 	echo " "
 	echo "--enable_gpu_tools <CUDA Prefix Directory> <Path to libcudart.$shared_ext> <Path to libcublas.$shared_ext>  Default: Off"
 	echo "Require CUDA Toolkit with minimum version 9.0 and a NVIDA CUDA capable GPU with minimum architecture version 3.5."
 	echo ""
 	echo "--build_static Default: Off"
 	echo "Note: GPU Tools are automatically disabled and cannot be used if operating system is not linux-gnu"
 	echo " " 
 	echo "--cxx_compiler <C++ compiler> Default: $(which g++)"
 	echo "Recommended versions are 6 and 7" 
 	echo ""
  	echo "--c_compiler <C compiler> Default: $(which gcc)"
 	echo "Recommended versions are 6 and 7"
 	echo " "
  	echo "--fc_compiler <Fortran compiler> Default: $(which gfortran)"
 	echo ""
 	echo "--use_mkl <Path to MKL include and lib folders> Default: Disabled"
 	echo "By default Intel Math Kernel Library is installed in /opt/intel which means"
 	echo "that the path would be /opt/intel/mkl. The path must include an include folder"
 	echo "and a lib folder containing: "
 	echo "  intel64/libmkl_core.a"
 	echo "  intel64/libmkl_gnu_thread.a," 
 	echo "  intel64/libmkl_intel_lp64.a" 
 	echo ""
 	echo "--build_fresh Default: Disabled"
 	echo "Builds a completely fresh installation of solar by deleting all previously built libraries and binaries" 
 	echo " "

}


for ((i = 0; i < $#; i++ )) 
do
    
     current_arg=${args[$i]}
    
 if [ "$current_arg" = "--help" ] || [ "$current_arg" = "--help" ]
 then
	help_display
	exit
 elif [ "$current_arg" = "--disable_build_third_party_binaries" ]  
 then
 	BUILD_THIRD_PARTY_BINARIES=0	
  elif [ "$current_arg" = "--build_static" ]  
 then
 	BUILD_DYNAMIC=0	
  elif [ "$current_arg" = "--build_fresh" ]  
 then
 	MAKE_CLEAN=1	 	
 elif [ "$current_arg" = "--enable_gpu_tools" ] && [ $(($i + 3)) -lt $# ] 
 then
 	USE_GPU_TOOLS=1
 	i=$(($i + 1))
 	CUDA_PREFIX=${args[$i]}
 	i=$(($i + 1))
 	LIBCUDART_PATH=${args[$i]}
  	i=$(($i + 1))
 	LIBCUBLAS_PATH=${args[$i]}	
 			
 elif [ "$current_arg" = "--script_name" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	SOLAR_SCRIPT_NAME=${args[$i]} 
  elif [ "$current_arg" = "--script_path" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	SOLAR_SCRIPT_PATH=${args[$i]}	
 	if [ ! -d  $SOLAR_SCRIPT_PATH ]
 	then
 	   mkdir $SOLAR_SCRIPT_PATH
 	fi
 elif [ "$current_arg" = "--release_path" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	SOLAR_RELEASE_PATH=${args[$i]} 
 	 	 	
 elif [ "$current_arg" = "--use_mkl" ] && [ $(($i + 1)) -lt $# ] 
 then
    USE_MKL=1
 	i=$(($i + 1))
 	MKL_PATH=${args[$i]} 	
 	
 elif [ "$current_arg" = "--c_compiler" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	C_COMPILER=${args[$i]} 
 	

 elif [ "$current_arg" = "--cxx_compiler" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	CXX_COMPILER=${args[$i]} 
 elif [ "$current_arg" = "--fc_compiler" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	FC_COMPILER=${args[$i]}
 else
 
 	echo "Invalid argument was entered: $current_arg"
 	echo "See --help for list of valid arguments"
 	exit 1
 fi

done
if [ ! -d  $SCRIPT_PATH/bin_src/solarmain-$SOLAR_RELEASE ]
then 
echo "Build directory for solarmain-$SOLAR_RELEASE not found at $SCRIPT_PATH/bin_src/"
exit 1
fi
mkdir -p $SOLAR_RELEASE_PATH
mkdir -p $SOLAR_RELEASE_PATH/lib
mkdir -p $SOLAR_RELEASE_PATH/bin
OS=$(uname)
if [ "$OS" != "Linux" ] && [ $BUILD_DYNAMIC -eq 0 ] 
then
	BUILD_DYNAMIC=1
	echo "Warning: static build was selected but since operating system is $OSTYPE and not linux-gnu a dynamic build will be used instead"
fi

if [ $BUILD_DYNAMIC -eq 0 ] && [ $USE_GPU_TOOLS -eq 1 ]
then
   echo "Warning: GPU tools will not be included since a static build was selected"
   $USE_GPU_TOOLS=0
fi


if ! command -v $CXX_COMPILER &>/dev/null
then
   echo "C++ compiler not found at $CXX_COMPILER"
   exit 1
fi
if ! command -v $C_COMPILER &>/dev/null
then
   echo "C compiler not found at $C_COMPILER"
   exit 1
fi
if ! command -v $FC_COMPILER &>/dev/null
then
   echo "Fortran compiler not found at $FC_COMPILER"
   exit 1
fi

if [ $USE_GPU_TOOLS -eq 1 ] 
then 

	if [ ! -d $CUDA_PREFIX ]
	then 
    	echo "Unable to located directory $CUDA_PREFIX"
    	exit 1
	fi

	if ! command -v $CUDA_PREFIX/bin/nvcc &>/dev/null
	then
   	echo "nvcc compiler not found at $CUDA_PREFIX/bin/nvcc"
   	exit 1
	fi

	if [ ! -d $CUDA_PREFIX/include ]
	then 
   	 echo "Unable to located directory $CUDA_PREFIX/include"
    	exit 1
	fi

	if ! command -v $LIBCUART_PATH &>/dev/null
	then
   	echo "libcudart.$shared_ext not found at $LIBCUDART_PATH"
   	exit 1
	fi

	if [[ ! -f $LIBCUBLAS_PATH ]]
	then
   	echo "libcublas.$shared_ext not found at $LIBCUBLAS_PATH"
   	exit 1
	fi



 	cd $SCRIPT_PATH/lib_src/libgpu-commands
 	chmod +x ./configure.sh
 	./configure.sh --install_prefix $SOLAR_RELEASE_PATH/lib --cxx_compiler $CXX_COMPILER  --cuda_prefix $CUDA_PREFIX --include_path $SCRIPT_PATH/bin_src/solarmain-$SOLAR_RELEASE/include
	if [  $? -ne 0 ] 
 	then
 	echo "Failed to generate makefile for libgpu-commands.$shared_ext"
 	exit 1
 	fi
 	
	if [ $MAKE_CLEAN -eq 1 ]
 	then
 	make clean
	fi
	
 	make all
 	if [  $? -ne 0 ]
 	then
 	echo "Makefile failed for libgpu-commands.$shared_ext"
 	exit 1
 	fi
 	
 	make install
 	if [  $? -ne 0 ]
 	then
 	echo "Makefile failed for libgpu-commands.$shared_ext"
 	exit 1	
 	fi
fi 


 cd  $SCRIPT_PATH/lib_src/libRicVolume
 chmod +x ./configure.sh 
./configure.sh --install_prefix $SCRIPT_PATH/bin_src/solarmain-$SOLAR_RELEASE --cxx_compiler $CXX_COMPILER --c_compiler $C_COMPILER
if [  $? -ne 0 ]
then
    echo "Failed to configure makefile for libRicVolume.a"
    exit 1
fi
if [ $MAKE_CLEAN -eq 1 ]
then
make clean
fi
make all 
if [  $? -ne 0 ]
then
    echo "Failed to build libRicVolume.a"
    exit 1
fi 
make install
cd  $SCRIPT_PATH/lib_src/libsafe
 chmod +x ./configure.sh
./configure.sh --install_prefix $SCRIPT_PATH/bin_src/solarmain-$SOLAR_RELEASE --c_compiler $C_COMPILER 
if [  $? -ne 0 ]
then
    echo "Failed to configure makefile for libsafe.a"
    exit 1
fi 
if [ $MAKE_CLEAN -eq 1 ]
then
	make clean
fi
make all 
if [  $? -ne 0 ]
then
    echo "Failed to build libsafe.a"
    exit 1
fi 
make install 

cd  $SCRIPT_PATH/lib_src/libtcl8.4
chmod +x ./build_tcl8.4.sh
if [ $MAKE_CLEAN  -eq 1 ]
then
	./build_tcl8.4.sh --install_prefix $SCRIPT_PATH/bin_src/solarmain-$SOLAR_RELEASE --c_compiler $C_COMPILER --build_fresh
	if [  $? -ne 0 ]
	then
		echo "Failed to build and install libtcl8.4.a"
		exit 1
	fi
else
	./build_tcl8.4.sh --install_prefix $SCRIPT_PATH/bin_src/solarmain-$SOLAR_RELEASE --c_compiler $C_COMPILER 
	if [  $? -ne 0 ]
	then
		echo "Failed to build and install libtcl8.4.a"
		exit 1
	fi	
fi	


cd $SCRIPT_PATH/lib_src/libplinkio
chmod +x ./configure.sh
./configure.sh --install_prefix $SCRIPT_PATH/bin_src/solarmain-$SOLAR_RELEASE --c_compiler $C_COMPILER

if [  $? -ne 0 ]
then
	echo "Failed to configure makefile for libplinkio.a"
	exit 1
fi

if [ $MAKE_CLEAN -eq 1 ]
then
	make clean
fi

make all
if [  $? -ne 0 ]
then
    echo "Failed to build libplinkio.a"
    exit 1
fi 
make install

cd $SCRIPT_PATH/bin_src/solarmain-$SOLAR_RELEASE 
if [  $? -ne 0 ]
then
	echo "solarmain-$SOLAR_RELEASE build folder not found at $SCRIPT_PATH/bin_src/"
	exit 1
fi
MKL_ARG=
if [ $USE_MKL -eq 1 ]
then
    MKL_ARG="--use_mkl $MKL_PATH"
fi
chmod +x ./configure.sh
if [ $BUILD_DYNAMIC -eq 0 ]
then
	./configure.sh --install_path $SOLAR_RELEASE_PATH/bin --cxx_compiler $CXX_COMPILER --c_compiler $C_COMPILER --fc_compiler $FC_COMPILER $MKL_ARG --static
	if [  $? -ne 0 ]
	then
		echo "Failed to configure makefile for solarmain"
		exit 1
	fi

elif [ $USE_GPU_TOOLS -eq 1 ]
then
	./configure.sh --install_path $SOLAR_RELEASE_PATH/bin --cxx_compiler $CXX_COMPILER --c_compiler $C_COMPILER --fc_compiler $FC_COMPILER $MKL_ARG --enable_gpu_tools $SOLAR_RELEASE_PATH/lib/libgpu-commands.$shared_ext $LIBCUDART_PATH $LIBCUBLAS_PATH 
	if [  $? -ne 0 ]
	then
		echo "Failed to configure makefile for solarmain"
		exit 1
	fi

else
	./configure.sh --install_path $SOLAR_RELEASE_PATH/bin --cxx_compiler $CXX_COMPILER --c_compiler $C_COMPILER --fc_compiler $FC_COMPILER $MKL_ARG
	if [  $? -ne 0 ]
	then
		echo "Failed to configure makefile for solarmain"
		exit 1
	fi

fi
	
if [ $MAKE_CLEAN -eq 1 ]
then
	make clean
fi
make all 
if [ $? -ne 0 ] 
then
    echo "Failed to build solarmain"
    exit 1
fi

make install

cd $SCRIPT_PATH

if [ $BUILD_THIRD_PARTY_BINARIES -eq 1 ]
then
  cd ./bin_src/third-party-bins
  chmod +x ./build_third_party_bins.sh
  if [ $MAKE_CLEAN -eq 1 ]
  then
  
  	./build_third_party_bins.sh --install_path $SOLAR_RELEASE_PATH/bin  --c_compiler $C_COMPILER --fc_compiler $FC_COMPILER --build_fresh --safelib_path $SCRIPT_PATH/bin_src/solarmain-$SOLAR_RELEASE
   	if [ $? -ne 0 ] 
   	then
   		echo "Failed to install third party binaries"
   	fi
  else
  
  	./build_third_party_bins.sh --install_path $SOLAR_RELEASE_PATH/bin  --c_compiler $C_COMPILER --fc_compiler $FC_COMPILER  --safelib_path $SCRIPT_PATH/bin_src/solarmain-$SOLAR_RELEASE
   	if [ $? -ne 0 ] 
   	then
   		echo "Failed to install third party binaries"
   	fi
  fi
 
else 
    echo "Using precompiled third party binaries that have been compiled for Linux"
    cp $SCRIPT_PATH/build_src/third-party-bins/precompiled_bins/* $SOLAR_RELEASE_PATH/bin
fi


cp -R $SCRIPT_PATH/lib/* $SOLAR_RELEASE_PATH/lib
cp $SCRIPT_PATH/bin/* $SOLAR_RELEASE_PATH/bin

mkdir -p $SOLAR_SCRIPT_PATH
echo "#!/bin/sh" >$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "SOLAR_BIN=$SOLAR_RELEASE_PATH/bin" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "SOLAR_LIB=$SOLAR_RELEASE_PATH/lib" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "PATH=\$SOLAR_BIN:\$PATH" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
if [ $USE_GPU_TOOLS -eq 1 ]
then
	CUDART_PATH=${LIBCUDART_PATH%/libcudart.*}
	CUBLAS_PATH=${LIBCUBLAS_PATH%/libcublas.*}
	if [ $OS = "Linux" ]
	then
	echo "LD_LIBRARY_PATH=\$SOLAR_LIB:\$LD_LIBRARY_PATH:$CUDA_PREFIX/lib64:$CUDART_PATH:$CUBLAS_PATH" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
	else
	echo "DYLD_LIBRARY_PATH=\$SOLAR_LIB:\$DYLD_LIBRARY_PATH:$CUDA_PREFIX/lib64:$CUDART_PATH:$CUBLAS_PATH" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
	fi
else
    if [ $OS = "Linux" ]
    then
	echo "LD_LIBRARY_PATH=\$SOLAR_LIB:\$LD_LIBRARY_PATH" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
	else
	echo "DYLD_LIBRARY_PATH=\$SOLAR_LIB:\$DYLD_LIBRARY_PATH" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
	fi  
fi
echo "TCL_LIBRARY=\$SOLAR_LIB/tcl8.4" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "TK_LIBRARY=\$SOLAR_LIB/tk8.4" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "SOLAR_PROGRAM_NAME=\$0"  >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "export SOLAR_BIN" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "export SOLAR_LIB" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "export PATH" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
if [ $OS = "Linux" ]
then
echo "export LD_LIBRARY_PATH" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
else
echo "export DYLD_LIBRARY_PATH" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
fi
echo "export TCL_LIBRARY" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "export TK_LIBRARY" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "export SOLAR_PROGRAM_NAME" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "if [ ! -s \$SOLAR_BIN/solarmain ]" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "then" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "    echo \"The directory \$SOLAR_BIN for solar binaries is not accessible.\"" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "else" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "    if [ \"\$1\" = \"-noce\"  -o \"\$SOLAR_noce\" != \"\" -o \"\$#\" != \"0\" ]" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "    then" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "        \$SOLAR_BIN/solarmain \$*" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "    else" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "        \$SOLAR_BIN/rlwrap -n \$SOLAR_BIN/solarmain \$*" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "        if [ \"\$?\" != \"0\" ]" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "	      then" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "	          \$SOLAR_BIN/solarmain \$*"  >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "        fi" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "    fi" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "fi" >>$SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
chmod +x $SOLAR_SCRIPT_PATH/$SOLAR_SCRIPT_NAME
echo "Solar Version $SOLAR_RELEASE has been successfully installed"
echo "Startup script $SOLAR_SCRIPT_NAME can be found at $SOLAR_SCRIPT_PATH"

