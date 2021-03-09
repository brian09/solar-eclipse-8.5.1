#!/bin/bash

SCRIPT=$0
SCRIPT_PATH=$(dirname $SCRIPT)
cd $SCRIPT_PATH
SCRIPT_PATH=$(pwd)
args=("$@")
CUDA_PREFIX=/usr/local/cuda
INSTALL_PREFIX=/usr/local/lib
INCLUDE_PATH=$SCRIPT_PATH/include
CXX_COMPILER=$(which g++)
OS="Linux"
shared_ext="so"
if [[ "$OSTYPE" = "linux-gnu"* ]]
then
        OS="Linux"
        shared_ext="so"
elif [[ "$OSTYPE" = "darwin"* ]]
then
        OS="Mac"
        shared_ext="dylib"
else
        echo "Not setup for $OSTYPE build"
        echo "Only linux-gnu and darwin are accepted operating system"
        exit 1
fi

function help_display {
 	echo "Usage: configure [Options]"
 	echo "Script Options:"
 	echo ""
 	echo "Require CUDA Toolkit with minimum version 9.0 and a NVIDA CUDA capable GPU with minimum architecture version 3.5"
 	echo ""
 	echo "--cuda_prefix <Path to CUDA distribution containing include and bin folders> Default: /usr/local/cuda"
 	echo "The include folder must contain CUDA header files and the bin must contain compiler nvcc"
 	echo "Header files will be read from CUDA_PREFIX/include and nvcc read from CUDA_PREFIX/bin/nvcc"
 	echo ""
 	echo "--install_prefix <Path to folder where libgpu-commands.$shared_ext is installed>" 
 	echo "Default: /usr/local/lib"
 	echo " "
 	echo "--include_path <Path to directory containing solar-eclipse header files>"
 	echo "Default: $SCRIPT_PATH/include"
 	echo "This folder should include the following headers template library Eigen, RicVolume library headers" 
 	echo "(e.g. RicVolume.h, nifti1.h, and RicMatrix.h), safelib.h, expression.h, libplinkio headers, libtcl8.4 headers," 
 	echo "solar.h, and tablefile.h" 
 	echo "By default $SCRIPT_PATH/include should include all these files and there is no need"
 	echo "to specify this folder"
 	echo " "
 	echo "--cxx_compiler <C++ compiler> Default: $(which g++)"
 	echo " " 
}


for ((i = 0; i < $#; i++ ))
 do
	current_arg=${args[$i]}

 if [ "$current_arg" = "--help" ] || [ "$current_arg" = "--help" ]
 then
	help_display
	exit
 elif [ "$current_arg" = "--cuda_prefix" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	CUDA_PREFIX=${args[$i]}
  elif [ "$current_arg" = "--cxx_compiler" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	CXX_COMPILER=${args[$i]}
 			
 elif [ "$current_arg" = "--install_prefix" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	INSTALL_PREFIX=${args[$i]} 
 elif [ "$current_arg" = "--include_path" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	INCLUDE_PATH=${args[$i]} 

 else
 	
 	echo "Invalid argument was entered: $current_arg"
 	echo "See --help for list of valid arguments"
 	exit 1
 fi

done
if  [ ! -d $INSTALL_PREFIX ]
then
   echo "Install path: $INSTALL_PREFIX not found"
   exit 1
fi
if [ ! command -v $CXX_COMPILER &>/dev/null ]
then
   echo "C++ compiler not found at $CXX_COMPILER"
   exit 1
fi

if  [ ! -d $CUDA_PREFIX/include ] || [  ! command -v $CUDA_PREFIX/bin/nvcc &>/dev/null ]
then
   echo "Folder containing CUDA include folder and nvcc compiler not found using prefix $CUDA_PREFIX"
   echo "Folder requires bin folder containing nvcc and include folder containing CUDA headers"
   exit 1
fi
TOOLKIT_VERSION_MAJOR=$($CUDA_PREFIX/bin/nvcc --version | grep "release" | awk '{print $6}' | cut -c2- | cut  -d . -f 1)
TOOLKIT_VERSION_MINOR=$($CUDA_PREFIX/bin/nvcc --version | grep "release" | awk '{print $6}' | cut -c2- | cut  -d . -f 2)

if [ $TOOLKIT_VERSION_MAJOR -lt 9 ]
then
	echo "CUDA Toolkit version cannot be less than 9.0"
	exit 1
fi


if [ $TOOLKIT_VERSION_MAJOR -eq 9 ] 
then
	GENCODE_FLAGS='-gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_37,code=sm_37 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_70,code=compute_70'

elif  [ $TOOLKIT_VERSION_MAJOR -eq 10 ] 
then 
	GENCODE_FLAGS='-gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_37,code=sm_37 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_75,code=compute_75'
	
	
elif [ $TOOLKIT_VERSION_MAJOR -eq 11 ] && [ $TOOLKIT_VERSION_MINOR -eq 0 ]
then
	GENCODE_FLAGS='-gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_37,code=sm_37 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=sm_61  -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80'

elif [ $TOOLKIT_VERSION_MAJOR -eq 11 ] && [ $TOOLKIT_VERSION_MINOR -eq 1 ]
then
	GENCODE_FLAGS='-gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_37,code=sm_37 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_86,code=sm_86 -gencode=arch=compute_86,code=compute_86'

else
	GENCODE_FLAGS='-gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_37,code=sm_37 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_86,code=sm_86'				
fi


if  [ ! -d $INCLUDE_PATH ]
then
   echo "Folder $INCLUDE_PATH not found"
   exit 1
fi
if [ ! -d $INCLUDE_PATH/Eigen ]
then 
    echo "Header file folder Eigen not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -d $INCLUDE_PATH/plinkio ]
then 
    echo "Header file folder plinkio not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/expression.h ]
then 
    echo "Header file expression.h not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/nifti1.h ]
then 
    echo "Header file nifti1.h not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/nifti1_io.h ]
then 
    echo "Header file nifti1_io.h not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/plinkio.h ]
then 
    echo "Header file plinkio.h not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/RicMatrix.h ]
then 
    echo "Header file RicMatrix.h not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/RicPoint.h ]
then 
    echo "Header file RicPoint.h not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/RicUtil.h ]
then 
    echo "Header file RicUtil.h not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/RicVolume.h ]
then 
    echo "Header file RicVolume.h not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/RicVolumeSet.h ]
then 
    echo "Header file RicVolumeSet.h not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/safelib.h ]
then 
    echo "Header file safelib.h not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/solar.h ]
then 
    echo "Header file solar.h not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/tablefile.h ]
then 
    echo "Header file tablefile.h not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/tcl.h ]
then 
    echo "Header file tcl.h not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/tclDecls.h ]
then 
    echo "Header file tclDecls.h not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/tclPlatDecls.h ]
then 
    echo "Header file tclPlatDecls.h not found in $INCLUDE_PATH"
    exit 1
fi
if [ ! -f $INCLUDE_PATH/znzlib.h ]
then 
    echo "Header file znzlib.h not found in $INCLUDE_PATH"
    exit 1
fi

rm -f Makefile
rm -f sources.mk
echo "INCLUDES= \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-data.h \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-exception.h \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-fphi-settings.h \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-fphi-variables.h \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-gwas-estimator-context.h \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-gwas-estimator.h \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-gwas.h \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-gwas-screen-variables.h \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-gwas-settings.h \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-pedigree-data.h \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-selection.h \\" >> sources.mk
echo "$SCRIPT_PATH/src/solar-trait-reader.h \\" >> sources.mk
echo " " >> sources.mk
echo "OBJECTS= \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-fphi.o \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-fphi-command.o \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-fphi-variables.o \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-gwas.o \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-gwas-command.o \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-gwas-estimator.o \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-gwas-estimator-context.o \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-gwas-screen.o \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-gwas-screen-variables.o \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-pedifromsnps.o \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-pedigree-data.o \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-selection.o \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-pedifromsnps-kernels.o \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-fphi-kernels.o \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-gwas-kernels.o \\" >> sources.mk
echo "$SCRIPT_PATH/src/gpu-gwas-screen-kernels.o \\" >> sources.mk
echo "include sources.mk" >> Makefile
echo ".SUFFIXES: .cu .cc .o .h" >> Makefile
echo ".cu.o:" >> Makefile
echo "	$CUDA_PREFIX/bin/nvcc -c   -maxrregcount 62  -std=c++11 -O3 -use_fast_math -Xcompiler -fPIC -Xcompiler -fexceptions -Xcompiler -w $GENCODE_FLAGS -Xcompiler -I$CUDA_PREFIX/include -Xcompiler -I$INCLUDE_PATH -o \$@ \$<" >> Makefile
echo ".cc.o:" >> Makefile
echo "	$CXX_COMPILER -c  -I$INCLUDE_PATH -pthread -fopenmp -O3 -std=c++11 -fPIC -w -I$CUDA_PREFIX/include -o \$@ \$<" >> Makefile
echo ".PHONY: all" >> Makefile
if [ $OS = "Linux" ]
then

echo "all: libgpu-commands.so" >> Makefile
echo "libgpu-commands.so :  \$(OBJECTS) \$(INCLUDES)" >> Makefile
echo "	$CUDA_PREFIX/bin/nvcc $GENCODE_FLAGS -Xlinker -fexceptions -Xlinker -O3  -shared -o \$@ \$(OBJECTS)" >> Makefile
echo "install:" >> Makefile	
echo "	mkdir -p $INSTALL_PREFIX"  >> Makefile
echo "	cp libgpu-commands.so $INSTALL_PREFIX" >> Makefile
echo "clean:" >> Makefile
echo "	rm -f \$(OBJECTS) libgpu-commands.so" >> Makefile

else

echo "all: libgpu-commands.dylib" >> Makefile
echo "libgpu-commands.dylib :  \$(OBJECTS) \$(INCLUDES)" >> Makefile
echo "	$CUDA_PREFIX/bin/nvcc $GENCODE_FLAGS -Xlinker -fexceptions -Xlinker -O3  -shared -o \$@ \$(OBJECTS)" >> Makefile
echo "install:" >> Makefile	
echo "	mkdir -p $INSTALL_PREFIX"  >> Makefile
echo "	cp libgpu-commands.dylib $INSTALL_PREFIX" >> Makefile
echo "clean:" >> Makefile
echo "	rm -f \$(OBJECTS) libgpu-commands.dylib" >> Makefile

fi

echo "libgpu-commands makefile successfully made"
