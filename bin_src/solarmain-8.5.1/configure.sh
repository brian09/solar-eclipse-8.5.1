#!/bin/bash
SCRIPT=$0
SCRIPT_PATH=$(dirname $SCRIPT)
cd $SCRIPT_PATH
SCRIPT_PATH=$(pwd)
args=("$@")
BUILD_DYNAMIC=1
USE_GPU_TOOLS=0
CUDA_LIBGPUCOMMANDS=0
CUDA_LIBCUBLAS=0
CUDA_LIBCUDART=0
INSTALL_PATH=$SCRIPT_PATH
LIBRARY_PATH=$SCRIPT_PATH/lib
INCLUDE_PATH=$SCRIPT_PATH/include
USE_MKL=0
MKL_PATH=/opt/intel/mkl
CXX_COMPILER=$(which g++)
C_COMPILER=$(which gcc)
FC_COMPILER=$(which gfortran)
shared_ext="so"
OS="Linux"
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
 	echo "Usage: configure.sh [Options]"
 	echo "Script Options:"
 	echo ""
 	echo "--install_path <Path where solarmain is installed> Default: $INSTALL_PATH"
 	echo " " 
 	echo "--cxx_compiler <C++ compiler> Default: $(which g++)"
 	echo " " 
  	echo "--c_compiler <C compiler> Default: $(which gcc)"
 	echo " "
  	echo "--fc_compiler <Fortran compiler> Default: $(which gfortran)"
 	echo " "
 	echo "--lib_path <Library search path> Default: $LIBRARY_PATH"
 	echo " "
 	echo "--include_path <Path to directory containing header files> Default: $INCLUDE_PATH" 	
 	echo ""
 	echo "--static Compiles solarmain static Default: Disabled"
 	echo "Note: If compiled static GPU tools are disabled or if not built on Linux"
 	echo ""
 	echo "--enable_gpu_tools <Path to libgpu-commands.$shared_ext> <Path to libcudart.$shared_ext <Path to libcublas.$shared_ext>  Default: Disabled "
 	echo "Enables CUDA GPU commands gpu_pedifromsnps,gpu_fphi, and gpu_gwas"
 	echo "Require installation of CUDA Toolkit and NVIDIA CUDA capable GPU with architecture version greater than or equal to 3.5"
 	echo ""
 	echo "--use_mkl <Path to MKL include and lib folders> Default: Disabled" 
 	echo "By default Intel Math Kernel Library is installed in /opt/intel which means"
 	echo "that the path would be /opt/intel/mkl.  The path must include an include folder"
 	echo "and a lib folder containing: "
 	echo "  intel64/libmkl_core.a"
 	echo "  intel64/libmkl_gnu_thread.a," 
 	echo "  intel64/libmkl_intel_lp64.a" 
 	echo ""

}

OS=$(uname)
for ((i = 0; i < $#; i++ )); do
   
     current_arg=${args[$i]}
     
 if [ "$current_arg" = "--help" ] || [ "$current_arg" = "--help" ]
 then
	help_display
	exit 0
 elif [ "$current_arg" = "--static" ]  
 then
 	if [ "$OS" != "Linux" ]
 	then	
 		echo "Operating System: $OSTYPE must be linux-gnu to use the -static option"
 		exit 1
 	fi
 	BUILD_DYNAMIC=0	
 elif [ "$current_arg" = "--enable_gpu_tools" ] && [ $(($i + 3)) -lt $# ] 
 then
 	USE_GPU_TOOLS=1
 	i=$(($i + 1))
 	CUDA_LIBGPUCOMMANDS=${args[$i]}
 	i=$(($i + 1))
 	CUDA_LIBCUDART=${args[$i]}
  	i=$(($i + 1))
 	CUDA_LIBCUBLAS=${args[$i]}	 		
 elif [ "$current_arg" = "--install_path" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	INSTALL_PATH=${args[$i]} 
  elif [ "$current_arg" = "--lib_path" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	LIBRARY_PATH=${args[$i]}	
 elif [ "$current_arg" = "--include_path" ] && [ $(($i + 1)) -lt $# ] 
 then
 	i=$(($i + 1))
 	INCLUDE_PATH=${args[$i]} 
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
if [ $BUILD_DYNAMIC -eq 0 ] 
then
    echo "GPU tools have been disabled because --static option was used"
	USE_GPU_TOOLS=0
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

cxx_version=$($CXX_COMPILER -dumpversion)	
if [ $cxx_version -lt 5 ]
then
 	echo "Warning: C++ compiler version less than 5.0"
 	echo "Compilation errors may occur"
elif [ $cxx_version -gt 7 ]
then
 	echo "Warning C++ compiler version is greater 7"
 	echo "Runtime errors may occur"
fi

c_version=$($C_COMPILER -dumpversion)
if [ $c_version -lt 5 ] 
then
 	echo "Warning: C compiler version is less than 5.0"
 	echo "Compilation errors may occur"
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

if [ $BUILD_DYNAMIC -eq 1 ] 
then
	if [ ! -f $LIBRARY_PATH/libplinkio.a ] && [ ! -f $LIBRARY_PATH/libplinkio.$shared_ext ]
	then 
    		echo "Library libplinkio.a or libplinkio.$shared_ext not found in $LIBRARY_PATH"
    		exit 1
	fi
else
	if [ ! -f $LIBRARY_PATH/libplinkio.a ]
	then 
    		echo "Library libplinkio.a not found in $LIBRARY_PATH"
    		exit 1
	fi
fi
	
if [ ! -f $LIBRARY_PATH/libtcl8.4.a ] 
then 
    echo "Library libtcl8.4.a not found in $LIBRARY_PATH"
    exit 1
fi
if [ ! -f $LIBRARY_PATH/libRicVolume.a ] 
then 
    echo "Library libRicVolume.a not found in $LIBRARY_PATH"
    exit 1
fi
if [ ! -f $LIBRARY_PATH/libsafe.a ] 
then 
    echo "Library libsafe.a not found in $LIBRARY_PATH"
    echo "If you're having trouble finding this library compile safelib.c"
    echo "into a static library.  By default it can be found in the src directory"
    exit 1
fi
if [ $USE_MKL -eq 1 ] 
then
    if  [ ! -d $MKL_PATH/include ]
    then
        echo "Folder $MKL_PATH/include not found"
        exit 1
    fi
    if [ ! -f $MKL_PATH/lib/intel64/libmkl_gnu_thread.a ] 
    then 
        echo "Library libmkl_gnu_thread.a not found in $MKL_PATH/lib/intel64"
        exit 1
    fi
    if [ ! -f $MKL_PATH/lib/intel64/libmkl_intel_lp64.a ] 
    then 
        echo "Library libmkl_intel_lp64.a not found in $MKL_PATH/lib/intel64"
        exit 1
    fi
    if [ ! -f $MKL_PATH/lib/intel64/libmkl_core.a ] 
    then 
        echo "Library libmkl_core.a not found in $MKL_PATH/lib/intel64"
        exit 1
    fi
fi
if [ ! -f $CUDA_LIBGPUCOMMANDS ] && [ $USE_GPU_TOOLS -eq 1 ]
then 
    echo "Library libgpu-commands.$shared_ext not found at $CUDA_LIBCOMMANDS"
    exit 1
fi
if [ ! -f $CUDA_LIBCUBLAS ] && [ $USE_GPU_TOOLS -eq 1 ]
then 
    echo "Library libcublas.$shared_ext not found at $CUDA_LIBCUBLAS"
    exit 1
fi
if [ ! -f $CUDA_LIBCUDART ] && [ $USE_GPU_TOOLS -eq 1 ]
then 
    echo "Library libcudart.$shared_ext not found at $CUDA_LIBCUDART"
    exit
fi
rm -f Makefile
rm -f sources.mk
echo "SOURCE_PATH=$SCRIPT_PATH/src" >> sources.mk
echo "" >> sources.mk
echo "CXX_OBJECTS= \\" >> sources.mk
echo "\$(SOURCE_PATH)/alnormc.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/annotate_gwas.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/ccsearch.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/calculate_snp_frequencies.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/chi.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/constraint.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/covariate.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/create_evd.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/create_fake_pedigree.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/define.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/drand.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/eqvar.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/evd.o \\" >> sources.mk
#echo "\$(SOURCE_PATH)/evdlikc.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/expression.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/field.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/fisherpedigree.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/fphi.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/freq.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/function.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/gwas.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/help.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/howclose.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/ibd.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/ibdoption.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/ibs.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/inorm_nifti.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/key.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/listmaker.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/loadsave.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/loglike.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/map.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/marker.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/mathmatrix.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/matrix.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/maximize.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/mibd.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/model.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/mu.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/nifti_assemble.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/nifti_to_csv.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/normal.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/omega.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/option.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/parameter.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/pedifromsnps.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/pedigree.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/phenotypes.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/plink_converter.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/power.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/reorder_phenotype.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/scale.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/scan.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/scratchfile.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/simqtl.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/snp.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/solar.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/solarfile.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/solarfilecmd.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/solar_mle_setup.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/solar-trait-reader.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/split_pheno_file.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/sporadic.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/tablefile.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/tablefilecmd.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/token.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/trait.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/verbosity.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/voxel.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/zscore.o \\" >> sources.mk
echo "" >> sources.mk
echo "C_OBJECTS= \\" >> sources.mk
echo "\$(SOURCE_PATH)/ccuserid.o \\" >> sources.mk
#echo "\$(SOURCE_PATH)/cpedifromsnps.o \\" >> sources.mk
#echo "\$(SOURCE_PATH)/cpedifromsnps_old.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/cplink_converter.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/eigstruc.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/fdate.o \\" >> sources.mk
#echo "\$(SOURCE_PATH)/getusername.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/mehd.o \\" >> sources.mk
#echo "\$(SOURCE_PATH)/mt.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/ncdf.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/npdf.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/nrutil.o \\" >> sources.mk
#echo "\$(SOURCE_PATH)/pipeback.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/plinkio.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/plotpipe.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/safelib.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/solarmain.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/tclgr.o \\" >> sources.mk
#echo "\$(SOURCE_PATH)/testuidt.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/w_mehd.o \\" >> sources.mk
echo "" >> sources.mk
echo "FC_OBJECTS= \\" >> sources.mk
echo "\$(SOURCE_PATH)/alnorm.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/astop.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/asycov.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/beta.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/calc.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/cdfchi.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/cdfnor.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/chisq.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/clear.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/copyped.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/covar.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/dasycov.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/dcopyped.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/dcovar.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/ddfun.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/delta7.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/dfun.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/digam.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/direct.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/dmean.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/dppfa.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/edftst.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/evdlik.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/evdout.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/fdist.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/fit2dp.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/fortfiles.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/fpedifromsnps.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/fun.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/fun_mehd.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/gamlog.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/gaucdf.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/hessin.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/hpsort.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/inital.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/ipmpar.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/kin.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/linpack.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/logo.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/loop.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/mblank.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/mvncdf.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/newlik.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/optima.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/output.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/pedtst.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/phidens.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/pinput.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/ppnd.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/preopt.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/preped.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/qdprog.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/random.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/resid.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/scor.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/search.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/settab.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/smpoutput.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/spmpar.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/sweep.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/symeig.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/tdist.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/trigam.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/upcase.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/writestr.o \\" >> sources.mk
echo "\$(SOURCE_PATH)/yesno.o \\" >> sources.mk
echo "" >> sources.mk
echo "HEADERS= \\" >> sources.mk
echo "\$(SOURCE_PATH)/config.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/expression.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/mvn.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/nrutil.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/plotpipe.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/solar.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/solar_mle_setup.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/solar-trait-reader.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/tablefile.h \\" >> sources.mk
echo "\$(SOURCE_PATH)/token.h \\" >> sources.mk
echo "" >> sources.mk
echo "LIBRARIES= \\" >> sources.mk
echo "$LIBRARY_PATH/libsafe.a \\" >> sources.mk
echo "$LIBRARY_PATH/libRicVolume.a \\" >> sources.mk
if [ $BUILD_DYNAMIC -eq 1 ] && [ ! -f $LIBRARY_PATH/libplinkio.a ]
then
	echo "$LIBRARY_PATH/libplinkio.$shared_ext \\" >> sources.mk
else
	echo "$LIBRARY_PATH/libplinkio.a \\" >> sources.mk
fi
echo "$LIBRARY_PATH/libtcl8.4.a \\" >> sources.mk
if [ $USE_GPU_TOOLS -eq 1 ]
then
	echo "$CUDA_LIBGPUCOMMANDS \\" >> sources.mk
	echo "$CUDA_LIBCUBLAS \\" >> sources.mk
	echo "$CUDA_LIBCUDART \\" >> sources.mk
fi
echo "-ldl \\" >> sources.mk
echo "-lstdc++ \\" >> sources.mk
if [ $USE_MKL -eq 1 ]
then
    echo "-lgomp \\" >> sources.mk
    echo "" >> sources.mk
    echo "MKL_LIBS= \\" >> sources.mk
    echo "$MKL_PATH/lib/intel64/libmkl_core.a \\" >> sources.mk
    echo "$MKL_PATH/lib/intel64/libmkl_gnu_thread.a \\" >> sources.mk
    echo "$MKL_PATH/lib/intel64/libmkl_intel_lp64.a \\" >> sources.mk
else
    echo "" >> sources.mk
    echo "MKL_LIBS= \\" >> sources.mk
fi

echo "include sources.mk" >> Makefile
echo ".SUFFIXES: .f .cc .o .h" >> Makefile
echo "FFLAGS=-m64  -O2 -fno-second-underscore -fexceptions" >> Makefile
echo "CFLAGS=-std=c99 -m64 -O2 -I$INCLUDE_PATH -DUSE_SAFELIB -fexceptions" >> Makefile
MKL_INCLUDE=
if [ $USE_MKL -eq 1 ]
then
    MKL_INCLUDE="-I$MKL_PATH/include -DEIGEN_USE_MKL_ALL"
fi
if [ $USE_GPU_TOOLS -eq 1 ] 
then
	echo "CXXFLAGS= -O2  -std=c++11 -fopenmp  -DGPU_TOOLS  -m64 -I../lib -I$INCLUDE_PATH -DUSE_SAFELIB -fexceptions -DTR1 -DLINUX64 -DUSE_SAFELIB  $MKL_INCLUDE" >> Makefile
else
	echo "CXXFLAGS=-O2  -std=c++11 -fopenmp  -m64 -I../lib -I$INCLUDE_PATH -DUSE_SAFELIB -fexceptions -DTR1 -DLINUX64 -DUSE_SAFELIB  $MKL_INCLUDE" >> Makefile
fi

if [ $BUILD_DYNAMIC -eq 1 ]
then
	echo "LDFLAGS=-std=c++11 -m64 -O2 -L$LIBRARY_PATH -pthread -fopenmp -fexceptions" >> Makefile
else
	echo "LDFLAGS=-static -std=c++11 -m64 -O2 -L$LIBRARY_PATH -pthread -fopenmp -fexceptions" >> Makefile
fi

echo ".cc.o:" >> Makefile
echo "	$CXX_COMPILER -c \$(CXXFLAGS) -o \$@ \$<" >> Makefile
echo ".f.o:" >> Makefile
echo "	$FC_COMPILER -c \$(FCFLAGS) -o \$@ \$<" >> Makefile
echo ".c.o:" >> Makefile
echo "	$C_COMPILER -c \$(CFLAGS) -o \$@ \$<" >> Makefile
echo ".PHONY: all" >> Makefile
echo "all: solarmain" >> Makefile
echo "solarmain : \$(FC_OBJECTS) \$(CXX_OBJECTS) \$(C_OBJECTS)  \$(HEADERS)" >> Makefile
echo "	$FC_COMPILER -o \$@  \$(FC_OBJECTS) \$(CXX_OBJECTS) \$(C_OBJECTS) \$(LDFLAGS) \$(LIBRARIES) -Wl,--start-group \$(MKL_LIBS) -Wl,--end-group" >> Makefile
echo "install:" >> Makefile
echo "	cp solarmain $INSTALL_PATH" >> Makefile
echo "clean:" >> Makefile
echo "	rm -f \$(FC_OBJECTS) \$(CXX_OBJECTS) \$(C_OBJECTS)  solarmain" >> Makefile
echo "solarmain makefile successfully made"
exit 0
