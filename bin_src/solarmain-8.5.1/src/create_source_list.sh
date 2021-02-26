#!/bin/bash
rm -f files_script.sh
array=($(ls *.cc))
echo "echo \"SOURCE_PATH=$SCRIPT_PATH/src\" >> sources.mk" >> files_script.sh
echo "echo \"CXX_OBJECTS= \\ \" >> sources.mk" >> files_script.sh
for i in "${array[@]}"
do
	echo "echo \"\$(SOURCE_PATH)/$i \\ \" >> sources.mk" >> files_script.sh
	
done
echo "echo \"C_OBJECTS= \\ \" >> sources.mk" >> files_script.sh
array=($(ls *.c))
for i in "${array[@]}"
do
	echo "echo \"\$(SOURCE_PATH)/$i \\ \" >> sources.mk" >> files_script.sh
done
echo "echo \"FC_OBJECTS= \\ \" >> sources.mk" >> files_script.sh
array=($(ls *.f))
for i in "${array[@]}"
do
	echo "echo \"\$(SOURCE_PATH)/$i \\ \" >> sources.mk" >> files_script.sh
done
echo "echo \"HEADERS= \\ \" >> sources.mk" >> files_script.sh
array=($(ls *.h))
for i in "${array[@]}"
do
	echo "echo \"\$(SOURCE_PATH)/$i \\ \" >> sources.mk" >> files_script.sh
done
