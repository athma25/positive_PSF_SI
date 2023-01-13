#!/bin/bash

# Input batch ID
if [ $# -ne 1 ]
then
	echo Give all inputs
	exit 1
fi

dir="./output/$1"

if [ ! -d ${dir} ]
then
	mkdir ${dir}
else
	echo Directory exists. Backup and clear
	exit 2
fi

# Constant (across the sweep) parameters
dx=5e-2
l=2
dt=5e-4
T=10

# Create full factorial list of parameter files based on the following ranges
eta_0=0
grd1=1
D1=0.1

n0=${#eta_0[@]}
n1=${#grd1[@]}
n2=${#D1[@]}

for i0 in $(seq 0 $((n0-1)))
do
	for i1 in $(seq 0 $((n1-1)))
	do
		for i2 in $(seq 0 $((n2-1)))
		do
			fname=$(printf "run%03d.txt" $((i0+n0*i1+n0*n1*i2)))
			echo -e "\"Model parameters\" \t \"\" \t \"\"" > $fname
			echo -e "n \t 2 \t \"Number of species\"" >> $fname
			echo -e "eta_0 \t ${eta_0[i0]} \t \"Strength of external driver\"" >> $fname
			echo -e "\n \"Initial conditions\" \t \"\" \t \"\"" >> $fname
			echo -e "grd \t 0 \t \"Type of gradient - 0: linear, 1: custom\"" >> $fname
			echo -e "grd1 \t ${grd1[i1]} \t \"Gradient parameter 1\"" >> $fname
			echo -e "\n \"Numerics - parameters\" \t \"\" \t \"\"" >> $fname
			echo -e "dx \t ${dx} \t \"Mesh size\"" >> $fname
			echo -e "l \t $l \t \"Domain\"" >> $fname
			echo -e "dt \t ${dt} \t \"Time step\"" >> $fname
			echo -e "T \t $T \t \"Total time\"" >> $fname
			
			fname=$(printf "run%03dm.txt" $((i0+n0*i1+n0*n1*i2)))
			echo -e "\"Model parameters\" \t \"\" \t \"\"" > $fname
			echo -e "ep \t 0 \t 2.5 \t \"Optimal soil conditions\"" >> $fname
			echo -e "sig \t 1 \t 1 \t \"Niche width\"" >> $fname
			echo -e "eta \t 0 \t 0 \t \"Strength of conditioning\"" >> $fname
			echo -e "D \t ${D1[i2]} \t 0.1 \t \"Disperal constants\"" >> $fname
			echo -e "sig_d \t 0.5 \t 0.5 \t \"Dispersal range\"" >> $fname
		done
	done
done

mv run[0-9]*.txt ${dir}

ls ${dir} > sweep.txt
mv sweep.txt ${dir}

echo Generated parameter files at ${dir}
echo Last file is $fname
exit 0
