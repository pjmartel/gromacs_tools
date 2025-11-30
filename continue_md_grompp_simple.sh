#!/bin/bash -e 
if ! command -v gmx &> /dev/null; then
    . /programs/gromacs-2025.2/bin/GMXRC.bash
fi

prefix=md
replica=0
base_name=${prefix}_${replica}

# Times are in ns
tstart=0
tend=500
dt=100  
# Steps for given time interval
steps=$((dt*500000))  # assumes a 2fs timestep

# Initial segment.  Run once, for tstart=0
if [ $tstart -eq 0 ]; then
    sed  -e "s/\(steps\s\+= \)[0-9]\+/\1${steps}/" ${base_name}_0_0.mdp > ${base_name}_0_${dt}.mdp
    gmx grompp -f ${base_name}_0_${dt}.mdp -c ${base_name}_0_0.gro -t ${base_name}_0_0.cpt -e ${base_name}_0_0.edr -p topol.top -o ${base_name}_0_${dt}.tpr
    gmx mdrun -deffnm ${base_name}_0_${dt} 
    tstart=$((tstart+dt))
fi

# following segments run on loop
for ((cur_time=$tstart ; cur_time<$tend ; cur_time+=$dt)) ; do

    echo -e "\n****** Running segment: ${cur_time} ns to $((cur_time+dt)) ns ******\n"

    prev=${base_name}_$((cur_time-dt))_${cur_time}

    if [[ ! -f ${prev}.mdp ]] ; then
        echo "Some job files are missing. Aborting."
        break
    fi
    cur=${base_name}_${cur_time}_$((cur_time+dt))

    timeps=$((1000*cur_time)) # time in mdp is in picoseconds ! 

    sed  -e "s/\(tinit\s\+= \)[0-9]\+/\1${timeps}/" ${prev}.mdp > ${cur}.mdp

    gmx grompp -f ${cur}.mdp -c ${prev}.gro -t ${prev}.cpt -e ${prev}.edr -p topol.top -o ${cur}.tpr
    gmx mdrun -deffnm $cur
    if [[ -f ./STOP ]] ; then
        echo "Found STOP file. Stopping."
        break
    fi    
done
