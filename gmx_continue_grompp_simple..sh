#!/bin/bash -e 

# This script continues a GROMACS MD simulation in segments, starting from an initial segment (0 to dt ns) 
# and then proceeding in increments of dt until tend ns is reached. IO files from the previous segment 
# are used as input for the next segment. 
# The intial segment is run only if tstart is 0, otherwise the script assumes that the initial segment 
# has already been run and starts from the continuation segments.
# The script also checks for a STOP file after each segment, and if found, it halts further processing.
# The script assumes that the initial segment's mdp file is named ${base_name}_0_0.mdp, and that the
# output files from each segment are named in the format ${base_name}_${start_time}_${end_time}.* 
# (e.g. .gro, .cpt, .edr, .mdp, .tpr)
# The initial segment is run only if the users sets the start tine to 0, 
# otherwise the script assumes that the initial segment has already been run and starts from the 
# continuation segments.
# This scritp does not use the "mdrun -cpi" option, but instead uses the .cpt file from the previous 
# segment as input for grompp, which is a common approach for continuing simulations in segments.
# (The Gromacs manual recommends using the -cpi option for continuing simulations, but it is not strictly 
# necessary as long as the .cpt file is used correctly in grompp.)
# The script also corrects the in each segment's mdp file to set the correct tinit value, which is important
# for proper continuation of the simulation.
# This script does not handle interrupted simulations or restarts, but it can be easily modified to 
# include such functionality if needed. 
if ! command -v gmx &> /dev/null; then
    . /programs/gromacs-2025.2/bin/GMXRC.bash
fi

# User-defined variables
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
    echo -e "\n****** Running initial segment: 0 ns to ${dt} ns ******\n"
    # Create mdp file for initial segment from template
    sed  -e "s/\(steps\s\+= \)[0-9]\+/\1${steps}/" ${base_name}_0_0.mdp > ${base_name}_0_${dt}.mdp
    gmx grompp -f ${base_name}_0_${dt}.mdp -c ${base_name}_0_0.gro -t ${base_name}_0_0.cpt -e ${base_name}_0_0.edr -p topol.top -o ${base_name}_0_${dt}.tpr
    gmx mdrun -deffnm ${base_name}_0_${dt}
    # Increment start time for following segments 
    tstart=$((tstart+dt))
fi

# Continuation segments
for ((cur_time=$tstart ; cur_time<$tend ; cur_time+=$dt)) ; do

    echo -e "\n****** Running segment: ${cur_time} ns to $((cur_time+dt)) ns ******\n"

    prev=${base_name}_$((cur_time-dt))_${cur_time}

    if [[ ! -f ${prev}.mdp ]] ; then
        echo "Some job files are missing. Aborting."
        break
    fi
    cur=${base_name}_${cur_time}_$((cur_time+dt))

    timeps=$((1000*cur_time)) # time in mdp is in picoseconds ! 

    # Fix tinit for mdp file of current segment
    sed  -e "s/\(tinit\s\+= \)[0-9]\+/\1${timeps}/" ${prev}.mdp > ${cur}.mdp

    gmx grompp -f ${cur}.mdp -c ${prev}.gro -t ${prev}.cpt -e ${prev}.edr -p topol.top -o ${cur}.tpr
    gmx mdrun -deffnm $cur

    # Check for STOP file to halt further processing
    if [[ -f ./STOP ]] ; then
        echo "Found STOP file. Stopping."
        break
    fi    
done
