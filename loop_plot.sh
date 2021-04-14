#!/bin/bash
# Use postprocess_v01_kinetic_energy_spectra.py to process multiple data files in a loop
# Wave types: ROS, IG, TOT, Era5
# Usage: 1) Change the input file path: inp_fn, according to the output file path set in loop_run_sph.sh or run_sph.py
#        2) Run the script in terminal: /bin/bash loop_plot.sh



### ROS, IG, TOT, Era5
date -d "2020-02-15 00  + 6 hour" +"%Y%m%d%H"
start="2020-02-15 00"
for wave in "ROS" "IG" "TOT" "era5"
do
    for i in {0..27}    
    do
        h=$[6*$i]
        echo $wave" loop "$i
        time_s=$(date -d "2020-02-15 00  + $h hour" +"%Y%m%d%H") 
        inp_fn="../../../../../../scratch/cen/mi-theo/u300924/sph/sph_pyspharm_oper/data/202010/sph_"$wave"_"$time_s".nc"    
        echo "Calling postprocess_v01_kinetic_energy_spectra.py for: "$inp_fn
        python postprocess_v01_kinetic_energy_spectra.py $inp_fn
        
    done
done
