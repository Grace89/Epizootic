#!/bin/bash

NSIM=1000
MAXINTFAM=9

cd ..

# Array of interacting families
array=(2 8)

# Sim ID Counter
sim_id=1

#for i in `seq 1 $MAXINTFAM`;

# Loop across interacting families
for i in "${array[@]}";
        do

      # Loop across super spreader mode
      for super in TRUE FALSE;
        do

        	# for abiotic in TRUE FALSE;
        	# 	do

                for abiotic_trans_rate in 0 0.1 0.01 .001;
                     do

              #  echo $i
              # First argument is number of interacting families
              # Second argument is number of simulations
              # Third argument is whether in super spreader mode or not
              # Fourth argument is sim ID

               nohup Rscript ./Code/sirSim_27_02_2020.R $i $NSIM $super $abiotic_trans_rate $sim_id &

               ## Increment counter for sim ID
               ((sim_id++))
            # done
        done
    done
done
