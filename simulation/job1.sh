#!/bin/bash
#SBATCH -J sim3job1
#SBATCH  --nodes=1
#SBATCH  --ntasks=30
#SBATCH --cpus-per-task=1
#SBATCH --time=720:00:00
#SBATCH -o simulation3/jobout/job1_out.log
#SBATCH -e simulation3/jobout/job1_err.log

cd /home/project07/Bencong/scratch/project07
source /opt/share/etc/miniconda3-py39.sh
conda activate R4.2


for delta in 0.5 1 1.5 
do
    for prop in 0 0.2 0.4 0.6 0.8  
    do
        srun --ntasks=1 --exclusive R --vanilla --slave --args $delta 500 $prop < simulation3/script/DataGeneratorMOB.R &
        srun --ntasks=1 --exclusive R --vanilla --slave --args $delta 1000 $prop < simulation3/script/DataGeneratorMOB.R &  
    done
done
wait 


for delta in "simulation3/patternMOB1" "simulation3/patternMOB2" "simulation3/patternMOB3"
do
    for prop in 0 0.2 0.4 0.6 0.8  
    do
        srun --ntasks=1 --exclusive R --vanilla --slave --args $delta 500 $prop < simulation3/script/BISON.R &
        srun --ntasks=1 --exclusive R --vanilla --slave --args $delta 500 $prop < simulation3/script/spartaco.R &
    done
done
wait 

for delta in "simulation3/patternMOB1" "simulation3/patternMOB2" "simulation3/patternMOB3"
do
    for prop in 0 0.2 0.4 0.6 0.8  
    do
        srun --ntasks=1 --exclusive R --vanilla --slave --args $delta 500 $prop < simulation3/script/OTHERS.R &
        srun --ntasks=1 --exclusive R --vanilla --slave --args $delta 1000 $prop < simulation3/script/OTHERS.R &  
    done
done
wait 



for delta in "simulation3/patternMOB1" "simulation3/patternMOB2" "simulation3/patternMOB3"
do
    for prop in 0 0.2 0.4 0.6 0.8  
    do
        srun --ntasks=1 --exclusive R --vanilla --slave --args $delta 1000 $prop < simulation3/script/BISON.R &
        srun --ntasks=1 --exclusive R --vanilla --slave --args $delta 1000 $prop < simulation3/script/spartaco.R &
    done
done
wait 


for delta in "simulation3/patternMOB1" "simulation3/patternMOB2" "simulation3/patternMOB3"
do
    for prop in 0 0.2 0.4 0.6 0.8  
    do
        srun --ntasks=1 --exclusive R --vanilla --slave --args $delta 500 $prop < simulation3/script/summary_spartaco.R &
        srun --ntasks=1 --exclusive R --vanilla --slave --args $delta 1000 $prop < simulation3/script/summary_spartaco.R &  
    done
done
wait 

for delta in "simulation3/patternMOB1" "simulation3/patternMOB2" "simulation3/patternMOB3"
do
    for prop in 0.2 0.4 0.6 0.8  
    do
        srun --ntasks=1 --exclusive R --vanilla --slave --args $delta 500 $prop < simulation3/script/accuracy_bison.R &
        srun --ntasks=1 --exclusive R --vanilla --slave --args $delta 1000 $prop < simulation3/script/accuracy_bison.R &  
    done
done
wait 




