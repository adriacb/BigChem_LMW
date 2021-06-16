#!/bin/tcsh
# @ job_name = analysis
# @ initialdir = .
# @ wall_clock_limit = 24:00:00
# @ output = analys.q.o
# @ error =  analys.q.e
# @ total_tasks = 1
# @ tasks_per_node = 1
# @ cpus_per_task = 1

# marenostrum
module purge
module load intel impi mkl netcdf hdf5 amber python/2.7.13
source /gpfs/projects/ub63/apps/python/bin/activate

set DIR = PYZ_1
#mdmix analyze align byname -s $DIR >& log/log_align_${DIR}.txt
mdmix analyze density byname -s $DIR -C8  >& log/log_density_${DIR}.txt

#set DIR = PYZ_2
#mdmix analyze align byname -s $DIR >& log/log_align_${DIR}.txt
#mdmix analyze density byname -s $DIR -C8  >& log/log_density_${DIR}.txt


#set DIR = PYZ_3
#mdmix analyze align byname -s $DIR >& log/log_align_${DIR}.txt
#mdmix analyze density byname -s $DIR -C8 >& log/log_density_${DIR}.txt


#set DIR = PYZ_1
#mdmix analyze align byname -s $DIR -C4  >& log_align_${DIR}.txt
#mdmix analyze density byname -s $DIR -C8 >& log/log_density_${DIR}.txt
#mdmix analyze energy byname -s $DIR >& log_energy_${DIR}.txt

#set DIR = PYZ_2
#mdmix analyze align byname -s $DIR -C4  >& log_align_${DIR}.txt
#mdmix analyze density byname -s $DIR -C4 >& log_density_${DIR}.txt
#mdmix analyze energy byname -s $DIR >& log_energy_${DIR}.txt

#set DIR = PYZ_3
#mdmix analyze align byname -s $DIR -C4  >& log_align_${DIR}.txt
#mdmix analyze density byname -s $DIR -C4 >& log_density_${DIR}.txt
#mdmix analyze energy byname -s $DIR >& log_energy_${DIR}.txt

#set DIR = WAT_1
#mdmix analyze align byname -s $DIR -C8  >& log/log_align_${DIR}.txt
#mdmix analyze density byname -s $DIR -C8 >& log/log_density_${DIR}.txt
#mdmix analyze energy byname -s $DIR >& log_energy_${DIR}.txt

#set DIR = WAT_2
#mdmix analyze align byname -s $DIR -C4  >& log_align_${DIR}.txt
#mdmix analyze density byname -s $DIR -C8 >& log/log_density_${DIR}.txt
#mdmix analyze energy byname -s $DIR >& log_energy_${DIR}.txt


#set DIR = WAT_3
#mdmix analyze align byname -s $DIR -C4  >& log_align_${DIR}.txt
#mdmix analyze density byname -s $DIR -C8 >& log/log_density_${DIR}.txt
#mdmix analyze energy byname -s $DIR >& log_energy_${DIR}.txt


#mdmix analyze energy bysolvent -s PYZ >& log/log_energy_bysolvent_PYZ.txt
#mdmix analyze energy bysolvent -s MAM >& log/log_energy_bysolvent_MAM.txt
#mdmix analyze energy bysolvent -s WAT >& log/log_energy_bysolvent_WAT.txt


exit 0
