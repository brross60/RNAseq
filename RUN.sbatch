###Submitting a job in the background
###make a name.sbatch file and add the line below to the very top (don't skip a line) 
#### to run do chmod +x to the file then sbatch name.sbatch to send it to the scheduler.

#!/bin/sh
#SBATCH -A gts-mwhiteley3
#SBATCH -J abscess
#SBATCH --ntasks-per-node=20
#SBATCH -oReport-%j.out 
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Bross60@gatech.edu


###Interactive
salloc -A gts-mwhiteley3 -N1 --ntasks-per-node=20 -t10:00:00


*** change gts-mwhiteley3 to gts-mwhiteley3-biocluster if out of resources for the month


The scheduler will reject a job if the job submission exceeds the 264,960 CPU-hours ÷ Number of Requested Processors.
If your job is listed as complete with no output, please check that the nodes * cores per node * walltime < 
264,960 processor-hours.
1 node * 20cores * 48hrs =960 <<<< 264,960 processor-hours

##https://docs.pace.gatech.edu/phoenix_cluster/slurm_guide_phnx/
##https://docs.pace.gatech.edu/software/anacondaEnv/
