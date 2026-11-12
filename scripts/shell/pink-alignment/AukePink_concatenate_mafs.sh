#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --job-name=AukePink_concat-mafs
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=patrick.barry@noaa.gov
#SBATCH --output=/home/pbarry/AukePink_lcWGS/job_outfiles/AukePink_concatenate-mafs_%A.out

module unload bio/angsd/0.933
module load bio/angsd/0.933

for i in /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060173.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060174.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060175.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060176.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060177.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060178.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060179.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060180.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060181.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060182.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060183.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060184.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060185.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060186.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060187.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060188.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060189.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060190.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060191.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060192.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060193.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060194.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060195.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060196.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060197.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060198.1_polymorphic.mafs.gz /home/pbarry/AukePink_lcWGS/gls/AukePink_NC_060199.1_polymorphic.mafs.gz
do zcat $i | tail -n +2 -q >> /home/pbarry/AukePink_lcWGS/gls/AukePink_wholegenome_polymorphic.mafs; done
cut -f 1,2,3,4 /home/pbarry/AukePink_lcWGS/gls/AukePink_wholegenome_polymorphic.mafs > /home/pbarry/AukePink_lcWGS/gls/AukePink_wholegenome_polymorphic.sites
gzip /home/pbarry/AukePink_lcWGS/gls/AukePink_wholegenome_polymorphic.mafs

angsd sites index /home/pbarry/AukePink_lcWGS/gls/AukePink_wholegenome_polymorphic.sites
