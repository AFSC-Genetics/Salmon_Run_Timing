#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --job-name=euclide_concat-mafs
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/runtiming/sock-sock/job_outfiles/euclide_concatenate-mafs_%A.out

module unload bio/angsd/0.933
module load bio/angsd/0.933

for i in /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042535.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042536.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042537.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042538.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042539.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042540.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042541.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042542.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042543.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042544.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042545.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042546.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042547.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042548.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042549.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042550.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042551.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042552.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042553.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042554.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042555.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042556.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042557.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042558.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042559.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042560.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042561.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042562.1_minInd0.3.mafs.gz /scratch2/nhowe/runtiming/sock-sock/gls/euclide_NC_042563.1_minInd0.3.mafs.gz
do zcat $i | tail -n +2 -q >> /scratch2/nhowe/runtiming/sock-sock/gls/euclide_wholegenome_minInd0.3.mafs; done
cut -f 1,2,3,4 /scratch2/nhowe/runtiming/sock-sock/gls/euclide_wholegenome_minInd0.3.mafs > /scratch2/nhowe/runtiming/sock-sock/gls/euclide_wholegenome_minInd0.3.sites
gzip /scratch2/nhowe/runtiming/sock-sock/gls/euclide_wholegenome_minInd0.3.mafs

angsd sites index /scratch2/nhowe/runtiming/sock-sock/gls/euclide_wholegenome_minInd0.3.sites
