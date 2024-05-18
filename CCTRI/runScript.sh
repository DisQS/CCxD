#!/bin/bash
#argument numbers
#    noOfSamples -> argv[1]
#    noOfSteps -> argv[2]
#    offsetVal -> argv[3]
#    multiply/divide -> argv[4]
#    spinangle -> argv[4]
#    symmetrise -> argv[5]
#    readIn -> argv[6]
#    readInAddress -> argv[7]
NO_OF_SAMPLES=7
NOOFSTEPS=30
OFFSETVAL=0
MULTIPLY_DIVIDE=0
SPINANGLESTART=0
SPINANGLEEND=0.49
SPIANGLEDTH=0.05
SYMMETRISE=0
READIN=0
READINADDRESS=0

currdir=`pwd`'/../data'
cd $currdir
jobdir="TRIRG-$NOOFSAMPLES-$NOOFSTEPS"
mkdir -p $jobdir

jobfile=`printf "$jobdir.slurm"`
logfile=`printf "$jobdir.log"`


cd $jobdir

cat > ${jobfile} << EOD
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3850
#SBATCH --time=08:00:00
#SBATCH --account=su007-rr
#SBATCH --array=0-9


module purge
module load GCC/12.2.0
module load CMake/3.22.1
module load Eigen/3.4.0

cmake ../CCTRI/
cmake --build .

srun ${jobdir}/TRIRG ${NOOFSAMPLES} ${NOOFSTEPS} ${OFFSETVAL} ${MULTIPLY_DIVIDE} $(($SPINANGLESTART+($SPINANGLEDTH*$SLURM_ARRAY_TASK_ID))) ${SYMMETRISE} ${READIN} ${READINADDRESS}

EOD

sbatch $jobfile
