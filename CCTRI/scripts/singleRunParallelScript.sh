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
#module load CMake/3.24.3

#Initialise variables
NOOFSAMPLES=7
NOOFSTEPS=33
OFFSETVAL=0
MULTIPLY_DIVIDE=0
SPIANGLEDTH=0
SINGLEANGLEDTH=0
SYMMETRISE=1
READIN=0
READINADDRESS=0
IDNO=0


currdir=`pwd`
cd $currdir
jobdir="TRIRG-$NOOFSAMPLES-$NOOFSTEPS-$SYMMETRISE"
mkdir -p $jobdir

jobfile=`printf "$jobdir.sh"`
logfile=`printf "$jobdir.log"`



cd $jobdir
cp ../TRIRG ./TRIRG

cat > ${jobfile} << EOD
#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3850
#SBATCH --time=08:00:00
#SBATCH --account=su007-rr


module purge
module load GCCcore/11.3.0
module load GCC/11.3.0
module load GCC/11.3.0 OpenMPI/4.1.4
module load Eigen/3.4.0
module load CMake/3.24.3
module load parallel/20220722


srun ./TRIRG "$NOOFSAMPLES $NOOFSTEPS $OFFSETVAL $SINGLEANGLEDTH $SPIANGLEDTH $SYMMETRISE $READIN $READINADDRESS" 

EOD

chmod 755 ${jobfile}


sbatch $jobfile


sleep 1

