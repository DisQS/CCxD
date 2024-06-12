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
NOOFSAMPLES=6
NOOFSTEPS=40
OFFSETVAL=0
MULTIPLY_DIVIDE=0
SPIANGLEDTH=0.01
SINGLEANGLEDTH=0.01
SYMMETRISE=0
READIN=0
READINADDRESS=0


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
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3850
#SBATCH --time=08:00:00
#SBATCH --account=su007-rr


module purge
module load GCCcore/12.2.0
module load GCC/12.2.0
module load Eigen/3.4.0
module load CMake/3.24.3

echo "srun ./TRIRG ${NOOFSAMPLES} ${NOOFSTEPS} ${OFFSETVAL} \$1 \$2 ${SYMMETRISE} ${READIN} ${READINADDRESS}"
srun ./TRIRG ${NOOFSAMPLES} ${NOOFSTEPS} ${OFFSETVAL} \$1 \$2 ${SYMMETRISE} ${READIN} ${READINADDRESS};

EOD

chmod 755 ${jobfile}

for i in {0..20}
do
for j in {0..20}
do
sbatch $jobfile $i $j
done;
done;


sleep 1

