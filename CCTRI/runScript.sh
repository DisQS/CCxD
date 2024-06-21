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
NOOFSTEPS=40
OFFSETVAL=0
MULTIPLY_DIVIDE=0
SPIANGLEDTH=0.01
SINGLEANGLEDTH=0.01
SYMMETRISE=1
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
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3850
#SBATCH --time=08:00:00
#SBATCH --account=su007-rr


module purge
module load GCCcore/11.3.0
module load GCC/11.3.0
module load Eigen/3.4.0
module load CMake/3.24.3
module load parallel/20220722


touch inputs.txt
echo " ${NOOFSAMPLES} ${NOOFSTEPS} ${OFFSETVAL} 0 0 ${SYMMETRISE} ${READIN} ${READINADDRESS} " > "inputs.txt"
for i in {1..20}
do
for j in {1..20}
do
echo " ${NOOFSAMPLES} ${NOOFSTEPS} ${OFFSETVAL} \$i \$j ${SYMMETRISE} ${READIN} ${READINADDRESS} " >> "inputs.txt"
done;
done;


MY_PARALLEL_OPTS="-N 1 --delay .2 -j \$SLURM_NTASKS --joblog parallel-\${SLURM_JOBID}.log -a inputs.txt"
MY_SRUN_OPTS="-N 1 -n 1 --exclusive"
echo "srun ./TRIRG ${NOOFSAMPLES} ${NOOFSTEPS} ${OFFSETVAL} 0 0 ${SYMMETRISE} ${READIN} ${READINADDRESS}"
parallel --dryrun \$MY_PARALLEL_OPTS srun \$MY_SRUN_OPTS ./TRIRG 

EOD

chmod 755 ${jobfile}


sbatch $jobfile


sleep 1

