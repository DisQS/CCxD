#!/bin/bash
#
# run as: for size in 10 12; do echo $size; ../../PBS/L31-conE-rVal.sh 1.0 $size 1000 0; done

# settings from input


configs=${1}

echo "C-Dists: making distributions for= " $configs " steps"

# settings for directories

currdir=`pwd`'/../data'
cd $currdir
jobdir="C-RG-$configs"
mkdir -p $jobdir

jobname=$jobdir-$maxRGSteps
echo $jobname

jobfile=`printf "$jobname.sh"`
wlsfile=`printf "$jobname.wls"`
logfile=`printf "$jobname.log"`

# settings for parallel submission

cd $jobdir

cat > ${wlsfile} << EOD
#!/usr/bin/env wolframscript 
Print["(*Preliminaries*)"];
maindir="$currdir";


Print[maindir];
Print["--- FINISHED!"];
EOD

cat > ${jobfile} << EOD
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2012
#SBATCH --time=48:00:00
module purge
module load Mathematica
module list
pwd
ls -al $jobdir/$wlsfile
$jobdir/$wlsfile
EOD

cd ..
chmod 755 ${jobdir}/${jobfile}
chmod 755 ${jobdir}/${wlsfile}
##(sbatch -q devel $jobdir/${jobfile}) # for queueing system
sbatch ${jobdir}/${jobfile} # for queueing system
#(source ${jobdir}/${jobfile} ) >& ${jobdir}/${logfile} & # for parallel shell execution

#echo "<return>"
sleep 1

cd ..
