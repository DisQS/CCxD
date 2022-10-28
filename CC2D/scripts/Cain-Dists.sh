#!/bin/bash
#
# run as: for size in 10 12; do echo $size; ../../PBS/L31-conE-rVal.sh 1.0 $size 1000 0; done

# settings from input



size=$1
echo "C-Dists: making distributions for= " $size " configs" 

# settings for directories

currdir=`pwd`'/../data'
cd $currdir
jobdir="C-RG-$size"


jobname="C-Dists-$size"
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
inputFiles = FileNames[RegularExpression[".*nc"]];
Print["Found ", Length[inputFiles], " nc files to import data from"]
inputFileNumbers = {}
Do[If[StringPart[inputFiles[[i]], -8] == "_", 
   AppendTo[inputFileNumbers, StringPart[inputFiles[[i]], -7]], 
   AppendTo[inputFileNumbers, 
    StringTake[inputFiles[[i]], {-8, -7}]]], {i, 1, 
   Length[inputFiles]}];
Print[inputFileNumbers]
Do[Print["Generating distributions for step number " <> 
   inputFileNumbers[[i]]], {i, 1, 1}]
SetDirectory["$jobdir"];
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
#sbatch ${jobdir}/${jobfile} # for queueing system
(source ${jobdir}/${jobfile} ) >& ${jobdir}/${logfile} & # for parallel shell execution

#echo "<return>"
sleep 1

cd ..
