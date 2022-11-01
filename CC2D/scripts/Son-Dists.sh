#!/bin/bash
#
# run as: for size in 10 12; do echo $size; ../../PBS/L31-conE-rVal.sh 1.0 $size 1000 0; done

# settings from input



size=$1
echo "S-Dists: making distributions for= " $size " configs" 

# settings for directories

currdir=`pwd`'/../data'
cd $currdir
jobdir="S-RG-$size"


jobname="S-Dists-$size"
echo $jobname

jobfile=`printf "$jobname.sh"`
wlsfile=`printf "$jobname.wls"`
logfile=`printf "$jobname.log"`

# settings for parallel submission

cd $jobdir

cat > ${wlsfile} << EOD
#!/usr/bin/env wolframscript 
Print["(*Preliminaries*)"];
maindir="$currdir"
SetDirectory["$currdir/$jobdir"]
Print[FileNames[]]
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
   inputFileNumbers[[i]]];
	thraw = Import["C-"<>"$size"<>"_"<>ToString[inputFileNumbers[[i]]]<>"raw.nc","Data"];
	thdist = BinCounts[thraw[[i]],{0,1,0.001}];
	traw = Cos[thraw];
	graw = traw^2;
	gdist = BinCounts[graw,{0,1,0.001}];
	zraw = Log[(1/traw)-1];
	qdist = BinCounts[zraw,{-8,8,0.004}];

	Export["S-THdist-" <> "$size" <> "-" <> ToString[inputFileNumbers[[i]]] <> ".txt", thdist];
	Export["S-Tdist-" <> "$size" <> "-" <> ToString[inputFileNumbers[[i]]] <> ".txt", tdist];
	Export["S-Gdist-" <> "$size" <> "-" <> ToString[inputFileNumbers[[i]]] <> ".txt", gdist];
	Export["S-Qdist-" <> "$size" <> "-" <> ToString[inputFileNumbers[[i]]] <> ".txt", qdist];
	thdistplot = ListPlot[thdist];
	tdistplot = ListPlot[tdist];
	gdistplot = ListPlot[gdist];
	qdistplot = ListPlot[qdist];
	Export["S-THdist-" <> "$size" <> "-" <> ToString[inputFileNumbers[[i]]] <> ".jpg", ListPlot[thdist]];
	Export["S-Tdist-" <> "$size" <> "-" <> ToString[inputFileNumbers[[i]]] <> ".jpg", ListPlot[tdist]];
	Export["S-Gdist-" <> "$size" <> "-" <> ToString[inputFileNumbers[[i]]] <> ".jpg", ListPlot[gdist]];
	Export["S-Qdist-" <> "$size" <> "-" <> ToString[inputFileNumbers[[i]]] <> ".jpg", ListPlot[qdist]];,
 {i, 1, Length[inputFileNumbers]}]

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
