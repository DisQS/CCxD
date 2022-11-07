#!/bin/bash
#
# run as: for size in 10 12; do echo $size; ../../PBS/L31-conE-rVal.sh 1.0 $size 1000 0; done

# settings from input


maxRGSteps=${2}
configs=${1}
outputfreq=$3

echo "C-RG: making for configs=" $configs "with maxRGSteps=" $maxRGSteps

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
mkdir raw
mkdir dists

cat > ${wlsfile} << EOD
#!/usr/bin/env wolframscript 
Print["(*Preliminaries*)"];
size=$configs;
maxRGSteps=$maxRGSteps;
tp[t1_, r1_, t2_, r2_, t3_, r3_, t4_, r4_, t5_, 
   r5_, \[Phi]1_, \[Phi]2_, \[Phi]3_, \[Phi]4_] = 
  Abs[(-Exp[I (\[Phi]1 + \[Phi]4 - \[Phi]2)] r1 r3 r5 t2 t4 + 
      Exp[I (\[Phi]1 + \[Phi]4)] t2 t4 - Exp[I \[Phi]4] t1 t3 t4 + 
      t1 t5 + Exp[I \[Phi]3] r2 r3 r4 t1 t5 - 
      Exp[I \[Phi]1] t2 t3 t5)/(-1 - Exp[I \[Phi]3] r2 r3 r4 + 
      Exp[I \[Phi]2] r1 r3 r5 + 
      Exp[I (\[Phi]2 + \[Phi]3)] r1 r2 r4 r5 + 
      Exp[I \[Phi]1] t1 t2 t3 - 
      Exp[I (\[Phi]1 + \[Phi]4)] t1 t2 t4 t5 + 
      Exp[I \[Phi]4] t3 t4 t5)];

Print["Creating initial distribution for t"]
initdist =  
  ProbabilityDistribution[
   Piecewise[{{2 z, 0 < z <= 1}}], {z, -\[Infinity], \[Infinity]}];
last = ParallelTable[RandomVariate[initdist], {i, 1, size}];
Print["Finished distribution creation, moving onto evaluating renormalized t for ", $maxRGSteps, "steps"]
zdata = {};
maxrgsteps = $maxRGSteps;
outputfreq = $outputfreq;
nprint = Floor[maxrgsteps/outputfreq];
Do[
  Print["--- working on rg step: ", ind];
  tpdata = Table[With[{
      x1 = last[[RandomInteger[{1, size}]]],
      x2 = last[[RandomInteger[{1, size}]]],
      x3 = last[[RandomInteger[{1, size}]]],
      x4 = last[[RandomInteger[{1, size}]]],
      x5 = last[[RandomInteger[{1, size}]]]}, 
     tp[x1, Sqrt[1 - x1^2], x2, Sqrt[1 - x2^2], x3, Sqrt[1 - x3^2], 
      x4, Sqrt[1 - x4^2], x5, Sqrt[1 - x5^2], 
      RandomReal[{0, 2 \[Pi]}], RandomReal[{0, 2 \[Pi]}], 
      RandomReal[{0, 2 \[Pi]}], RandomReal[{0, 2 \[Pi]}]]], {i, 1, 
     size}];
	zdata = Log[(1/tpdata)-1];
	tdist = BinCounts[tpdata,{0,1,0.001}];
	gdist = Table[tdist[[i]]/(2*(i/1000)),{i,1,1000}];
	qdist = BinCounts[zdata,{-8,8,0.004}];
	
  If[Mod[ind, nprint] == 0,
   filename = 
    "C-" <> ToString[size] <> "_" <> ToString[ind] <> "raw.nc";
   Print["   saving ", filename];
   Export["$jobdir/raw/" <> filename, tpdata];
	Export["$jobdir/dists/C-Tdist-" <> "$configs" <> "-" <> ToString[ind] <> ".txt", tdist];
        Export["$jobdir/dists/C-Gdist-" <> "$configs" <> "-" <> ToString[ind] <> ".txt", gdist];
        Export["$jobdir/dists/C-Qdist-" <> "$configs" <> "-" <> ToString[ind] <> ".txt", qdist];
        Export["$jobdir/raw/C-Zraw-" <> "$configs" <> "-" <> ToString[ind] <> ".nc", zdata];
   ];
	last = tpdata;,
  {ind, 1, maxrgsteps}];
maindir="$currdir";


Print[maindir];
Print["--- FINISHED!"];
EOD

cat > ${jobfile} << EOD
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
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
