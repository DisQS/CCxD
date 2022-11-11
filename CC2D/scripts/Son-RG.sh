#!/bin/bash
#
# run as: for size in 10 12; do echo $size; ../../PBS/L31-conE-rVal.sh 1.0 $size 1000 0; done

# settings from input


maxRGSteps=${2}
configs=${1}
outputfreq=$3

echo "S-RG: making for configs=" $configs "with maxRGSteps=" $maxRGSteps

# settings for directories

currdir=`pwd`'/../data'
cd $currdir
jobdir="S-RG-$configs"
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
tp1[\[Theta]1_, \[Theta]2_, \[Theta]3_, \[Theta]4_, \[Theta]5_, \
\[Phi]1_, \[Phi]2_, \[Phi]3_, \[Phi]4_] := 
 Abs[(-Cos[\[Theta]1] (Exp[I \[Phi]4] Cos[\[Theta]3] Cos[\[Theta]4] + 
        I Cos[\[Theta]5] (Exp[
             I \[Phi]3] Sin[\[Theta]2] Sin[\[Theta]3] Sin[\[Theta]4] \
- 1)) - Exp[I \[Phi]1] Cos[\[Theta]2] Cos[\[Theta]3] Cos[\[Theta]5] - 
     I Exp[I (\[Phi]1 + \[Phi]4)] Cos[\[Theta]4] Cos[\[Theta]2] + 
     I Exp[I (\[Phi]1 + \[Phi]4 - \[Phi]2)] Cos[\[Theta]2] Cos[\
\[Theta]4] Sin[\[Theta]1] Sin[\[Theta]3] Sin[\[Theta]5])/(Exp[
       I \[Phi]4] Cos[\[Theta]3] Cos[\[Theta]4] Cos[\[Theta]5] + 
     Exp[I \[Phi]1] Cos[\[Theta]1] Cos[\[Theta]2] (Cos[\[Theta]3] +  
        I Exp[I \[Phi]4] Cos[\[Theta]4] Cos[\[Theta]5]) + 
     I (-1 + Exp[
          I \[Phi]3] Sin[\[Theta]2] Sin[\[Theta]3] Sin[\[Theta]4] + 
        Exp[I \[Phi]2] Sin[\[Theta]1] (Sin[\[Theta]3] - 
           Exp[I \[Phi]3] Sin[\[Theta]2] Sin[\[Theta]4]) \
Sin[\[Theta]5]))]


Print["Creating initial distribution for t"]
initdist =  
  ProbabilityDistribution[
   Piecewise[{{2 z, 0 < z <= \[Pi]/2}}], {z, -\[Infinity], \[Infinity]}];
last = ParallelTable[RandomReal[{0,\[Pi]/2}], {i, 1, size}];
Print["Finished distribution creation, moving onto evaluating renormalized t for ", $maxRGSteps, "steps"]
maxrgsteps = $maxRGSteps;
outputfreq = $outputfreq;
nprint = Floor[maxrgsteps/outputfreq];
Do[
  Print["--- working on rg step: ", ind];
  thpdata = Table[With[{
     	\[Theta]1 = last[[RandomInteger[{1, size}]]], \[Theta]2 = 
   last[[RandomInteger[{1, size}]]], \[Theta]3 = 
   last[[RandomInteger[{1, size}]]], \[Theta]4 = 
   last[[RandomInteger[{1, size}]]], \[Theta]5 = 
   last[[RandomInteger[{1, size}]]]}, 
 ArcCos[tp1[\[Theta]1, \[Theta]2, \[Theta]3, \[Theta]4, \[Theta]5, 
   RandomReal[{0, 2 \[Pi]}], RandomReal[{0, 2 \[Pi]}], 
   RandomReal[{0, 2 \[Pi]}], RandomReal[{0, 2 \[Pi]}]]]], {i, 1, size}];
  If[Mod[ind, nprint] == 0,
	tdata = Cos[thpdata];
	gdata = tdata^2;
	zdata = Log[(1/tdata^2)-1];
	thdist = BinCounts[thpdata,{0,\[Pi]/2,0.001}];
	tdist = BinCounts[tdata,{0,1,0.001}];
	gdist = BinCounts[gdata,{0,1,0.001}];
	qdist = BinCounts[zdata,{-8,8,0.004}];
   filename =  "S-" <> ToString[size] <> "_" <> ToString[ind] <> "raw.nc";
   Print["   saving ", filename];
   Export["$jobdir/raw/" <> filename, thpdata];
	        Export["$jobdir/raw/S-Traw-"<>"$configs"<>"-"<>ToString[ind]<>".nc",tdata];
        Export["$jobdir/raw/S-Graw-"<>"$configs"<>"-"<>ToString[ind]<>".nc",gdata];
        Export["$jobdir/raw/S-Zraw-"<>"$configs"<>"-"<>ToString[ind]<>".nc",zdata];
        Export["$jobdir/dists/S-Thdist-"<>"$configs"<>"-"<>ToString[ind]<>".txt",thdist];
        Export["$jobdir/dists/S-Tdist-"<>"$configs"<>"-"<>ToString[ind]<>".txt",tdist];
        Export["$jobdir/dists/S-Gdist-"<>"$configs"<>"-"<>ToString[ind]<>".txt",gdist];
        Export["$jobdir/dists/S-Qdist-"<>"$configs"<>"-"<>ToString[ind]<>".txt",qdist];
   ];
	last = thpdata,
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
