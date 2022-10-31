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
initlist = ParallelTable[RandomVariate[initdist], {i, 1, size}];
Print["Finished distribution creation, moving onto evaluating renormalized t for ", $maxRGSteps, "steps"]
thpdata = {initlist};
maxrgsteps = $maxRGSteps;
outputfreq = $outputfreq;
nprint = Floor[maxrgsteps/outputfreq];
Do[
  Print["--- working on rg step: ", ind];
  AppendTo[thpdata, Table[With[{
     	\[Theta]1 = thpdata[[ind]][[RandomInteger[{1, size}]]], \[Theta]2 = 
   thpdata[[ind]][[RandomInteger[{1, size}]]], \[Theta]3 = 
   thpdata[[ind]][[RandomInteger[{1, size}]]], \[Theta]4 = 
   thpdata[[ind]][[RandomInteger[{1, size}]]], \[Theta]5 = 
   thpdata[[ind]][[RandomInteger[{1, size}]]]}, 
 ArcCos[tp1[\[Theta]1, \[Theta]2, \[Theta]3, \[Theta]4, \[Theta]5, 
   RandomReal[{0, 2 \[Pi]}], RandomReal[{0, 2 \[Pi]}], 
   RandomReal[{0, 2 \[Pi]}], RandomReal[{0, 2 \[Pi]}]]]], {i, 1, size}]];
  If[Mod[ind, nprint] == 0,
   filename = 
    "S-" <> ToString[size] <> "_" <> ToString[ind] <> "raw.nc";
   Print["   saving ", filename];
   Export["$jobdir/" <> filename, thpdata[[ind]]];
   ],
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
#sbatch ${jobdir}/${jobfile} # for queueing system
(source ${jobdir}/${jobfile} ) >& ${jobdir}/${logfile} & # for parallel shell execution

#echo "<return>"
sleep 1

cd ..
