#!/bin/bash




# Parameter initialization with 'default' values
MAXSINGLETHREADSIZE=10000000 # Maximum size allowed in a single process, int
SYMMETRIZE=0 # Boolean flag indicating if the distribution after each step is to be symmetrised
MATRIX=0 # Boolean flag indicating the scattering matrix to be used in calcualtions, 0=Cain, 1=Son (just to find the file)
CONFIGS=500000000 # The amount of values generated to create the distributions (just to find the file)
BINCOUNT=1000 # Number of bins in the distributions 
ZRANGE=25 # Number of bins in the distributions of z, set differently as it is a much larger range
FINDITERATIONS=3 # essentially an indexed derived from the iter number given to the master script, used to pick which FP you want 
ITERATIONS=5 #How many iterations to run for
ACTUALITER=15 # The iteration to start calculating the RG flow from
TASKFARM=0 # Boolean flag indicating whether the process will be run on the taskfarm
READIN=0 # Boolean flag indicating that the calculation is not starting from a standard distribution
BATCHSIZE=20 #e taskfarm at once in a single batch, int
BOOTSTRAPSIZE=1 # Amount of times to do everything in efforts to boostrap further down the line
GEOMETRICDISORDERAMT=0 #Proportion of geometric disorder to include in the calculation [0,0.5]
TRYAGAIN=0;
# Flag and argument capture loop

CURRBATCH=1
CURRPROC=1
CURRITER=0
NEWRG=0
STARTCURRITER=0
MUTABLECONFIGSIZE=0
MUTABLELOOPSIZE=0
PUSHAMT=0
while getopts 'sc:i:p:m:ht:o:r:n:g:b:t:' OPTION; do
        case "$OPTION" in
        t)
		TRYAGAIN=$OPTARG
		;;
	h)
                echo "*****************************************************"
                echo "*            Renormalisation Script Help            *"
                echo "*                                                   *"
                echo "*    -m -> Scattering matrix in use('C' or 'S')     *"
                echo "*    -c -> Number of configs (default 1M)           *"
                echo "*    -i -> Number of RG Steps (default 30)          *"
                echo "*    -r -> Restart process (argument takes iter no) *"
		echo "*    -n -> Number of RG steps for continuing run    *"
                echo "*                                                   *"
                echo "*****************************************************"
                exit 1
                ;;
	g)
		GEOMETRICDISORDERAMT=$OPTARG
		;;
	b)
		BOOTSTRAPSIZE=$OPTARG
		;;
	p)
		PUSHAMT=$OPTARG
		;;
	r)
		READIN=1
		CURRITER=$OPTARG
		STARTCURRITER=$OPTARG
		;;
	n)
		NEWRG=$OPTARG
		;;
	?)
                echo "Invalid input, consult -h for help"
                exit 1
                ;;

        esac
done
shift "$(($OPTIND -1))"


if [[ $CONFIGS -gt $MAXSINGLETHREADSIZE ]]
then
	ONLYZ=1
fi

# Batch and process amount calculation

CONFIGREMAINDER=$(( $CONFIGS % $MAXSINGLETHREADSIZE)) 
echo $CONFIGREMAINDER
NOOFPROCS=$((($CONFIGS- 1 + $MAXSINGLETHREADSIZE)/$MAXSINGLETHREADSIZE))
echo $NOOFPROCS
PROCREMAINDER=$(($NOOFPROCS % $BATCHSIZE))
echo $PROCREMAINDER
NOOFBATCH=$((($NOOFPROCS-1 + $BATCHSIZE)/$BATCHSIZE))
echo $NOOFBATCH



currdir=`pwd`'/../data'
cd $currdir
jobdir="$MATRIX-RG-$CONFIGS-$FINDITERATIONS-$GEOMETRICDISORDERAMT"

generatorjobfile=`printf "$jobdir-generator-script.sh"`
generatorfile=`printf "$jobdir-generator.wls"`
distavgjobfile=`printf "$jobdir-distavg-script.sh"`
distavgfile=`printf "$jobdir-distavg.wls"`
rgjobfile=`printf "$jobdir-rg-script.sh"`
rgfile=`printf "$jobdir-rg.wls"`
symmetrizejobfile=`printf "$jobdir-symmetrize-script.sh"`
symmetrizefile=`printf "$jobdir-symmetrize.wls"`
statsjobfile=`printf "$jobdir-stats-script.wls"`
statsfile=`printf "$jobdir-stats.wls"`

cd $jobdir
if [[ $READIN -eq 1 ]]; then
	ITERATIONS=$NEWRG
fi
for ((j=1; j<=$(($BOOTSTRAPSIZE));j++));
do
	for (( i=$CURRITER;i<=$(($CURRITER+$ITERATIONS));i++ ));
 	do
		mkdir -p push/BS$j/$MATRIX-$CONFIGS-$PUSHAMT/$i  ;
	done;
	if [[ $MATRIX -eq 1 ]]; then
        cp dists/BS$j/$ACTUALITER/$MATRIX-theta-$CONFIGS-averaged-symmetrized.txt push/BS$j/$MATRIX-$CONFIGS-$PUSHAMT/0/$MATRIX-theta-$CONFIGS-averaged.txt
	fi
	cp dists/BS$j/$ACTUALITER/$MATRIX-z-$CONFIGS-averaged-symmetrized.txt push/BS$j/$MATRIX-$CONFIGS-$PUSHAMT/0/$MATRIX-z-$CONFIGS-averaged.txt
	cp dists/BS$j/$ACTUALITER/$MATRIX-t-$CONFIGS-averaged.txt push/BS$j/$MATRIX-$CONFIGS-$PUSHAMT/0/$MATRIX-t-$CONFIGS-averaged.txt
	cp dists/BS$j/$ACTUALITER/$MATRIX-g-$CONFIGS-averaged.txt push/BS$j/$MATRIX-$CONFIGS-$PUSHAMT/0/$MATRIX-g-$CONFIGS-averaged.txt

done;




cat > ${distavgfile} << EOD
#!/usr/bin/env wolframscript 
curriter = ToExpression[\$ScriptCommandLine[[2]]];
pushamt = ToExpression[\$ScriptCommandLine[[3]]];
currstrap = ToExpression[\$ScriptCommandLine[[4]]];
Print["---Averaging distributions for step no ", curriter]
If[$CONFIGREMAINDER === 0, remainder = $MAXSINGLETHREADSIZE, remainder = $CONFIGREMAINDER];
Print["---Importing from "<>"$jobdir/push/BS"<>ToString[currstrap]<>"/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<> ToString[i] <>".txt"];
If[$MATRIX === 1,
	importedth = Table[Flatten[Import["$jobdir/push/BS"<>ToString[currstrap]<>"/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-theta-"<>ToString[i] <>".txt", "Data"]],{i,1,$BATCHSIZE-1}];
	remainderimportth = Import["$jobdir/push/BS"<>ToString[currstrap]<>"/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-theta-"<>"$BATCHSIZE" <>".txt","Data"];
];	

importedt = Table[Flatten[Import["$jobdir/push/BS"<>ToString[currstrap]<>"/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-t-"<> ToString[i] <>".txt","Data"]],{i,1,$BATCHSIZE-1}];
remainderimportt = Import["$jobdir/push/BS"<>ToString[currstrap]<>"/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-t-"<>"$BATCHSIZE" <>".txt","Data"];
	
importedg = Table[Flatten[Import["$jobdir/push/BS"<>ToString[currstrap]<>"/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-g-"<> ToString[i] <>".txt","Data"]],{i,1,$BATCHSIZE-1}];
remainderimportg = Import["$jobdir/push/BS"<>ToString[currstrap]<>"/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-g-"<>"$BATCHSIZE" <>".txt","Data"];
	
importedz = Table[Flatten[Import["$jobdir/push/BS"<>ToString[currstrap]<>"/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<> ToString[i] <>".txt","Data"]],{i,1,$BATCHSIZE-1}];
remainderimportz = Import["$jobdir/push/BS"<>ToString[currstrap]<>"/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>"$BATCHSIZE" <>".txt","Data"];

Print[Length[importedz[[1]]]];
Print[Take[importedz[[1]],1000]];
Print["---Averaging over ", $BATCHSIZE, " distributions"];
If[$MATRIX ===1,
	averagedth = Table[total=0;For[j=1,j<$BATCHSIZE,j++,
			total+=importedth[[j]][[i]]];total = total/$BATCHSIZE; N[total+(remainderimportth[[i]]/($BATCHSIZE * (remainder/$MAXSINGLETHREADSIZE)))]
	,{i,1,$BINCOUNT}];
];
averagedt = Table[total=0;For[j=1,j<$BATCHSIZE,j++,
		total+=importedt[[j]][[i]]];total = total/$BATCHSIZE; N[total+(remainderimportt[[i]]/($BATCHSIZE * (remainder/$MAXSINGLETHREADSIZE)))]
,{i,1,$BINCOUNT}];

averagedg = Table[total=0;For[j=1,j<$BATCHSIZE,j++,
		total+=importedg[[j]][[i]]];total = total/$BATCHSIZE; N[total+(remainderimportg[[i]]/($BATCHSIZE * (remainder/$MAXSINGLETHREADSIZE)))]
,{i,1,$BINCOUNT}];

averagedz = Table[total=0;For[j=1,j<$BATCHSIZE,j++,
		total+=importedz[[j]][[i]]];total = total/$BATCHSIZE; N[total+(remainderimportz[[i]]/($BATCHSIZE * (remainder/$MAXSINGLETHREADSIZE)))]
,{i,1,Length[importedz[[1]]]}];
Print[Length[averagedz]];
Print["---Exporting"]
If[$MATRIX === 1,

	Export["$jobdir/push/BS"<>ToString[currstrap]<>"/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-theta-"<>"$CONFIGS"<>"-averaged.txt",Flatten[averagedth]];
];	
Export["$jobdir/push/BS"<>ToString[currstrap]<>"/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>"$CONFIGS"<>"-averaged.txt",Flatten[averagedz]];
Export["$jobdir/push/BS"<>ToString[currstrap]<>"/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-g-"<>"$CONFIGS"<>"-averaged.txt",Flatten[averagedg]];
Export["$jobdir/push/BS"<>ToString[currstrap]<>"/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-t-"<>"$CONFIGS"<>"-averaged.txt",Flatten[averagedt]];


EOD

cat > ${distavgjobfile} << EOD
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2012
#SBATCH --time=1:00:00
module purge
module load Mathematica
module list
pwd
ls -al $jobdir/$distavgfile
$jobdir/$distavgfile \$CURRITER \$PUSHAMT \$CURRSTRAP



EOD

cat > ${rgfile} << EOD
#!/usr/bin/env wolframscript

curriter = ToExpression[\$ScriptCommandLine[[2]]];
loopsize = ToExpression[\$ScriptCommandLine[[3]]];
configsize = ToExpression[\$ScriptCommandLine[[4]]];
currproc = ToExpression[\$ScriptCommandLine[[5]]];
pushamt = ToExpression[\$ScriptCommandLine[[6]]];
currstrap = ToExpression[\$ScriptCommandLine[[7]]];

p0=$GEOMETRICDISORDERAMT;
p1=$GEOMETRICDISORDERAMT;
Print["Current Iteration: ", curriter];
Print["Current Process: ", currproc];
Print["Disorder Amount: ", p0];
Print["Push Amount: ", pushamt];
Print["Current Strap: ", currstrap];

Launder[sop_, min_, max_, laundersize_] := (
        binWidth = (max - min)/(Length[sop]);
        normed = sop/(Total[sop])*binWidth;
        hmax = Max[normed];
        newData = Table[(
                prospectpoint = {RandomReal[{min, max}], RandomReal[{0, hmax}]};
                binNo = Ceiling[(prospectpoint[[1]] - min)/binWidth];
                While[normed[[binNo]] < prospectpoint[[2]], 
                        prospectpoint = {RandomReal[{min, max}], RandomReal[{0, hmax}]}; 
                        binNo = Ceiling[(prospectpoint[[1]] - min)/binWidth] ];
        prospectpoint[[1]]), laundersize];
        Return[newData]
);

If[$MATRIX === 1,
Print["Test1"];
tp[theta1_, theta2_, theta3_, theta4_, theta5_, phi1_, phi2_, phi3_, phi4_] := 
        Abs[(-Cos[theta1] (Exp[I phi4] Cos[theta3] Cos[theta4] + I Cos[theta5] (Exp[I phi3] Sin[theta2] Sin[theta3] Sin[theta4] - 1)) - Exp[I phi1] Cos[theta2] Cos[theta3] Cos[theta5] - I Exp[I (phi1 + phi4)] Cos[theta4] Cos[theta2] + I Exp[I (phi1 + phi4 - phi2)] Cos[theta2] Cos[theta4] Sin[theta1] Sin[theta3] Sin[theta5])/
        (Exp[I phi4] Cos[theta3] Cos[theta4] Cos[theta5] + Exp[I phi1] Cos[theta1] Cos[theta2] (Cos[theta3] +  I Exp[I phi4] Cos[theta4] Cos[theta5]) + I (-1 + Exp[I phi3] Sin[theta2] Sin[theta3] Sin[theta4] + Exp[I phi2] Sin[theta1] (Sin[theta3] - Exp[I phi3] Sin[theta2] Sin[theta4]) Sin[theta5]))];
,
Print["Test2"];
tp[t1_, r1_, t2_, r2_, t3_, r3_, t4_, r4_, t5_, r5_, phi1_, phi2_, phi3_, phi4_] := 
        Abs[(-Exp[I (phi1 + phi4 - phi2)] r1 r3 r5 t2 t4 + Exp[I (phi1 + phi4)] t2 t4 - Exp[I phi4] t1 t3 t4 + t1 t5 + Exp[I phi3] r2 r3 r4 t1 t5 - Exp[I phi1] t2 t3 t5)/
        (-1 - Exp[I phi3] r2 r3 r4 + Exp[I phi2] r1 r3 r5 + Exp[I (phi2 + phi3)] r1 r2 r4 r5 + Exp[I phi1] t1 t2 t3 - Exp[I (phi1 + phi4)] t1 t2 t4 t5 + Exp[I phi4] t3 t4 t5)];
];
Print["Importing from -- "<>"$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/BS"<>ToString[currstrap]<>"/"<>ToString[curriter-1]<>"/"<>"$MATRIX"<>"-z-"<>"$CONFIGS"<>"-averaged.txt"];
import = Flatten[Import["$jobdir/push/BS"<>ToString[currstrap]<>"/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter-1]<>"/"<>"$MATRIX"<>"-z-"<>"$CONFIGS"<>"-averaged.txt","Data"]];
min=pushamt-$ZRANGE;
max=pushamt+$ZRANGE;

generatedfromdist=Launder[import,min,max,configsize];
Print[import[[1]]];
generatedfromdist=Sqrt[1/(E^generatedfromdist+1)];
If[$MATRIX===1,
	generatedfromdist=ArcCos[generatedfromdist];
];


Print["---Generated distribution to perform RG step upon"]
finalt = Table[0,{i,1,$BINCOUNT}];
finalg = Table[0,{i,1,$BINCOUNT}];
finalz = Table[0,{i,1,2 * $ZRANGE * $BINCOUNT}];
finalth = Table[0,{i,1,$BINCOUNT}];

If[$MATRIX===1,       
	For[i=0,i<loopsize,i++,
		Print["---RG step "<>ToString[(i/loopsize) * 100]<>"% done"];
		rgstep = Table[With[{
                	t1 = generatedfromdist[[RandomInteger[{1, configsize}]]],
                	t2 = generatedfromdist[[RandomInteger[{1, configsize}]]],
                	t3 = generatedfromdist[[RandomInteger[{1, configsize}]]],
                	t4 = generatedfromdist[[RandomInteger[{1, configsize}]]],
                	t5 = generatedfromdist[[RandomInteger[{1, configsize}]]]},
                	ArcCos[tp[t1, t2, t3,t4, t5, RandomReal[{0, 2 Pi}], RandomReal[{0, 2 Pi}], RandomReal[{0, 2 Pi}], RandomReal[{0, 2 Pi}]]]],
                {j, 1, configsize}];
		Print["RG step complete"];
		
		tvalues = Cos[rgstep];
		gvalues = tvalues^2;
		zvalues = Log[(1/tvalues^2)-1];
		finalt+=BinCounts[tvalues,{0,1,1/$BINCOUNT}]/loopsize;
                finalg+=BinCounts[gvalues,{0,1,1/$BINCOUNT}]/loopsize;
                finalz+=BinCounts[zvalues,{-$ZRANGE,$ZRANGE,1/$BINCOUNT}]/loopsize;
                finalth+=BinCounts[rgstep,{0,Pi/2,(Pi/2)/$BINCOUNT}]/loopsize;

	];
	filename = "$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/BS"<>ToString[currstrap]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-theta-"<>ToString[currproc]<>".nc";
        Print["---Exporting to : ", filename];

	Export["$jobdir/push/BS"<>ToString[currstrap]<>"/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-theta-"<>ToString[currproc]<>".txt",N[finalth]];
        Export["$jobdir/push/BS"<>ToString[currstrap]<>"/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-t-"<>ToString[currproc]<>".txt",N[finalt]];
        Export["$jobdir/push/BS"<>ToString[currstrap]<>"/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-g-"<>ToString[currproc]<>".txt",N[finalg]];
	Export["$jobdir/push/BS"<>ToString[currstrap]<>"/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>ToString[currproc]<>".txt",N[finalz]];

	Print["---RG step 100% done"];
	
,	
        For[i=0,i<loopsize,i++,
		Print["---RG step "<>ToString[Round[(i/loopsize)*100]]<>"% done"];
                rgstep = Table[With[{
                t1=With[{test = RandomReal[]}, If[test < p0, 0.001, If[test < (p0 + p1), 0.999,generatedfromdist[[RandomInteger[{1, configsize}]]]]]],
                        t2=With[{test = RandomReal[]}, If[test < p0, 0.001, If[test < (p0 + p1), 0.999,generatedfromdist[[RandomInteger[{1, configsize}]]]]]],
                        t3=With[{test = RandomReal[]}, If[test < p0, 0.001, If[test < (p0 + p1), 0.999,generatedfromdist[[RandomInteger[{1, configsize}]]]]]],
                        t4=With[{test = RandomReal[]}, If[test < p0, 0.001, If[test < (p0 + p1), 0.999,generatedfromdist[[RandomInteger[{1, configsize}]]]]]],
                        t5=With[{test = RandomReal[]}, If[test < p0, 0.001, If[test < (p0 + p1), 0.999,generatedfromdist[[RandomInteger[{1, configsize}]]]]]]},
        
		tp[t1, Sqrt[1 - t1^2], t2, Sqrt[1 - t2^2], t3, Sqrt[1 - t3^2], t4, Sqrt[1 - t4^2], t5, Sqrt[1 - t5^2], RandomReal[{0, 2 Pi}], RandomReal[{0, 2 Pi}], RandomReal[{0, 2 Pi}], RandomReal[{0, 2 Pi}]]], 
                {i, 1, configsize}];
		gvalues = rgstep^2;
		zvalues = Log[(1/rgstep^2)-1];
                finalt+=(BinCounts[rgstep,{0,1,1/$BINCOUNT}]/loopsize);
                finalg+=(BinCounts[gvalues,{0,1,1/$BINCOUNT}]/loopsize);
                finalz+=(BinCounts[zvalues,{-$ZRANGE,$ZRANGE,1/$BINCOUNT}]/loopsize);
		Print[finalz[[1]]];

	];
        filename = "$jobdir/push/BS"<>ToString[currstrap]<>"/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>ToString[currproc]<>".txt";
	Print["---Exporting to: ", filename];
        Export["$jobdir/push/BS"<>ToString[currstrap]<>"/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-t-"<>ToString[currproc]<>".txt",N[finalt]];
        Export["$jobdir/push/BS"<>ToString[currstrap]<>"/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-g-"<>ToString[currproc]<>".txt",N[finalg]];
        Export["$jobdir/push/BS"<>ToString[currstrap]<>"/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>ToString[currproc]<>".txt",N[finalz]];

	
	
	Print["---RG step 100% done"];
];






EOD

cat > ${rgjobfile} << EOD
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2012
#SBATCH --time=10:00:00
module purge
module load Mathematica
module list
pwd
ls -al $jobdir/$rgfile
$jobdir/$rgfile \$CURRITER \$MUTABLELOOPSIZE \$MUTABLECONFIGSIZE \$CURRPROC \$PUSHAMT \$CURRSTRAP




EOD


cd ..

chmod 755 ${jobdir}/${generatorfile}
chmod 755 ${jobdir}/${generatorjobfile}
chmod 755 ${jobdir}/${distavgfile}
chmod 755 ${jobdir}/${distavgjobfile}
chmod 755 ${jobdir}/${rgfile}
chmod 755 ${jobdir}/${rgjobfile}
chmod 755 ${jobdir}/${symmetrizefile}
chmod 755 ${jobdir}/${symmetrizejobfile}
chmod 755 ${jobdir}/${statsfile}
chmod 755 ${jobdir}/${statsjobfile}
CURRSTRAP=0;
for ((z=1;z<=$BOOTSTRAPSIZE; z++ ));
do
CURRSTRAP=$(($CURRSTRAP+1));
if [[ $TRYAGAIN -ge 1 ]]; then
	CURRSTRAP=$(($TRYAGAIN));
fi
CURRITER=$STARTCURRITER;
if [[ $READIN -eq 1 ]]; then
CURRITER=$(($CURRITER-1));
fi

endproc=$(sbatch --wrap="echo $CURRITER" | cut -d ' ' -f4);
for  (( i=1;i<=$ITERATIONS;i++ ));
do
	CURRPROC=1
	CURRITER=$(($CURRITER+1));
	joined=""
	delim=""
	if [[ $NOOFPROCS -ge $BATCHSIZE ]]; then
		for  (( j=1;j<=$BATCHSIZE;j++ ));
		do
			jid=$(sbatch --dependency=afterok:$endproc --export=CURRITER=$CURRITER,MUTABLELOOPSIZE=$(($NOOFBATCH)),MUTABLECONFIGSIZE=$MAXSINGLETHREADSIZE,CURRPROC=$CURRPROC,PUSHAMT=$PUSHAMT,CURRSTRAP=$CURRSTRAP ${jobdir}/${rgjobfile} | cut -d ' ' -f4);
			joined="$joined$delim$jid";
			delim=","
			CURRPROC=$(($CURRPROC+1));
		done;
	fi
	
	CURRPROC=$(( ($BATCHSIZE * $NOOFBATCH) + 1 ))
	for  (( j=2;j<=$PROCREMAINDER;j++ ));
	do
		jid=$(sbatch --dependency=afterok:$endproc --export=CURRITER=$CURRITER,MUTABLELOOPSIZE=$((1)),MUTABLECONFIGSIZE=$MAXSINGLETHREADSIZE,CURRPROC=$CURRPROC,PUSHAMT=$PUSHAMT,CURRSTRAP=$CURRSTRAP ${jobdir}/${rgjobfile} | cut -d ' ' -f4);
		joined="$joined$delim$jid";
		CURRPROC=$(($CURRPROC + 1));
	done;
	if [[ $CONFIGREMAINDER -ge 1 ]]; 
	then
		jid=$(sbatch --dependency=afterok:$endproc --export=CURRITER=$CURRITER,MUTABLELOOPSIZE=$((1)),MUTABLECONFIGSIZE=$CONFIGREMAINDER,CURRPROC=$CURRPROC,PUSHAMT=$PUSHAMT,CURRSTRAP=$CURRSTRAP ${jobdir}/${rgjobfile} | cut -d ' ' -f4);
		joined="$joined$delim$jid";
	fi

	endproc=$(sbatch --dependency=afterok:$joined --export=CURRITER=$CURRITER,PUSHAMT=$PUSHAMT,CURRSTRAP=$CURRSTRAP ${jobdir}/${distavgjobfile} | cut -d ' ' -f4);
	
done;
done;





