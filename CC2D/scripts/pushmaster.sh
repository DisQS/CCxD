#!/bin/bash




# Parameter initialization with 'default' values
MAXSINGLETHREADSIZE=10000000 # Maximum size allowed in a single process, int
SYMMETRIZE=0 # Boolean flag indicating if the distribution after each step is to be symmetrised
MATRIX=0 # Boolean flag indicating the scattering matrix to be used in calcualtions, 0=Cain, 1=Son (just to find the file)
CONFIGS=500000000 # The amount of values generated to create the distributions (just to find the file)
BINCOUNT=1000 # Number of bins in the distributions 
ZBINCOUNT=8000 # Number of bins in the distributions of z, set differently as it is a much larger range
FINDITERATIONS=10 # The amount of RG steps 
ITERATIONS=10
ACTUALITER=13 # The iteration to start calculating the RG flow from
TASKFARM=0 # Boolean flag indicating whether the process will be run on the taskfarm
READIN=0 # Boolean flag indicating that the calculation is not starting from a standard distribution
BATCHSIZE=5 #e taskfarm at once in a single batch, int
# Flag and argument capture loop
CURRBATCH=1
CURRPROC=1
CURRITER=0
NEWRG=0
MUTABLECONFIGSIZE=0
MUTABLELOOPSIZE=0
PUSHAMT=0
while getopts 'sc:i:p:m:ht:o:r:n:' OPTION; do
        case "$OPTION" in
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
	p)
		PUSHAMT=$OPTARG
		;;
	r)
		READIN=1
		CURRITER=$OPTARG
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
jobdir="$MATRIX-RG-$CONFIGS-$FINDITERATIONS"

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
for (( i=$CURRITER;i<=$(($CURRITER+$ITERATIONS));i++ ));
 do
	mkdir -p push/$MATRIX-$CONFIGS-$PUSHAMT/$i  ;
done
if [[ $MATRIX -eq 1 ]]; then
	cp dists/$ACTUALITER/$MATRIX-theta-$CONFIGS-averaged-symmetrized.txt push/$MATRIX-$CONFIGS-$PUSHAMT/0/$MATRIX-theta-$CONFIGS-averaged.txt
fi
cp dists/$ACTUALITER/$MATRIX-z-$CONFIGS-averaged-symmetrized.txt push/$MATRIX-$CONFIGS-$PUSHAMT/0/$MATRIX-z-$CONFIGS-averaged.txt
cp dists/$ACTUALITER/$MATRIX-t-$CONFIGS-averaged.txt push/$MATRIX-$CONFIGS-$PUSHAMT/0/$MATRIX-t-$CONFIGS-averaged.txt
cp dists/$ACTUALITER/$MATRIX-g-$CONFIGS-averaged.txt push/$MATRIX-$CONFIGS-$PUSHAMT/0/$MATRIX-g-$CONFIGS-averaged.txt


# Each instance of generator file generates a maxsinglethreadsize amount batchsize amount of times and saves them in separate files. This is to reduce overhead.
cat > ${generatorfile} << EOD
#!/usr/bin/env wolframscript
(*Importing mutable variables from the master job file*)
loopsize = ToExpression[\$ScriptCommandLine[[2]]];
configsize = ToExpression[\$ScriptCommandLine[[3]]];
currproc = ToExpression[\$ScriptCommandLine[[4]]];
Print["---Generating initial distribution for process number ",currproc];
If[$MATRIX === 1,
                Print["Creating initial distribution for theta"];
		For[i=0,i<loopsize,i++,
			Print["---Generating step "<>ToString[Round[(i/loopsize) * 100]]<>"% done"];
			(*Currently generated theta uniformly, at this point changes to the initial distribution will have to be somewhat manual, although a manually created distribution could be supplied as a 0th iteration*)
                	generate = Table[RandomReal[{0,Pi/2}],{i,1,configsize}];
			tvalues = Cos[generate];
			zvalues = Log[(1/tvalues^2)-1];
			gvalues = tvalues^2;
			Print["---Exporting to: "<> "$jobdir/raw/0/"<>"$MATRIX"<>"-theta-"<>"$CONFIGS"<>"-"<> ToString[(($BATCHSIZE * i)+currproc)] <>".nc"];
			Export["$jobdir/dists/0/"<>"$MATRIX"<>"-theta-"<>ToString[(($BATCHSIZE* i)+ currproc)] <>".txt",BinCounts[generate,{0,Pi/2,(Pi/2)/$BINCOUNT}]];
			Export["$jobdir/dists/0/"<>"$MATRIX"<>"-t-"<>ToString[(($BATCHSIZE* i)+ currproc)] <>".txt",BinCounts[tvalues,{0,1,1/$BINCOUNT}]];
			Export["$jobdir/dists/0/"<>"$MATRIX"<>"-g-"<>ToString[(($BATCHSIZE* i)+ currproc)] <>".txt",BinCounts[gvalues,{0,1,1/$BINCOUNT}]];
			Export["$jobdir/dists/0/"<>"$MATRIX"<>"-z-"<>ToString[(($BATCHSIZE* i)+ currproc)] <>".txt",BinCounts[zvalues,{-8,8,1/$BINCOUNT}]];
		];

,
                Print["Creating initial distribution for t"];
		For[i=0,i<loopsize,i++,
			Print["---Generating step "<>ToString[Round[(i/loopsize)*100]]<>"% done"];
                	generate = Table[Sqrt[RandomReal[]],{i,1,configsize}];
			gvalues = generate^2;
			zvalues = Log[(1/generate^2)-1];
			Print["---Exporting to: "<>"$jobdir/raw/0/"<>"$MATRIX"<>"-t-"<>ToString[(($BATCHSIZE * i)+ currproc)]<>".nc"];
                	Export["$jobdir/dists/0/"<>"$MATRIX"<>"-t-"<>ToString[(($BATCHSIZE * i) + currproc)]<>".txt",BinCounts[generate,{0,1,1/$BINCOUNT}]];
                	Export["$jobdir/dists/0/"<>"$MATRIX"<>"-g-"<>ToString[(($BATCHSIZE * i) + currproc)]<>".txt",BinCounts[gvalues,{0,1,1/$BINCOUNT}]];
                	Export["$jobdir/dists/0/"<>"$MATRIX"<>"-z-"<>ToString[(($BATCHSIZE * i) + currproc)]<>".txt",BinCounts[zvalues,{(-8 - 1/$BINCOUNT),8,1/$BINCOUNT}]];
		];
];
Print["---Generating step 100% done"]
EOD
MUTABLELOOPSIZE=5
cat > ${generatorjobfile} << EOD
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2012
#SBATCH --time=1:00:00
module purge
module load Mathematica
module list
pwd
echo $1
echo $2
echo $3
ls -al $jobdir/$generatorfile
$jobdir/$generatorfile \$MUTABLELOOPSIZE \$MUTABLECONFIGSIZE \$CURRPROC




EOD


cat > ${distavgfile} << EOD
#!/usr/bin/env wolframscript 
curriter = ToExpression[\$ScriptCommandLine[[2]]];
pushamt = ToExpression[\$ScriptCommandLine[[3]]];
Print["---Averaging distributions for step no ", curriter]
If[$CONFIGREMAINDER === 0, remainder = $MAXSINGLETHREADSIZE, remainder = $CONFIGREMAINDER];
Print["---Importing from "<>"$jobdir/dists/"<>ToString[curriter]<>"/"];
If[$MATRIX === 1,
	importedth = Table[Flatten[Import["$jobdir/push/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-theta-"<>ToString[i] <>".txt", "Data"]],{i,1,$NOOFPROCS-1}];
	remainderimportth = Import["$jobdir/push/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-theta-"<>"$NOOFPROCS" <>".txt","Data"];
	
	importedt = Table[Flatten[Import["$jobdir/push/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-t-"<> ToString[i] <>".txt", "Data"]],{i,1,$NOOFPROCS-1}];
	remainderimportt = Import["$jobdir/push/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-t-"<>"$NOOFPROCS" <>".txt","Data"];
	
	importedg = Table[Flatten[Import["$jobdir/push/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-g-"<> ToString[i] <>".txt", "Data"]],{i,1,$NOOFPROCS-1}];
	remainderimportg = Import["$jobdir/push/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-g-"<>"$NOOFPROCS" <>".txt","Data"];
	
	importedz = Table[Flatten[Import["$jobdir/push/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<> ToString[i] <>".txt", "Data"]],{i,1,$NOOFPROCS-1}];
	remainderimportz = Import["$jobdir/push/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>"$NOOFPROCS" <>".txt","Data"];
,
	importedt = Table[Flatten[Import["$jobdir/push/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-t-"<> ToString[i] <>".txt","Data"]],{i,1,$NOOFPROCS-1}];
	remainderimportt = Import["$jobdir/push/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-t-"<>"$NOOFPROCS" <>".txt","Data"];
	
	importedg = Table[Flatten[Import["$jobdir/push/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-g-"<> ToString[i] <>".txt","Data"]],{i,1,$NOOFPROCS-1}];
	remainderimportg = Import["$jobdir/push/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-g-"<>"$NOOFPROCS" <>".txt","Data"];
	
	importedz = Table[Flatten[Import["$jobdir/push/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<> ToString[i] <>".txt","Data"]],{i,1,$NOOFPROCS-1}];
	remainderimportz = Import["$jobdir/push/"<>"$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>"$NOOFPROCS" <>".txt","Data"];
];
Print[Length[importedz[[1]]]];
Print[Take[importedz[[1]],1000]];
Print["---Averaging over ", $NOOFPROCS, " distributions"];
If[$MATRIX ===1,
	averagedth = Table[total=0;For[j=1,j<$NOOFPROCS,j++,
			total+=importedth[[j]][[i]]];total = total/$NOOFPROCS; N[total+(remainderimportth[[i]]/($NOOFPROCS * (remainder/$MAXSINGLETHREADSIZE)))]
	,{i,1,$BINCOUNT}];
];
averagedt = Table[total=0;For[j=1,j<$NOOFPROCS,j++,
		total+=importedt[[j]][[i]]];total = total/$NOOFPROCS; N[total+(remainderimportt[[i]]/($NOOFPROCS * (remainder/$MAXSINGLETHREADSIZE)))]
,{i,1,$BINCOUNT}];

averagedg = Table[total=0;For[j=1,j<$NOOFPROCS,j++,
		total+=importedg[[j]][[i]]];total = total/$NOOFPROCS; N[total+(remainderimportg[[i]]/($NOOFPROCS * (remainder/$MAXSINGLETHREADSIZE)))]
,{i,1,$BINCOUNT}];

averagedz = Table[total=0;For[j=1,j<$NOOFPROCS,j++,
		total+=importedz[[j]][[i]]];total = total/$NOOFPROCS; N[total+(remainderimportz[[i]]/($NOOFPROCS * (remainder/$MAXSINGLETHREADSIZE)))]
,{i,1,Length[importedz[[1]]]}];
Print[Length[averagedz]];
Print["---Exporting"]
If[$MATRIX === 1,

	Export["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-theta-"<>"$CONFIGS"<>"-averaged.txt",Flatten[averagedth]];
	Export["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-t-"<>"$CONFIGS"<>"-averaged.txt",Flatten[averagedt]];
	Export["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-g-"<>"$CONFIGS"<>"-averaged.txt",Flatten[averagedg]];
	Export["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>"$CONFIGS"<>"-averaged.txt",Flatten[averagedz]];
,
	Export["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>"$CONFIGS"<>"-averaged.txt",Flatten[averagedz]];
	Export["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-g-"<>"$CONFIGS"<>"-averaged.txt",Flatten[averagedg]];
	Export["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-t-"<>"$CONFIGS"<>"-averaged.txt",Flatten[averagedt]];
];

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
$jobdir/$distavgfile \$CURRITER \$PUSHAMT



EOD

cat > ${rgfile} << EOD
#!/usr/bin/env wolframscript

curriter = ToExpression[\$ScriptCommandLine[[2]]];
loopsize = ToExpression[\$ScriptCommandLine[[3]]];
configsize = ToExpression[\$ScriptCommandLine[[4]]];
currproc = ToExpression[\$ScriptCommandLine[[5]]];
pushamt = ToExpression[\$ScriptCommandLine[[6]]];
Print["Current Iteration: ", curriter];
Print["Current Process: ", currproc];

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
max = 1;
min = 0;

If[$MATRIX === 1,
	max = Pi/2;
,
	max = 1;];
If[$MATRIX === 1, 
	import = Flatten[Import["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter-1]<>"/"<>"$MATRIX"<>"-theta-"<>"$CONFIGS"<>"-averaged.txt","Data"]]
,
	import = Flatten[Import["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter-1]<>"/"<>"$MATRIX"<>"-t-"<>"$CONFIGS"<>"-averaged.txt","Data"]]
];
If[curriter===1,
	import = Flatten[Import["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter-1]<>"/"<>"$MATRIX"<>"-z-"<>"$CONFIGS"<>"-averaged.txt","Data"]];
	min=pushamt-8;
	max=pushamt+8;
,
];

generatedfromdist=Launder[import,min,max,configsize];
If[curriter===1,
	generatedfromdist=Sqrt[1/(E^generatedfromdist+1)];
	If[$MATRIX===1,
		generatedfromdist=ArcCos[generatedfromdist];
	];
];


Print["---Generated distribution to perform RG step upon"]

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
		filename = "$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-theta-"<>ToString[(($BATCHSIZE * i)+ currproc)]<>".nc";
		Print["---Exporting to : ", filename];
                Export["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-theta-"<>ToString[(($BATCHSIZE * i) + currproc)]<>".txt",BinCounts[rgstep,{0,Pi/2,(Pi/2)/$BINCOUNT}]];
                Export["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-t-"<>ToString[(($BATCHSIZE * i) + currproc)]<>".txt",BinCounts[tvalues,{0,1,1/$BINCOUNT}]];
                Export["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-g-"<>ToString[(($BATCHSIZE * i) + currproc)]<>".txt",BinCounts[gvalues,{0,1,1/$BINCOUNT}]];
                Export["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>ToString[(($BATCHSIZE * i) + currproc)]<>".txt",BinCounts[zvalues,{-12,12,1/$BINCOUNT}]];

	];
	Print["---RG step 100% done"];
	
,	
        For[i=0,i<loopsize,i++,
		Print["---RG step "<>ToString[Round[(i/loopsize)*100]]<>"% done"];
                rgstep = Table[With[{
                        t1 = generatedfromdist[[RandomInteger[{1, configsize}]]],
                        t2 = generatedfromdist[[RandomInteger[{1, configsize}]]],
                        t3 = generatedfromdist[[RandomInteger[{1, configsize}]]],
                        t4 = generatedfromdist[[RandomInteger[{1, configsize}]]],
                        t5 = generatedfromdist[[RandomInteger[{1, configsize}]]]},
                        tp[t1, Sqrt[1 - t1^2], t2, Sqrt[1 - t2^2], t3, Sqrt[1 - t3^2], t4, Sqrt[1 - t4^2], t5, Sqrt[1 - t5^2], RandomReal[{0, 2 Pi}], RandomReal[{0, 2 Pi}], RandomReal[{0, 2 Pi}], RandomReal[{0, 2 Pi}]]], 
                {i, 1, configsize}];
		gvalues = rgstep^2;
		zvalues = Log[(1/rgstep^2)-1];
		filename = "$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-t-"<>ToString[(($BATCHSIZE * i)+ currproc)]<>".nc";
		Print["---Exporting to: ", filename];
                Export["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-t-"<>ToString[(($BATCHSIZE * i) + currproc)]<>".txt",BinCounts[rgstep,{0,1,1/$BINCOUNT}]];
                Export["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-g-"<>ToString[(($BATCHSIZE * i) + currproc)]<>".txt",BinCounts[gvalues,{0,1,1/$BINCOUNT}]];
                Export["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>ToString[(($BATCHSIZE * i) + currproc)]<>".txt",BinCounts[zvalues,{-12,12,1/$BINCOUNT}]];

	
	];
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
$jobdir/$rgfile \$CURRITER \$MUTABLELOOPSIZE \$MUTABLECONFIGSIZE \$CURRPROC \$PUSHAMT




EOD


cat > ${symmetrizefile} << EOD
#!/usr/bin/env wolframscript 
curriter = ToExpression[\$ScriptCommandLine[[2]]];
Print["---Current iteration: "<>ToString[curriter]];
Print["---Symmetrising z distribution for current iterations"]
importedz = Flatten[Import["$jobdir/dists/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>"$CONFIGS"<>"-averaged.txt", "Data"]];
If[$MATRIX === 1,
importedtheta = Flatten[Import["$jobdir/dists/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-theta-"<>"$CONFIGS"<>"-averaged.txt", "Data"]];
thetasym = ((importedtheta + Reverse[importedtheta])/2);
Export["$jobdir/dists/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-theta-"<>"$CONFIGS"<>"-averaged-symmetrized.txt",thetasym];
];
zsym = ((importedz + Reverse[importedz])/2);
Print["---Exporting symmetrised distribution to: "<>"$jobdir/dists/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>"$CONFIGS"<>"-averaged-symmetrized.txt"];
Export["$jobdir/dists/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>"$CONFIGS"<>"-averaged-symmetrized.txt",zsym];
Print["---Done!"]
EOD

cat > ${symmetrizejobfile} << EOD
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2012
#SBATCH --time=1:00:00
module purge
module load Mathematica
module list
pwd
ls -al $jobdir/$symmetrizefile
$jobdir/$symmetrizefile \$CURRITER





EOD

cat > ${statsfile} << EOD
#!/usr/bin/env wolframscript
curriter = ToExpression[\$ScriptCommandLine[[2]]];
pushamt = ToExpression[\$ScriptCommandLine[[3]]];
If[$SYMMETRIZE===1,
	zdist = Flatten[Import["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>"$CONFIGS"<>"-averaged-symmetrized.txt","Data"]];
	prevzdist = Flatten[Import["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter-1]<>"/"<>"$MATRIX"<>"-z-"<>"$CONFIGS"<>"-averaged-symmetrized.txt","Data"]];
,
	zdist = Flatten[Import["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/"<>"$MATRIX"<>"-z-"<>"$CONFIGS"<>"-averaged.txt","Data"]];
	prevzdist = Flatten[Import["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter-1]<>"/"<>"$MATRIX"<>"-z-"<>"$CONFIGS"<>"-averaged.txt","Data"]];
];
normed = zdist*$BINCOUNT/Total[zdist];
normedprev = prevzdist*$BINCOUNT/Total[prevzdist];
centeredzdist = Table[{(i/1000)-12.0005,normed[[i]]},{i,1,24000}];
top = Take[Ordering[normed],-(Length[normed]/20)];
tip = Table[{(i/1000) - 12.0005, normed[[i]]}, {i, top}];
model[x_] = ampl Evaluate[PDF[NormalDistribution[x0, sigma], x]];
fit = FindFit[tip, model[x], {ampl, x0, sigma}, x];
zmax = x0/.fit[[2]];
centeredprev = Table[{(i/1000)-12.0005,normedprev[[i]]},{i,1,16000}]
ex = Table[centeredzdist[[i]][[1]]*centeredzdist[[i]][[2]]*(1/$BINCOUNT), {i, 1, 24000}];
exsquared = Table[centeredzdist[[i]][[1]]^2*centeredzdist[[i]][[2]]*(1/$BINCOUNT), {i, 1,24000}];
var = Total[exsquared] - Total[ex]^2;
stdev = Sqrt[var];
If[curriter != 1,
diff = Table[Abs[centeredzdist[[i]][[2]]-centeredprev[[i]][[2]]],{i,1,24000}];
,
diff = 0;];
sumdif = Total[diff];
output = "---Current iteration: "<>ToString[curriter]<>"\n---------------------------------\nMean: "<>ToString[Total[ex]]<>"\nVariance: "<>ToString[var]<>"\nStandard Deviation: "<>ToString[stdev]<>"\nDifference between this and last iteration: "<>ToString[sumdif]<>"\nzmax: "<>ToString[zmax];
Print[output];
Export["$jobdir/push/$MATRIX-$CONFIGS-"<>ToString[pushamt]<>"/"<>ToString[curriter]<>"/stats.txt",output]; 
EOD

cat > ${statsjobfile} << EOD
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2012
#SBATCH --time=1:00:00
module purge
module load Mathematica
module list
pwd
ls -al $jobdir/$statsfile
$jobdir/$statsfile \$CURRITER \$PUSHAMT


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
			jid=$(sbatch --dependency=afterok:$endproc --export=CURRITER=$CURRITER,MUTABLELOOPSIZE=$(($NOOFBATCH)),MUTABLECONFIGSIZE=$MAXSINGLETHREADSIZE,CURRPROC=$CURRPROC,PUSHAMT=$PUSHAMT ${jobdir}/${rgjobfile} | cut -d ' ' -f4);
			joined="$joined$delim$jid";
			delim=","
			CURRPROC=$(($CURRPROC+1));
		done;
	fi
	
	CURRPROC=$(( ($BATCHSIZE * $NOOFBATCH) + 1 ))
	for  (( j=2;j<=$PROCREMAINDER;j++ ));
	do
		jid=$(sbatch --dependency=afterok:$endproc --export=CURRITER=$CURRITER,MUTABLELOOPSIZE=$((1)),MUTABLECONFIGSIZE=$MAXSINGLETHREADSIZE,CURRPROC=$CURRPROC,PUSHAMT=$PUSHAMT ${jobdir}/${rgjobfile} | cut -d ' ' -f4);
		joined="$joined$delim$jid";
		CURRPROC=$(($CURRPROC + 1));
	done;
	if [[ $CONFIGREMAINDER -ge 1 ]]; 
	then
		jid=$(sbatch --dependency=afterok:$endproc --export=CURRITER=$CURRITER,MUTABLELOOPSIZE=$((1)),MUTABLECONFIGSIZE=$CONFIGREMAINDER,CURRPROC=$CURRPROC,PUSHAMT=$PUSHAMT ${jobdir}/${rgjobfile} | cut -d ' ' -f4);
		joined="$joined$delim$jid";
	fi

	endproc=$(sbatch --dependency=afterok:$joined --export=CURRITER=$CURRITER,PUSHAMT=$PUSHAMT ${jobdir}/${distavgjobfile} | cut -d ' ' -f4);
	if [[ $SYMMETRIZE -eq 1  ]]; then endproc=$(sbatch --dependency=afterok:$endproc --export=CURRITER=$CURRITER ${jobdir}/${symmetrizejobfile} | cut -d ' ' -f4); fi
	endproc=$(sbatch --dependency=afterok:$endproc --export=CURRITER=$CURRITER,PUSHAMT=$PUSHAMT ${jobdir}/${statsjobfile} | cut -d ' ' -f4);
	
done;






