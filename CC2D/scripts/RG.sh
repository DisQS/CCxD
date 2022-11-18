#!/bin/bash




# Parameter initialization with 'default' values

SYMMETRISE=0
MATRIX="C"
LAUNDER=0
CONFIGS=1000000
ITERATIONS=30
PRINTOUTFREQ=30
TASKFARM=0
READIN=0
# Flag and argument capture loop

while getopts 'slc:i:p:m:hto:r:' OPTION; do
	case "$OPTION" in
	r)
		READIN=1
		FILENAME=$OPTARG
		echo "Reading in data from $OPTARG"
		;;
	o)
		echo "Offset feature not available yet"
		;;
	s)
		SYMMETRISE=1
		echo "Output will be symmetrised"
		;;
	h)
		echo "*****************************************************"
		echo "*            Renormalisation Script Help            *"
		echo "*                                                   *"
		echo "*    -m -> Scattering matrix in use('C' or 'S')     *"
		echo "*    -c -> Number of configs (default 1M)           *"
		echo "*    -i -> Number of RG Steps (default 30)          *"
		echo "*    -p -> Print out frequency (default 30)         *"
		echo "*    -s -> Symmetrise (no argument)                 *"
		echo "*    -l -> Launder (no argument)                    *"
		echo "*    -t -> Run on the taskfarm (default is local)   *"
		echo "*    -r -> Read in file (not available yet)         *"
		echo "*                                                   *"
		echo "*****************************************************"
		exit 1
		;;
	l)
		LAUNDER=1
		echo "Output will be laundered"
		;;
	c)
		CONFIGS=$OPTARG
		echo "Running for $OPTARG configurations"
		;;
	m)
		echo $OPTARG
		if [[ $OPTARG = "C" ]]
		then
			echo "Calculating RG steps using the Cain S-Matrix"
		elif [[ $OPTARG = "S" ]]
		then
			MATRIX=$OPTARG
			echo "Calculating RG steps using the Son S-Matrix"
		else
			echo "Invalid argument for -m, consult -h for help"
			exit 1
		fi
		;;
	i)
		ITERATIONS=$OPTARG
		echo "Performing $OPTARG RG steps"
		;;
	p)
		PRINTOUTFREQ=$OPTARG
		echo "Program will export data $OPTARG times"
		;;
	t)
		TASKFARM=1
		echo "Program will be run on the taskfarm"
		;;
	?)
		echo "Invalid input, consult -h for help"
		exit 1
		;;
		
	esac
done
shift "$(($OPTIND -1))"



# Settings for directories

currdir=`pwd`'/../data'
cd $currdir
jobdir="$MATRIX-RG-$CONFIGS"
mkdir -p $jobdir

jobname=$jobdir-$ITERATIONS
echo $jobname

jobfile=`printf "$jobname.sh"`
wlsfile=`printf "$jobname.wls"`
logfile=`printf "$jobname.log"`

# Settings for parallel submission

cd $jobdir

# Generate folders for raw data and distribution data to be exported to
 
mkdir -p raw
mkdir -p dists


cat > ${wlsfile} << EOD
#!/usr/bin/env wolframscript 
Print["(*Preliminaries*)"];
outputfreq = $PRINTOUTFREQ
size=$CONFIGS;
maxRGSteps=$ITERATIONS;
matrix=$MATRIX;
launder=$LAUNDER;
symmetrise=$SYMMETRISE;
Print["Symmetrise it : ",$SYMMETRISE];

Launder[sop_, min_, max_] := (
	binWidth = (max - min)/(Length[sop]);
	normed = sop/(Total[sop])*  binWidth;
	hmax = Max[normed];
	newData = {};
	While[Length[newData] < size,
		prospectpoint = {RandomReal[{min, max}], 
		RandomReal[{0, hmax}]};
		binNo = Ceiling[(prospectpoint[[1]] - min)/binWidth];
		If[normed[[binNo]] >= prospectpoint[[2]], 
			AppendTo[newData, prospectpoint[[1]]]
		];
	];
	Return[newData];
); 
If[$MATRIX === "S",

tp[theta1_, theta2_, theta3_, theta4_, theta5_, phi1_, phi2_, phi3_, phi4_] := 
	Abs[(-Cos[theta1] (Exp[I phi4] Cos[theta3] Cos[theta4] + I Cos[theta5] (Exp[I phi3] Sin[theta2] Sin[theta3] Sin[theta4] - 1)) - Exp[I phi1] Cos[theta2] Cos[theta3] Cos[theta5] - I Exp[I (phi1 + phi4)] Cos[theta4] Cos[theta2] + I Exp[I (phi1 + phi4 - phi2)] Cos[theta2] Cos[theta4] Sin[theta1] Sin[theta3] Sin[theta5])/
	(Exp[I phi4] Cos[theta3] Cos[theta4] Cos[theta5] + Exp[I phi1] Cos[theta1] Cos[theta2] (Cos[theta3] +  I Exp[I phi4] Cos[theta4] Cos[theta5]) + I (-1 + Exp[I phi3] Sin[theta2] Sin[theta3] Sin[theta4] + Exp[I phi2] Sin[theta1] (Sin[theta3] - Exp[I phi3] Sin[theta2] Sin[theta4]) Sin[theta5]))]
,

tp[t1_, r1_, t2_, r2_, t3_, r3_, t4_, r4_, t5_, r5_, phi1_, phi2_, phi3_, phi4_] := 
	Abs[(-Exp[I (phi1 + phi4 - phi2)] r1 r3 r5 t2 t4 + Exp[I (phi1 + phi4)] t2 t4 - Exp[I phi4] t1 t3 t4 + t1 t5 + Exp[I phi3] r2 r3 r4 t1 t5 - Exp[I phi1] t2 t3 t5)/
	(-1 - Exp[I phi3] r2 r3 r4 + Exp[I phi2] r1 r3 r5 + Exp[I (phi2 + phi3)] r1 r2 r4 r5 + Exp[I phi1] t1 t2 t3 - Exp[I (phi1 + phi4)] t1 t2 t4 t5 + Exp[I phi4] t3 t4 t5)];
];
If[$SYMMETRISE===1,If[$READIN===0,adjustedsize = 2 * size],adjustedsize=size];
If[$READIN===1,
	last = Flatten[Import["$FILENAME","Data"]]
,
	If[$MATRIX === "S",
		Print["Creating initial distribution for theta"];
		last = ParallelTable[RandomReal[{0,Pi/2}],{i,1,adjustedsize}];
	,
		Print["Creating initial distribution for t"];
		last = ParallelTable[Sqrt[RandomReal[]],{i,1,adjustedsize}]
	];
	Print["Finished distribution creation, moving onto evaluating renormalized coefficients for ", $ITERATIONS, " steps"]
];
nprint = Floor[$ITERATIONS/outputfreq];
Do[
	Print["--- working on rg step: ", ind];
	If[$MATRIX==="S",
		tpdata = Table[With[{
			t1 = last[[RandomInteger[{1, adjustedsize}]]],
			t2 = last[[RandomInteger[{1, adjustedsize}]]],
			t3 = last[[RandomInteger[{1, adjustedsize}]]],
			t4 = last[[RandomInteger[{1, adjustedsize}]]],
			t5 = last[[RandomInteger[{1, adjustedsize}]]]},
			tp[t1, t2, t3,t4, t5, RandomReal[{0, 2 Pi}], RandomReal[{0, 2 Pi}], RandomReal[{0, 2 Pi}], RandomReal[{0, 2 Pi}]]],
		{i, 1, size}];
	,
	
		tpdata = Table[With[{
			t1 = last[[RandomInteger[{1, adjustedsize}]]],
			t2 = last[[RandomInteger[{1, adjustedsize}]]],
			t3 = last[[RandomInteger[{1, adjustedsize}]]],
			t4 = last[[RandomInteger[{1, adjustedsize}]]],
			t5 = last[[RandomInteger[{1, adjustedsize}]]]},
			tp[t1, Sqrt[1 - t1^2], t2, Sqrt[1 - t2^2], t3, Sqrt[1 - t3^2], t4, Sqrt[1 - t4^2], t5, Sqrt[1 - t5^2], RandomReal[{0, 2 Pi}], RandomReal[{0, 2 Pi}], RandomReal[{0, 2 Pi}], RandomReal[{0, 2 Pi}]]], 
		{i, 1, size}]
	];

	If[$MATRIX==="S",
		last = tpdata; thdata = tpdata; tpdata = Cos[thdata]; thdist = BinCounts[thdata,{0,Pi/2,0.001}];
	,
		last = tpdata;
	];
	
	zdata = Log[(1/tpdata^2)-1];
	gdata = tpdata^2;
	tdist = BinCounts[tpdata,{0,1,0.001}];
	gdist = BinCounts[gdata,{0,1,0.001}];
	qdist = BinCounts[zdata,{-8,8,0.004}];
	If[Mod[ind, nprint] == 0,
		Print["   saving "];
		If[$MATRIX==="S",
			Export["$jobdir/raw/"<>"$MATRIX"<>"-Thraw-"<>"$CONFIGS"<>"-"<>ToString[ind]<>".nc", thdata];
			Export["$jobdir/dists/"<>"$MATRIX"<>"-Thdist-" <> "$CONFIGS" <> "-" <> ToString[ind] <> ".txt", thdist];
		];
		Export["$jobdir/raw/"<>"$MATRIX"<>"-Traw-"<>"$CONFIGS"<>"-"<>ToString[ind]<>".nc", tpdata];
		Export["$jobdir/dists/"<>"$MATRIX"<>"-Tdist-" <> "$CONFIGS" <> "-" <> ToString[ind] <> ".txt", tdist];

		Export["$jobdir/raw/"<>"$MATRIX"<>"-Graw-" <> "$CONFIGS" <> "-" <> ToString[ind] <> ".nc", gdata];
		Export["$jobdir/dists/"<>"$MATRIX"<>"-Gdist-" <> "$CONFIGS" <> "-" <> ToString[ind] <> ".txt", gdist];

		Export["$jobdir/raw/"<>"$MATRIX"<>"-Zraw-" <> "$CONFIGS" <> "-" <> ToString[ind] <> ".nc", zdata];
		Export["$jobdir/dists/"<>"$MATRIX"<>"-Qdist-" <> "$CONFIGS" <> "-" <> ToString[ind] <> ".txt", qdist];
	];
	If[$SYMMETRISE===1,
		zsym = Join[zdata,-zdata];
		If[$MATRIX==="S",
			last = ArcCos[Sqrt[1/(Exp[zsym] + 1)]];
		,
			last = Sqrt[1/(Exp[zsym]+1)]
		]
	];
	If[launder==1,
		If[$MATRIX==="S",
			bins = BinCounts[last,{0,Pi/2,0.001}];
			last = Launder[bins,0,Pi/2];
		,
			bins = BinCounts[last,{0,1,0.001}];
			 last = Launder[bins,0,1]
		];
	];
	If[$MATRIX==="S",

		]

,{ind, 1, $ITERATIONS}];
maindir="$currdir";


Print[maindir];
Print["--- FINISHED!"];
EOD

# Generate job file to be run

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

# Make files executable

chmod 755 ${jobdir}/${jobfile}
chmod 755 ${jobdir}/${wlsfile}
##(sbatch -q devel $jobdir/${jobfile}) # for queueing system

if [[ TASKFARM -eq 1 ]]
then
	sbatch ${jobdir}/${jobfile} # for queueing system
else
	(source ${jobdir}/${jobfile} ) >& ${jobdir}/${logfile} & # for parallel shell execution
fi
#echo "<return>"
sleep 1

cd ..

