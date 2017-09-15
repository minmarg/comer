#!/bin/bash

dirname=$( dirname $0 )
[[ "${dirname:0:1}" != "/" ]] && dirname="$( pwd )/$dirname"
basename=$( basename $0 )


DATABASE=""
INPUTOUT=""
OUTFILE="hdpclust.`date +'%y%b%d'`"
PARPROC=3
SCANS=1000
MHUPS=1
IRGSSCANS=-1
TAU=-1
GAMMA=-1
KAPPA0=-1
NU0=-1
DEGFA=0
SCALE=1
MASTERNODENAME=$( hostname )
NODEFILE="$dirname/../var/node.lst"
CPUS=4
HDPCLUST="$dirname/hdpclust"

usage="
A helper script to launch \`hdpclust'.
2011(C)Mindaugas Margelevicius,VU IBT,Vilnius

$basename <Options>

Options:

-r <input>     Name pattern of restart-with files;
               <input>.grp -- grouped frequencies,
               <input>.par -- file of parameters.
-o <output>    Name pattern of output files.
       default=$OUTFILE
-e             A cluster for each initial sample vector.
-a [1|2|3]     Clustering of parallel processing in turn of
               1 -- single vector/table,
               2 -- single group,
               3 -- all groups.
       default=$PARPROC
-t             Do not apply whole table sampling in a scan.
-k <number>    Number of sampling scans (iterations).
       default=$SCANS
-u <number>    Number of M-H updates in a single scan.
       default=$MHUPS
-m <number>    Number of M-H interm. restricted sampling scans;
               -1 means likelihood sampling instead.
       default=$IRGSSCANS
-l [1|2]       1/2 -- Do split/merge proposals only in M-H updates.
-j             Use modified JN MCMC algorithm in M-H updates.
-d             Select clusters uniformly in M-H updates.
-f             Select cluster vectors uniformly in M-H updates.
-y             Prior value of concentration parameter tau.
               -1 means to sample parameter tau.
       default=$TAU
-g             Prior value of concentration parameter gamma.
               -1 means to sample parameter gamma.
       default=$GAMMA
-b             Prior value of scale factor kappa_0.
               -1 means to sample scale factor.
       default=$KAPPA0
-c             Prior value of degrees of freedom.
               -1 means to sample degrees of freedom.
       default=$NU0
-n <number>    Adjustment to the degrees of freedom in terms of
               dimensions over 2. (Valid values from -1).
       default=$DEGFA
-s <scale>     Scale factor for uninformative prior parameters.
       default=$SCALE
-p             Calculate prior parameters from file of grouped
               frequencies.

-H <hostname>  Master node name
       default=$MASTERNODENAME
-N <filename>  pathname to file containing list (1/line) of node names
               NOTE: The list should not include this computer
       default=$NODEFILE
-C <number>    number of CPUs to use on nodes
       default=$CPUS
-P <program>   absolute pathname to \`hdpclust' program
               (should be accessible from other nodes)
       default=$HDPCLUST
-h             short description

"


while getopts "r:o:ea:tk:u:m:l:jdfy:g:b:c:n:s:pH:N:C:P:h" Option
do
    case $Option in
        r ) INPUTOUT=${OPTARG} ;;
        o ) OUTFILE=${OPTARG} ;;
        e ) CL4EACH="-e" ;;
        a ) PARPROC=${OPTARG} ;;
        t ) NOTSMPL="-t" ;;
        k ) SCANS=${OPTARG} ;;
        u ) MHUPS=${OPTARG} ;;
        m ) IRGSSCANS=${OPTARG} ;;
        l ) SMPROPS="-l${OPTARG}" ;;
        j ) MJNMCMC="-j" ;;
        d ) UNFDISH="-d" ;;
        f ) UNFVECT="-f" ;;
        y ) TAU=${OPTARG} ;;
        g ) GAMMA=${OPTARG} ;;
        b ) KAPPA0=${OPTARG} ;;
        c ) NU0=${OPTARG} ;;
        n ) DEGFA=${OPTARG} ;;
        s ) SCALE=${OPTARG} ;;
        p ) GRPRIOR="-p" ;;
        H ) MASTERNODENAME=${OPTARG} ;;
        N ) NODEFILE=${OPTARG} ;;
        C ) CPUS=${OPTARG} ;;
        P ) HDPCLUST=${OPTARG} ;;
        h ) echo "$usage"; exit 0 ;;
        * ) echo Error: Unrecognized argument.; exit 1 ;;
    esac
done
shift $(( $OPTIND - 1 ))

ERRFILE="$( dirname $OUTFILE )/$( basename $OUTFILE .out ).err"

if [[ ! -f $NODEFILE ]]; then echo Error: File with the list of nodes does not exist.; exit 1; fi
if [[ ! -f $HDPCLUST ]]; then echo Error: Program $HDPCLUST does not exist.; exit 1; fi
if [[ -z $CPUS || $CPUS -le 0 ]]; then echo Error: Wrong number of CPUs specified.; exit 1; fi

[[ -n "$INPUTOUT" && "${INPUTOUT:0:1}" != "/" ]] && INPUTOUT="$( pwd )/$INPUTOUT"
[[ -n "$OUTFILE" && "${OUTFILE:0:1}" != "/" ]] && OUTFILE="$( pwd )/$OUTFILE"

NODES=$(< $NODEFILE )

mastercpus=1

hdpclust_cmd="$HDPCLUST -v -r $INPUTOUT -o $OUTFILE \
$CL4EACH -a $PARPROC $NOTSMPL -k $SCANS -u $MHUPS -m $IRGSSCANS $SMPROPS \
$MJNMCMC $UNFDISH $UNFVECT -y $TAU -g $GAMMA -b $KAPPA0 -c $NU0 \
-n $DEGFA -s $SCALE $GRPRIOR"

oldIFS=$IFS
IFS=$'\n'

for node in $NODES; do
    [[ "${node:0:1}" == "#" ]] && continue
    nodecpus=$CPUS
    nodecpus=$( echo $node | awk -F ':' '{if($2) print $2; else print "1"}')
    nodename=$( echo $node | awk -F ':' '{print $1}')
    nops=$(( nops + nodecpus ))
done

IFS=$oldIFS

cmd="mpiexec -f $NODEFILE -n $nops $hdpclust_cmd >$ERRFILE 2>&1 &"

## to run by PBS engine
## cmd="mpiexec -l $deployment"

echo $cmd
eval $cmd

