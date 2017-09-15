#!/bin/bash

dirname=$( dirname $0 )
[[ "${dirname:0:1}" != "/" ]] && dirname="$( pwd )/$dirname"
basename=$( basename $0 )

INF=0.00

DATABASE="$dirname/../db/prodb/PRODB"
OUTFILE="$dirname/../db/mpiscaler.I${INF}.`date +'%y%b%d'`.out"
SCOMETH="pro"
MASTERNODENAME=$( hostname )
NODEFILE="$dirname/../var/node.lst"
CPUS=4
MPISCALER="$dirname/mpiscaler"

usage="
This is a helper script to launch \`mpiscaler'.
2008(C)Mindaugas Margelevicius,VU IBT,Vilnius

$basename <Options>

Options:

-d <database>  absolute pathname to database to be read for scaling
       default=$DATABASE
-o <output>    absolute pathname to file to write output to
       default=$OUTFILE
-s (pro|LSO|HDP) Scoring method to use
       default=$SCOMETH

-H <hostname>  Master node name
       default=$MASTERNODENAME
-N <filename>  pathname to file containing list (1/line) of node names
               NOTE: The list should not include this computer
       default=$NODEFILE
-C <number>    number of CPUs to use on nodes
       default=$CPUS
-P <program>   absolute pathname to mpiscaler program
               (must be enabled to be accessed from other nodes)
       default=$MPISCALER
-h             short description

"


while getopts "d:o:s:H:N:C:P:h" Option
do
    case $Option in
        d ) DATABASE=${OPTARG} ;;
        o ) OUTFILE=${OPTARG} ;;
        s ) SCOMETH=${OPTARG} ;;
        H ) MASTERNODENAME=${OPTARG} ;;
        N ) NODEFILE=${OPTARG} ;;
        C ) CPUS=${OPTARG} ;;
        P ) MPISCALER=${OPTARG} ;;
        h ) echo "$usage"; exit 0 ;;
        * ) echo Error: Unrecognized argument.; exit 1 ;;
    esac
done
shift $(( $OPTIND - 1 ))

ERRFILE="$( dirname $OUTFILE )/$( basename $OUTFILE .out ).err"

if [[ -z "$( ls -1 ${DATABASE}* )" ]]; then echo Error: Database does not exist.; exit 1; fi
if [[ -z $NODEFILE || ! -f $NODEFILE ]]; then echo Error: File with the list of nodes does not exist.; exit 1; fi
if [[ -z $MPISCALER || ! -f $MPISCALER ]]; then echo Error: Program $MPISCALER does not exist.; exit 1; fi
if [[ -z $CPUS || $CPUS -le 0 ]]; then echo Error: Wrong number of CPUs specified.; exit 1; fi

NODES=$(< $NODEFILE )


mastercpus=1

##mpiscaler_cmd="$MPISCALER -d $DATABASE -I $INF -o $OUTFILE -s $SCOMETH"
mpiscaler_cmd="$MPISCALER -d $DATABASE -o $OUTFILE -s $SCOMETH"

deployment="-n $mastercpus -host $MASTERNODENAME $mpiscaler_cmd"



for node in $NODES; do
    [[ "${node:0:1}" == "#" ]] && continue
    nodecpus=$CPUS
    nodecpus=$( echo $node | awk -F ':' '{print $2}')
    nodename=$( echo $node | awk -F ':' '{print $1}')
    [[ -z "$nodecpus" ]] && nodecpus=$CPUS
    deployment="$deployment : -n $nodecpus -host $nodename $mpiscaler_cmd"
done

cmd="mpiexec -l $deployment >$ERRFILE 2>&1 &"

## to run by PBS engine
## cmd="mpiexec -l $deployment"

echo $cmd
eval $cmd

