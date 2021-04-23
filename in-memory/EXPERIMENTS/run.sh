-#!/usr/bin/env bash

# Bash script for running ANN algorithms (beautified by `beautysh`)

# Fail on error
set -e

# Echo on
# set -x

function usage {
    printf '\nUsage: %s [-ahr]\n\n' $(basename $0) >&2

    printf 'This script attempts to run in-memory ANN experiments\n' >&2

    printf 'options:\n' >&2
    printf -- ' -a: (default) run (A)ll algorithms - good luck!\n' >&2
    printf -- ' -c: clean all!\n' >&2
    printf -- ' -h: print this (H)elp message\n' >&2
    printf -- ' -r <alg>: only run <alg>, where <alg>=LinearScan|iDEC|iDEC-STEAL\n' >&2
    exit 2
}

# set DEBUG to 1 if you only want to debug this script
DEBUG=0

# Global variables
WORKING_DIR=$(pwd)
LINEARSCAN_DIR=../../external-memory/LinearScan
DATASETS=( $( cat dataset_info.txt ) )


# Parameters for algorithms
MAX_K=100
APP_RATIO=2

# A special iDEC T for high dimensional dataset

# # iDEC parameters
iDEC_M=6
iDEC_T=(0.0000000006 0.0000000012 0.000000025 0.000000050 0.000000060 0.000000072 0.000000086 0.000000104 0.000000125 0.000000150 0.000000180 0.000000216 0.000000259 0.000000311 0.000000373 0.000000448 0.000000538 0.000000645 0.000000775 0.000000930 0.000001116 0.000001340 0.000001609 0.000001931 0.000002319 0.000002783 0.000003341 0.000004011 0.000004815 0.000005780 0.000006939 0.000008330 0.000010000 0.000020000)
# iDEC_T=(0.0000000006)
# iDEC_T=(0.000040 0.000057 0.000082 0.000097 0.000116 0.000139 0.000166 0.000199 0.000237 0.000284 0.000339 0.000405 0.000484 0.000579 0.000691 0.000826 0.000987 0.001180 0.001410 0.001684 0.002013 0.002405 0.002874 0.003435 0.004104 0.004904 0.005861 0.007003 0.008368 0.010000)
# A more density t for plotting figure
# iDEC_T=(0.000010 0.000011 0.000013 0.000014 0.000016 0.000018 0.000020 0.000023 0.000026 0.000029 0.000032 0.000036 0.000041 0.000046 0.000052 0.000058 0.000065 0.000073 0.000082 0.000092 0.000104 0.000117 0.000131 0.000148 0.000166 0.000187 0.000210 0.000236 0.000265 0.000298 0.000335 0.000377 0.000424 0.000476 0.000536 0.000602 0.000677 0.000761 0.000855 0.000962 0.001081 0.001215 0.001366 0.001536 0.001727 0.001941 0.002183 0.002454 0.002759 0.003101 0.003486 0.003919 0.004406 0.004954 0.005569 0.006261 0.007038 0.007912 0.008895 0.010000)
# Some extra iDEC points
# iDEC_T=(0.011 0.012 0.013 0.014)


# A special iDEC T for high dimensional dataset
# iDEC_T=(0.025 0.03 0.035 0.040 0.045)

# An error exit function
function error_exit
{
    echo "$1" 1>&2
    exit 1
}

# Compiling LinearScan
function CompileLinearScan {
    cd ${WORKING_DIR}
    cd ${LINEARSCAN_DIR}
    echo "Compile LinearScan ..."

    make clean
    make linear-scan
    cp linear-scan ${WORKING_DIR}

    cd ${WORKING_DIR}
    chmod +x ./linear-scan
    echo "Done"
}


# Run LinearScan
function RunLinearScan {
    cd ${WORKING_DIR}
    if [ ! -f ./linear-scan ]; then
        error_exit "Executable file ${WORKING_DIR}/linear-scan does not exist, please compile it and copy it to ${WORKING_DIR} first ..."
    fi

    dsname=$1
    dsh5=$2
    n=$3
    qn=$4
    dim=$5

    printf "Run linear scan on dataset %s ...\n" "${dsname}"

    if [ -f ${dsname}/gnd.txt ]; then
        error_exit "Result file ${dsname}/gnd.txt already exists (indicating that you have already run linear scan on dataset ${dsname}) ..."
    fi

    mkdir -p ${dsname}

    echo "./linear-scan -k ${MAX_K} -ds ./datasets/${dsh5} -rf ${dsname}/gnd.txt -n ${n} -qn ${qn} -d ${dim} -pf ${dsname}/linear-scan.txt"
    if [ $DEBUG == 0 ]
    then
        ./linear-scan -k ${MAX_K} -ds ./datasets/${dsh5} -rf ${dsname}/gnd.txt -n ${n} -qn ${qn} -d ${dim} -pf ${dsname}/linear-scan.txt
    fi

    echo "Done\n"
}

# Clean LinearScan
function CleanLinearScan {
    echo "Clean LinearScan ..."
    cd ${WORKING_DIR}
    cd ${LINEARSCAN_DIR}
    make clean
    echo "Done"
}

# Compile iDEC
function CompileiDEC {
    cd ${WORKING_DIR}
    echo "Compile iDEC ..."
    cd ../iDEC
    make clean
    make
    chmod +x ./idec
    cp ./idec ${WORKING_DIR}
    cd ${WORKING_DIR}
    echo "Done."
}

# Run iDEC
function RuniDEC {

    cd ${WORKING_DIR}
    if [ ! -f ./idec ]; then
        error_exit "Executable file ${WORKING_DIR}/idec does not exist, please compile it and copy it to ${WORKING_DIR} first ..."
    fi

    dsname=$1
    dstrainb=$2
    dstestb=$3
    n=$4
    qn=$5
    dim=$6

    printf "Run iDEC on dataset %s ...\n" "${dsname}"

    mkdir -p iDEC/${dsname}/index
    mkdir -p iDEC/${dsname}/results
    cp ./idec iDEC/${dsname}
    cd iDEC/${dsname}

    echo "./idec -d $dim -n $n -ds ../../${dsname}/${dstrainb} -if index -m ${iDEC_M}"
    if [ $DEBUG == 0 ]
    then
        ./idec -d $dim -n $n -ds ../../${dsname}/${dstrainb} -if index -m ${iDEC_M}
    fi

    for t in "${iDEC_T[@]}"
    do
        echo "./idec -d $dim -n $n -ds ../../${dsname}/${dstrainb} -qs ../../${dsname}/${dstestb} -if index -rf ./results/idec-$t.txt -gt ../../${dsname}/gnd.txt -qn $qn -t $t"
        if [ $DEBUG == 0 ]
        then
            ./idec -d $dim -n $n -ds ../../${dsname}/${dstrainb} -qs ../../${dsname}/${dstestb} -if index -rf ./results/idec-$t.txt -gt ../../${dsname}/gnd.txt -qn $qn -t $t
        fi
    done
}

# Clean iDEC
function CleaniDEC {
    echo "Clean iDEC ..."
    cd ${WORKING_DIR}
    cd ../iDEC
    make clean
    echo "Done"
}

# Compile iDEC
function CompileiDECSteal {
    cd ${WORKING_DIR}
    echo "Compile iDEC-STEAL ..."
    cd ../iDEC
    make clean
    make idec-steal
    chmod +x ./idec-steal
    cp ./idec-steal ${WORKING_DIR}
    cd ${WORKING_DIR}
    echo "Done."
}

# Run iDEC
function RuniDECSteal {

    cd ${WORKING_DIR}
    if [ ! -f ./idec-steal ]; then
        error_exit "Executable file ${WORKING_DIR}/idec-steal does not exist, please compile it and copy it to ${WORKING_DIR} first ..."
    fi

    dsname=$1
    dstrainb=$2
    dstestb=$3
    n=$4
    qn=$5
    dim=$6

    printf "Run iDECSteal on dataset %s ...\n" "${dsname}"

    mkdir -p iDECSteal/${dsname}/index
    mkdir -p iDECSteal/${dsname}/results
    cp ./idec-steal iDECSteal/${dsname}
    cd iDECSteal/${dsname}

    echo "./idec-steal -d $dim -n $n -ds ../../${dsname}/${dstrainb} -if index -m ${iDEC_M}"
    if [ $DEBUG == 0 ]
    then
        ./idec-steal -d $dim -n $n -ds ../../${dsname}/${dstrainb} -if index -m ${iDEC_M}
    fi

    for t in "${iDEC_T[@]}"
    do
        echo "./idec-steal -d $dim -n $n -ds ../../${dsname}/${dstrainb} -qs ../../${dsname}/${dstestb} -if index -rf ./results/idec-$t.txt -gt ../../${dsname}/gnd.txt -qn $qn -t $t"
        if [ $DEBUG == 0 ]
        then
            ./idec-steal -d $dim -n $n -ds ../../${dsname}/${dstrainb} -qs ../../${dsname}/${dstestb} -if index -rf ./results/idec-$t.txt -gt ../../${dsname}/gnd.txt -qn $qn -t $t
        fi
    done
}

# Clean iDEC
function CleaniDECSteal {
    echo "Clean iDEC-STEAL ..."
    cd ${WORKING_DIR}
    cd ../iDEC
    make clean
    echo "Done"
}

# Clean all
function clean_all
{
    CleanLinearScan
    CleaniDEC
    CleaniDECSteal
}

# Run all algorithms
function all
{

    # CompileLinearScan
    CompileiDEC
    CompileiDECSteal

    for ds in "${DATASETS[@]}"
    do
        # read dataset information
        IFS='%'
        read -ra ds_info <<< "$ds" # str is read into an array as tokens separated by IFS

        # TODO: change accordingly
        dsname=${ds_info[0]}
        dstrainfvecs=${ds_info[1]}
        dstestfvecs=${ds_info[2]}
        n=${ds_info[3]}
        dim=${ds_info[4]}
        qn=${ds_info[5]}
        dsh5=${ds_info[6]}
        dstrainflat=${ds_info[7]}
        dstestflat=${ds_info[8]}

        # RunLinearScan $dsname $dsh5 $n $qn $dim

        RuniDEC $dsname $dstrainflat $dstestflat $n $qn $dim

        RuniDECSteal $dsname $dstrainflat $dstestflat $n $qn $dim

    done

    clean_all
}

# Compile a single algorithm
function compile_single
{
    case "$1" in
        LinearScan)
            CompileLinearScan
            ;;
        iDEC)
            CompileiDEC
            ;;
        iDEC-STEAL)
            CompileiDECSteal
        *)
            echo "Unknown option $1."
            usage
            exit 1
    esac
}

function single
{
    # CompileLinearScan
    compile_single $1

    for ds in "${DATASETS[@]}"
    do
        # read dataset information
        IFS='%'
        # TODO: change accordingly
        read -ra ds_info <<< "$ds" # str is read into an array as tokens separated by IFS
        dsname=${ds_info[0]}
        dstrainfvecs=${ds_info[1]}
        dstestfvecs=${ds_info[2]}
        n=${ds_info[3]}
        dim=${ds_info[4]}
        qn=${ds_info[5]}
        dsh5=${ds_info[6]}
        dstrainflat=${ds_info[7]}
        dstestflat=${ds_info[8]}

        cd ${WORKING_DIR}

        if [ ! -f ./$dsname/gnd.txt ]; then
            RunLinearScan $dsname $dsh5 $n $qn $dim
        fi

        case "$1" in
            LinearScan)
                RunLinearScan $dsname $dsh5 $n $qn $dim
                ;;
            iDEC)
                RuniDEC $dsname $dstrainflat $dstestflat $n $qn $dim
                ;;
            iDEC-STEAL)
                RuniDECSteal $dsname $dstrainflat $dstestflat $n $qn $dim
                ;;
            *)
                echo "Unknown option $1."
                usage
        esac
    done
    clean_all
}

# ---
if [ $# -eq 0 ]
then
    all
else
    while getopts 'achr:' OPTION
    do
        case $OPTION in
            a)    all ;;
            c)    clean_all ;;
            h)    usage ;;
            r)    single $OPTARG ;;
            ?)    usage ;;
        esac
    done
    shift $(($OPTIND - 1))
fi