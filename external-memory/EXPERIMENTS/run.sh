#!/usr/bin/env bash

# Bash script for running ANN algorithms

# Fail on error
set -e

# Echo on
# set -x

function usage {
    printf '\nUsage: %s [-ahr]\n\n' $(basename $0) >&2

    printf 'This install script attempts to run ANN experiments\n' >&2

    printf 'options:\n' >&2
    printf -- ' -a: (default) run (A)ll algorithms - good luck!\n' >&2
    printf -- ' -c: clean all!\n' >&2
    printf -- ' -h: print this (H)elp message\n' >&2
    printf -- ' -r <alg>: only run <alg>, where <alg>=LinearScan|iDEC\n' >&2
    exit 2
}

DEBUG=0
# read -p "Run this bash script ('y' for running|'n' for just debugging)? (y/n) " -n 1 -r
# echo    # (optional) move to a new line
# if [[ ! $REPLY =~ ^[Yy]$ ]]
# then
#     DEBUG=1
# fi

# Global variables
WORKING_DIR=$(pwd)
DATASETS=( $( cat dataset_info.txt ) )
# printf "%s\n" "${DATASETS[@]}"

# common parameters
MAX_K=10
APP_RATIO=2
BLK_SIZE=4096


# iDEC
iDEC_M=6
iDEC_ES=0 # 0: Disable early stop 1: enable
# Original parameters
# iDEC_T=(0.000020 0.000029 0.000041 0.000060 0.000086 0.000123 0.000177 0.000255 0.000367 0.000527 0.000759 0.001091 0.001570 0.002258 0.003248 0.004671 0.006720 0.009666 0.013904 0.020000)
# Special parameters in case iDEC_T is not large enough
iDEC_T=(0.00004000 0.00006000 0.00008000 0.00010000 0.00015000 0.00020000 0.00030000)
# iDEC_T=(0.00000001 0.00000003 0.00000009 0.00000029 0.00000087 0.00000267 0.00000818 0.00002500)

# QALSH
QALSH_C=(1.5)

# C2LSH
C_t=0 # C_t: 0 Thereotical Gurantee version 1: Threshold version

# An error exit function
function error_exit
{
    echo "$1" 1>&2
    exit 1
}

# Compiling LinearScan
function CompileLinearScan {
    cd ${WORKING_DIR}
    echo "Compile LinearScan ..."
    cd ../LinearScan

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

    # if [ -f ${dsname}/gnd.txt ]; then
    #     error_exit "Result file ${dsname}/gnd.txt already exists (indicating that you have already run linear scan on dataset ${dsname}) ..."
    # fi

    mkdir -p ${dsname}

    # echo "./linear-scan -k ${MAX_K} -ds ./datasets/${dsh5} -rf ${dsname}/gnd.txt -n ${n} -qn ${qn} -d ${dim} -pf ${dsname}/linear-scan.txt"
    # if [ $DEBUG == 0 ]
    # then
    #     ./linear-scan -k ${MAX_K} -ds ./datasets/${dsh5} -rf ${dsname}/gnd.txt -n ${n} -qn ${qn} -d ${dim} -pf ${dsname}/linear-scan.txt
    # fi

    echo "./linear-scan -k ${MAX_K} -ds ./datasets/${dsh5} -rf ${dsname}/gnd-e.txt -n ${n} -qn ${qn} -d ${dim} -b ${BLK_SIZE} -pf ${dsname}/linear-scan-e.txt"
    if [ $DEBUG == 0 ]
    then
        ./linear-scan -k ${MAX_K} -ds ./datasets/${dsh5} -rf ${dsname}/gnd-e.txt -n ${n} -qn ${qn} -d ${dim} -b ${BLK_SIZE} -pf ${dsname}/linear-scan-e.txt
    fi
    echo "Done\n"
}

# Clean LinearScan
function CleanLinearScan {
    echo "Clean LinearScan ..."
    cd ${WORKING_DIR}
    cd ../LinearScan
    make clean
    echo "Done"
}


# Compile iDEC
function CompileiDEC {
    echo "Compile iDEC ..."
    cd ${WORKING_DIR}
    cd ../iDEC
    make clean
    make
    chmod +x ./idec-compact
    cp ./idec-compact ${WORKING_DIR}
    echo "Done"
}

# Run iDEC
function RuniDEC {
    cd ${WORKING_DIR}
    if [ ! -f ./idec-compact ]; then
        error_exit "Executable file ${WORKING_DIR}/idec-compact does not exist, please compile it and copy it to ${WORKING_DIR} first ..."
    fi

    dsname=$1
    dstrainb=$2
    dstestt=$3
    n=$4
    qn=$5
    dim=$6

    printf "Run iDEC on dataset %s ...\n" "${dsname}"


    mkdir -p iDEC/${dsname}/index
    mkdir -p iDEC/${dsname}/results
    cp ./idec-compact iDEC/${dsname}
    cd iDEC/${dsname}

    m=${iDEC_M}

    printf "Build index for iDEC on dataset %s ...\n" "${dsname}"
    echo "./idec-compact -b ${BLK_SIZE} -d ${dim} -ds ../../${dsname}/${dstrainb} -n ${n} -m ${m} -f ./index"
    if [ $DEBUG == 0 ]
    then
        ./idec-compact -b ${BLK_SIZE} -d ${dim} -ds ../../${dsname}/${dstrainb} -n ${n} -m ${m} -f ./index
    fi

    printf "Build binary files iDEC on dataset %s ...\n" "${dsname}"
    echo "./idec-compact -b ${BLK_SIZE} -d ${dim} -ds ../../${dsname}/${dstrainb} -n ${n} -f ./index"
    if [ $DEBUG == 0 ]
    then
        ./idec-compact -b ${BLK_SIZE} -d ${dim} -ds ../../${dsname}/${dstrainb} -n ${n} -f ./index
    fi

    printf "Query for iDEC on dataset %s ...\n" "${dsname}"
    for t in "${iDEC_T[@]}"
    do

        echo "./idec-compact -f ./index -k ${MAX_K} -o ./results/idec-${m}-${t}.txt -r ../../${dsname}/gnd-e.txt  -q ../../${dsname}/${dstestt} -c ${APP_RATIO} -es ${iDEC_ES} -t ${t}"
        if [ $DEBUG == 0 ]
        then
            ./idec-compact -f ./index -k ${MAX_K} -o ./results/idec-${m}-${t}.txt -r ../../${dsname}/gnd-e.txt  -q ../../${dsname}/${dstestt} -c ${APP_RATIO} -es ${iDEC_ES} -t ${t}
        fi
    done
    echo "Done"
}

# Clean iDEC
function CleaniDEC {
    echo "Clean iDEC ..."
    cd ${WORKING_DIR}
    cd ../iDEC
    make clean
    echo "Done"
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
        *)
            echo "Unknown option $1."
            usage
            exit 1
    esac
}

function single
{
    compile_single $1

    for ds in "${DATASETS[@]}"
    do
        # read dataset information
        IFS='%'
        read -ra ds_info <<< "$ds" # str is read into an array as tokens separated by IFS
        dsname=${ds_info[0]}
        dstrainb=${ds_info[1]}
        dstestb=${ds_info[2]}
        dstraint=${ds_info[3]}
        dstestt=${ds_info[4]}
        n=${ds_info[5]}
        dim=${ds_info[6]}
        qn=${ds_info[7]}
        dsh5=${ds_info[8]}

        cd ${WORKING_DIR}

        if [ ! -f ./$dsname/gnd.txt ]; then
            RunLinearScan $dsname $dsh5 $n $qn $dim
        fi

        case "$1" in
            LinearScan)
                RunLinearScan $dsname $dsh5 $n $qn $dim
                CleanLinearScan
                ;;
            iDEC)
                RuniDEC $dsname $dstrainb $dstestt $n $qn $dim
                CleaniDEC
                ;;
            *)
                echo "Unknown option $1."
                usage
        esac
    done

}


# if [ $# -eq 0 ]
# then
#     all
# else
#     while getopts 'achr:' OPTION
#     do
#         case $OPTION in
#             a)    all ;;
#             c)    clean_all ;;
#             h)    usage ;;
#             r)    single $OPTARG ;;
#             ?)    usage ;;
#         esac
#     done
#     shift $(($OPTIND - 1))
# fi
