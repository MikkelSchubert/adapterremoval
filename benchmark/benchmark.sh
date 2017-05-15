#!/bin/bash
#
# Copyright (c) 2015 Mikkel Schubert <MikkelSch@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

set -o nounset # Fail on unset variables
set -o errexit # Fail on uncaught non-zero returncodes
set -o pipefail # Fail is a command in a chain of pipes fails


###############################################################################
## BENCHMARK PARAMETERS

# Number of replicates per test
NUM_REPLICATES=10

# Read lengths to examine
# Lengths > 100 use an interpolated profile, and should therefore not be used
# to estimate anything other than runtime (see 'simulate_reads').
READ_LENGTHS=(100 200)

# Insert sizes of reads to simulate for adapter ID
ADAPTER_ID_INSERT_SIZES=($(seq 250 5 350))


# Maximum number of threads to use (testing 1 .. max) where supported
MAX_THREADS=4

# Number of read (pairs) to simulate using pIRS for each replicate
SIMULATED_NREADS=1000000

REFSEQ="results/reference.fasta"

SIMULATED_PREFIX="results/simulated/reads"
SIMULATED_MIXED_PREFIX="results/simulated/mixed"
SIMULATED_ADAPTER_ID_PREFIX="results/simulated/adapter_id"


###############################################################################
## PRE-BENCHMARK CHECKS

function check_for_executable()
{
    echo -n "Checking for $1 executable '$2': " > /dev/stderr
    if [ -x "$2" ];
    then
        echo -e "OK" > /dev/stderr
    elif which "$2" &> /dev/null;
    then
        echo -e "'$(which "$2" | head -n1)'" > /dev/stderr
    else
        echo -e "ERROR, NOT FOUND!" > /dev/stderr
        exit 1;
    fi
}


function check_for_jar()
{
    echo -n "Checking for $1 jar at '$2': " > /dev/stderr
    if [ -e "$2" ];
    then
        echo -e "OK" > /dev/stderr
    else
        echo -e "ERROR, NOT FOUND!" > /dev/stderr
        exit 1;
    fi
}


EXEC_ADAPTERREMOVAL1x="bin/AdapterRemoval-1.5.4"
check_for_executable "AdapterRemoval v1.x" "${EXEC_ADAPTERREMOVAL1x}"
# https://github.com/slindgreen/AdapterRemoval/raw/master/AdapterRemoval-1.5.4.tar.gz

EXEC_ADAPTERREMOVAL2x="bin/AdapterRemoval-2.1.3"
check_for_executable "AdapterRemoval v2.x" "${EXEC_ADAPTERREMOVAL2x}"
# https://github.com/MikkelSchubert/adapterremoval

EXEC_LEEHOM="bin/leeHom_patched"
check_for_executable "leeHom" "${EXEC_LEEHOM}"
# https://github.com/grenaud/leeHom
# -- Apply patches/leeHom.patch

EXEC_SKEWER="bin/skewer"
check_for_executable "Skewer" "${EXEC_SKEWER}"
# https://github.com/relipmoc/skewer

EXEC_ALIENTRIMMER="bin/AlienTrimmer"
check_for_executable "AlienTrimmer" "${EXEC_ALIENTRIMMER}"
# ftp://ftp.pasteur.fr/pub/GenSoft/projects/AlienTrimmer/

EXEC_SCYTHE="bin/scythe"
check_for_executable "Scythe" "${EXEC_SCYTHE}"
# https://github.com/vsbuffalo/scythe

EXEC_CUTADAPT="bin/cutadapt"
check_for_executable "Cutadapt" "${EXEC_CUTADAPT}"
# https://code.google.com/p/cutadapt/

EXEC_FLEXBAR="bin/flexbar"
check_for_executable "Flexbar" "${EXEC_FLEXBAR}"
# http://sourceforge.net/projects/flexbar/

EXEC_PIRS="bin/pirs_patched"
check_for_executable "pIRS (with adapters)" "${EXEC_PIRS}"
# ftp://ftp.genomics.org.cn/pub/pIRS/
# -- Apply patches/pirs.patch

JAR_TRIMMOMATIC="bin/trimmomatic-0.33.jar"
check_for_jar "Trimmomatic" ${JAR_TRIMMOMATIC}
# http://www.usadellab.org/cms/?page=trimmomatic

EXEC_PEAT="./bin/PEAT"
check_for_executable "PEAT" ${EXEC_PEAT}
# https://github.com/jhhung/PEAT

EXEC_PEAR="bin/pear_patched"
check_for_executable "PEAR" ${EXEC_PEAR}
# http://sco.h-its.org/exelixis/web/software/pear/
# -- Apply patches/pear.patch

EXEC_MINION="bin/minion"
check_for_executable "minion (kraken)" ${EXEC_MINION}
# http://www.ebi.ac.uk/research/enright/software/kraken

EXEC_FASTQ_MCF="bin/fastq-mcf"
check_for_executable "fastq-mcf" ${EXEC_FASTQ_MCF}
# https://code.google.com/p/ea-utils/


EXEC_TIME=/usr/bin/time
check_for_executable "GNU time" ${EXEC_TIME}
# Needed for time / RAM usage

EXEC_JAVA="bin/java"
check_for_executable "Java JRE" ${EXEC_JAVA}
# Needed for Trimmomatic


# Script for evaluating trimming / collapsing RESULTS
SCRIPT_EVALUATE="./scripts/evaluate.py"
check_for_executable "Evaluation script for read trimming / merging" ${SCRIPT_EVALUATE}

# Script for evaluating trimming / collapsing RESULTS
SCRIPT_EVALUATE_ID="./scripts/evaluate_id.py"
check_for_executable "Evaluation script for adapter identification" ${SCRIPT_EVALUATE_ID}


# Script for merging tables generated by ${SCRIPT_EVALUATE}
SCRIPT_MERGE="./scripts/merge_tables.py"
check_for_executable "Table merging script" ${SCRIPT_MERGE}


###############################################################################

function shuffle_and_run()
{
    echo > /dev/stderr
    echo "Shuffling batch ..." > /dev/stderr

    python -c "import sys, random
lines = sys.stdin.readlines()
random.shuffle(lines)
sys.stdout.write(''.join(lines))" |
    while read command;
    do
        ${command};
    done
}


function do_run_piped()
{
    DST=$1
    DST_TIME=$2
    DST_STDOUT=$3
    DST_STDERR=$4
    shift 4

    if [ -e "${DST}.table" ];
    then
        echo "Skipping ${DST}.table" > /dev/stderr
    else
        echo "Building ${DST}.table" > /dev/stderr
        rm -rf "${DST:?}/"
        mkdir -p "${DST}/"

        if ! ${EXEC_TIME} --verbose --output "${DST_TIME}" "$@" \
            > "${DST_STDOUT}" 2> "${DST_STDERR}";
        then
            echo "Error running command!" > /dev/stderr
            exit 1
        fi
    fi
}


function do_run()
{
    DST=$1
    shift 1

    do_run_piped "${DST}" "${DST}/time" "${DST}/stdout.txt" "${DST}/stderr.txt" "$@"
}


function do_evaluate()
{
    DST=$1

    if [ ! -e "${DST}.table" ];
    then
        ${SCRIPT_EVALUATE} "$@"
    fi
}


###############################################################################


function fetch_reference()
{
    echo "------------------------------------------------------------" > /dev/stderr
    echo "Fetching reference sequence ..."  > /dev/stderr

    if [ -e "results/reference.fasta" ];
    then
        echo "Reference sequence found; skipping ..."  > /dev/stderr
        echo ""  > /dev/stderr
    else
        echo ""  > /dev/stderr
        wget -O results/reference.fasta.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz
        gunzip results/reference.fasta.gz
    fi


}


function simulate_reads()
{
    echo "------------------------------------------------------------" > /dev/stderr
    echo "Simulating reads ..."  > /dev/stderr
    echo ""  > /dev/stderr

    for readlen in ${READ_LENGTHS[*]};
    do
        INSERT_MEAN=$(((${readlen} * 3) / 2))
        INSERT_SD=$((${INSERT_MEAN} / 2))

        if [ "${readlen}" -gt 100 ];
        then
            # Use fake profiles, built using scripts/extend_profile.py
            PROFILE_CLI="-b profiles/phixv2.InDel.matrix -s profiles/humNew.PE100.matrix.gz"
        else
            PROFILE_CLI=
        fi

        for run_n in ${REPLICATES};
        do
            DST="${SIMULATED_PREFIX}_${run_n}_${readlen}"

            if [ -e "${DST}.read.info" ];
            then
                echo "    Skipping ${DST}.*" > /dev/stderr
            else
                echo "    Simulating reads run=${run_n}, l=${readlen}, m=${INSERT_MEAN}, v=${INSERT_SD} ..." > /dev/stderr
                rm -rf "${DST:?}/"
                mkdir -p "${DST}/"

                # -c 0 = uncompressed output
                if ! ${EXEC_PIRS} simulate ${PROFILE_CLI} -x "${SIMULATED_NREADS}" -l "${readlen}" -i "${REFSEQ}" -c 0 -m "${INSERT_MEAN}" -v "${INSERT_SD}" -Q 33 -o "${DST}/reads" \
                    > "${DST}/reads.stdout" 2> "${DST}/reads.stderr";
                then
                    echo "Error simulated reads ..." > /dev/stderr
                    exit 1
                fi

                gunzip "${DST}/reads_${readlen}_${INSERT_MEAN}.read.info.gz"
                ln -sf "$(basename "${DST}")/reads_${readlen}_${INSERT_MEAN}_1.fq" "${DST}_1.fq"
                ln -sf "$(basename "${DST}")/reads_${readlen}_${INSERT_MEAN}_2.fq" "${DST}_2.fq"
                ln -sf "$(basename "${DST}")/reads_${readlen}_${INSERT_MEAN}.read.info" "${DST}.read.info"
            fi
        done
    done

    echo > /dev/stderr
}



function simulate_mixed_reads()
{
# Random adapter pairs generated as follows:
#
# $ scripts/shuffle_fasta.py adapters/adapter_1.fasta
# Seed = 209548449294681565
# Seq  = AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCACCTAATCTCGTATGCCGTCTTCTGCTTG
# New  = AAACTTGCTCTGTGCCCGCTCCGTATGTCACAACAGTGCGTGTATCACCTCAATGCAGGACTCA
# New  = CTAATTTGCCGTAGCGACGTACTTCAGCCTCCAGGAATTGGACCCTTACGCACACGCATTCATG
# New  = GTTCATACGACGACGACCAATGGCACACTTATCCGGTACTTGCGTTTCAATGCGCATGCCCCAT
# New  = CCATGCCCCGAAGATTCCTATACCCTTAAGGTCGCAATTGTTCGAGTAAGCTGTACGCGCCCAT
#
# $ scripts/shuffle_fasta.py adapters/adapter_2.fasta
# Seed = 1852992042931739018
# Seq  = AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
# New  = GATCGGGAGTAATTTGGAGGCAGTAGTTCGTCGAAACTCGGAGCGTCTTTAGCAGGAG
# New  = TACCGTGAAAGGTGCGCTTAGTGGCATATGCGTTAAGAGCTAGGTAACGGTCTGGAGG
# New  = TAAGAAACTCGGAGTTTGGCCTGCGAGGTAGCTTGGGTGTTATGAAGAACGGCATGCG
# New  = GTTGCATTGACCCGAAGGGCTCGATGTTTAGGGAGGTCAGAAGTTGAGCGGGTTCAAA

    echo "------------------------------------------------------------" > /dev/stderr
    echo "Simulating mixed reads ..."  > /dev/stderr
    echo ""  > /dev/stderr

    for readlen in ${READ_LENGTHS[*]};
    do
	    INSERT_MEAN=$(((${readlen} * 3) / 2))
	    INSERT_SD=$((${INSERT_MEAN} / 2))

        if [ "${readlen}" -gt 100 ];
        then
            # Use fake profiles, built using scripts/extend_profile.py
            PROFILE_CLI="-b profiles/phixv2.InDel.matrix -s profiles/humNew.PE100.matrix.gz"
        else
            PROFILE_CLI=
        fi

	    for run_n in ${REPLICATES};
	    do
	        DST="${SIMULATED_MIXED_PREFIX}_${run_n}_${readlen}"

	        if [ -e "${DST}.read.info" ];
	        then
	            echo "    Skipping ${DST}.*" > /dev/stderr
	        else
	            echo "    Simulating mixed reads run=${run_n}, l=${readlen}, m=${INSERT_MEAN}, v=${INSERT_SD} ..." > /dev/stderr
	            rm -rf "${DST:?}/"
	            mkdir -p "${DST}/"

	            # -c 0 = uncompressed output
	            # -x ${NREADS}
	            if ! ${EXEC_PIRS} simulate ${PROFILE_CLI} \
	                -x "${SIMULATED_NREADS}" -l "${readlen}" -i "${REFSEQ}" \
	                -c 0 -m "${INSERT_MEAN}" -v "${INSERT_SD}" -Q 33 \
	                -o "${DST}/reads" \
	                -1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCACCTAATCTCGTATGCCGTCTTCTGCTTG \
	                -1 AAACTTGCTCTGTGCCCGCTCCGTATGTCACAACAGTGCGTGTATCACCTCAATGCAGGACTCA \
	                -1 CTAATTTGCCGTAGCGACGTACTTCAGCCTCCAGGAATTGGACCCTTACGCACACGCATTCATG \
	                -1 GTTCATACGACGACGACCAATGGCACACTTATCCGGTACTTGCGTTTCAATGCGCATGCCCCAT \
	                -1 CCATGCCCCGAAGATTCCTATACCCTTAAGGTCGCAATTGTTCGAGTAAGCTGTACGCGCCCAT \
	                -2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
	                -2 GATCGGGAGTAATTTGGAGGCAGTAGTTCGTCGAAACTCGGAGCGTCTTTAGCAGGAG \
	                -2 TACCGTGAAAGGTGCGCTTAGTGGCATATGCGTTAAGAGCTAGGTAACGGTCTGGAGG \
	                -2 TAAGAAACTCGGAGTTTGGCCTGCGAGGTAGCTTGGGTGTTATGAAGAACGGCATGCG \
	                -2 GTTGCATTGACCCGAAGGGCTCGATGTTTAGGGAGGTCAGAAGTTGAGCGGGTTCAAA \
	                > "${DST}/stdout.txt" 2> "${DST}/stderr.txt";
	            then
	                echo "Error simulated reads ..." > /dev/stderr
	                exit 1
	            fi

	            gunzip "${DST}/reads_${readlen}_${INSERT_MEAN}.read.info.gz"
	            ln -sf "$(basename "${DST}")/reads_${readlen}_${INSERT_MEAN}_1.fq" "${DST}_1.fq"
	            ln -sf "$(basename "${DST}")/reads_${readlen}_${INSERT_MEAN}_2.fq" "${DST}_2.fq"
	            ln -sf "$(basename "${DST}")/reads_${readlen}_${INSERT_MEAN}.read.info" "${DST}.read.info"
	        fi
	    done
    done

    echo > /dev/stderr
}


function simulate_adapter_id_reads()
{
    echo "------------------------------------------------------------" > /dev/stderr
    echo "Simulating reads for adapter identification ..."  > /dev/stderr
    echo ""  > /dev/stderr

    readlen=100
    for INSERT_MEAN in ${ADAPTER_ID_INSERT_SIZES[*]};
    do
        for run_n in ${REPLICATES};
        do
            DST="${SIMULATED_ADAPTER_ID_PREFIX}_${run_n}_${readlen}_${INSERT_MEAN}"

            if [ -e "${DST}.read.info" ];
            then
                echo "    Skipping ${DST}.*" > /dev/stderr
            else
                echo "    Simulating reads run=${run_n}, l=${readlen}, m=${INSERT_MEAN}, v=75 ..." > /dev/stderr
                rm -rf "${DST:?}/"
                mkdir -p "${DST}/"

                if ! ${EXEC_TIME} --verbose --output "${DST}/time" ${EXEC_PIRS} simulate -x "${SIMULATED_NREADS}" -l "${readlen}" -i "${REFSEQ}" -c 0 -m "${INSERT_MEAN}" -v 75 -Q 33 -o "${DST}/reads" \
                    > "${DST}/reads.stdout" 2> "${DST}/reads.stderr";
                then
                    echo "Error simulated reads ..." > /dev/stderr
                    exit 1
                fi

                gunzip "${DST}/reads_${readlen}_${INSERT_MEAN}.read.info.gz"
                ln -sf "$(basename "${DST}")/reads_${readlen}_${INSERT_MEAN}_1.fq" "${DST}_1.fq"
                ln -sf "$(basename "${DST}")/reads_${readlen}_${INSERT_MEAN}_2.fq" "${DST}_2.fq"
                ln -sf "$(basename "${DST}")/reads_${readlen}_${INSERT_MEAN}.read.info" "${DST}.read.info"
            fi
        done
    done

    echo > /dev/stderr
}



###############################################################################


function run_pear()
{
    NTHREADS=$1
    DST=$2
    FQ_INFO=$3
    FQ_MATE1=$4
    FQ_MATE2=$5
    shift 5

    do_run "${DST}" ${EXEC_PEAR} paired -j "${NTHREADS}" -t 0 \
        -f "${FQ_MATE1}" -r "${FQ_MATE2}" \
        -o "${DST}/reads" "$@"

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode PE --collapsed
}


function run_adapterremoval()
{
    MODE=$1
    EXECUTABLE=$2
    NTHREADS=$3
    DST=$4
    FQ_INFO=$5
    shift 5

    if [ "$NTHREADS" -eq 1 ];
    then
        NTHREADS_CLI=""
    else
        NTHREADS_CLI="--threads ${NTHREADS}"
    fi

    if [ "${MODE}" = "COLLAPSE" ];
    then
        MODE="PE"
        COLLAPSE_CLI="--collapsed"
    else
        COLLAPSE_CLI=""
    fi

    do_run "${DST}" "${EXECUTABLE}" --basename "${DST}/reads" ${NTHREADS_CLI} "$@"

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode "${MODE}" ${COLLAPSE_CLI}
}


function run_peat_se()
{
    NTHREADS=$1
    DST=$2
    FQ_INFO=$3
    FQ_MATE1=$4
    shift 4

    do_run "${DST}" ${EXEC_PEAT} single -q SANGER -n "${NTHREADS}" \
        -a "$(cat adapters/adapter_1.txt)" \
        -i "${FQ_MATE1}" -o "${DST}/reads.fastq" "$@"

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode SE
}


function run_peat_pe()
{
    NTHREADS=$1
    DST=$2
    FQ_INFO=$3
    FQ_MATE1=$4
    FQ_MATE2=$5
    shift 5

    do_run "${DST}" ${EXEC_PEAT} paired -n "${NTHREADS}" \
        -1 "${FQ_MATE1}" -2 "${FQ_MATE2}" \
        -o "${DST}/reads" "$@"

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode PE
}


function run_trimmomatic_se()
{
    MODE=$1
    NTHREADS=$2
    DST=$3
    FQ_INFO=$4
    FQ_MATE1=$5
    FQ_EXT=$6
    shift 6

    if [ "${MODE}" = "MIXED" ];
    then
        ADAPTERS="adapters/mixed.fasta"
    else
        ADAPTERS="adapters/adapters.fasta"
    fi

    do_run "${DST}" ${EXEC_JAVA} -jar ${JAR_TRIMMOMATIC} SE -phred33 -threads "${NTHREADS}" \
            "${FQ_MATE1}" "${DST}/reads${FQ_EXT}" "$@" \
            "ILLUMINACLIP:${ADAPTERS}:2:7:10"

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode SE \
        --dimers-are-discarded
}


function run_trimmomatic_pe()
{
    MODE=$1
    NTHREADS=$2
    DST=$3
    FQ_INFO=$4
    FQ_MATE1=$5
    FQ_MATE2=$6
    FQ_EXT=$7
    shift 7

    if [ "${MODE}" = "MIXED" ];
    then
        ADAPTERS="adapters/mixed.fasta"
    else
        ADAPTERS="adapters/adapters.fasta"
    fi

    do_run "${DST}" ${EXEC_JAVA} -jar ${JAR_TRIMMOMATIC} PE \
        -threads "${NTHREADS}" -phred33 \
        "${FQ_MATE1}" "${FQ_MATE2}" \
        "${DST}/reads.1${FQ_EXT}" "${DST}/reads.1.singletons.${FQ_EXT}" \
        "${DST}/reads.2${FQ_EXT}" "${DST}/reads.2.singletons.${FQ_EXT}" \
        "ILLUMINACLIP:${ADAPTERS}:2:30:7:1:true" "$@"

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode PE \
        --dimers-are-discarded
}



function run_flexbar()
{
    MODE=$1
    NTHREADS=$2
    DST=$3
    FQ_INFO=$4
    shift 4

    do_run "${DST}" ${EXEC_FLEXBAR} -m 0 -t "${DST}/reads" -n "${NTHREADS}" "$@"

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode "${MODE}"
}


function run_skewer()
{
    MODE=$1
    NTHREADS=$2
    DST=$3
    FQ_INFO=$4
    shift 4

    do_run "${DST}" ${EXEC_SKEWER} -o "${DST}/reads" -l 0 -t "${NTHREADS}" "$@" \

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode "${MODE}"
}


function run_scythe()
{
    DST=$1
    FQ_INFO=$2
    FQ_MATE1=$3
    shift 3

    do_run "${DST}" ${EXEC_SCYTHE} -M 0 "$@" \
                -a "adapters/adapter_1.fasta" \
                -o "${DST}/reads.trimmed" \
                "${FQ_MATE1}"

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode SE
}


function run_alientrimmer_se()
{
    DST=$1
    FQ_INFO=$2
    FQ_MATE1=$3
    shift 3

    do_run "${DST}" ${EXEC_ALIENTRIMMER} -l 0 "$@" \
        -i "${FQ_MATE1}" \
        -o "${DST}/reads.trimmed"

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode SE
}


function run_alientrimmer_pe()
{
    DST=$1
    FQ_INFO=$2
    FQ_MATE1=$3
    FQ_MATE2=$4
    shift 4

    do_run "${DST}" ${EXEC_ALIENTRIMMER} -l 0 "$@" \
        -if "${FQ_MATE1}" \
        -ir "${FQ_MATE2}" \
        -os "${DST}/reads.trimmed" \
        -or "${DST}/reads.trimmed.mate1" \
        -of "${DST}/reads.trimmed.mate2"

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode PE
}


function run_cutadapt_se()
{
    DST=$1
    FQ_INFO=$2
    FQ_MATE1=$3
    shift 3

    do_run "${DST}" ${EXEC_CUTADAPT} "$@" \
        -o "${DST}/reads.trimmed" \
        "${FQ_MATE1}"

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode SE
}


function run_cutadapt_pe()
{
    DST=$1
    FQ_INFO=$2
    FQ_MATE1=$3
    FQ_MATE2=$4
    shift 4

    do_run "${DST}" ${EXEC_CUTADAPT} "$@" \
        -o "${DST}/reads.trimmed.1" \
        -p "${DST}/reads.trimmed.2" \
        "${FQ_MATE1}" \
        "${FQ_MATE2}"

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode PE
}


function run_leeHom()
{
    MODE=$1
    DST=$2
    FQ_INFO=$3
    shift 3

    if [ "${MODE}" = "COLLAPSE" ];
    then
        MODE="PE"
        COLLAPSE_CLI="--collapsed"
    else
        COLLAPSE_CLI=""
    fi

    do_run "${DST}" "${EXEC_LEEHOM}" "$@" -fqo "${DST}/reads"

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode "${MODE}" ${COLLAPSE_CLI}
}


function run_fastq_mcf_se()
{
    DST=$1
    FQ_INFO=$2
    ADAPTERS=$3
    FQ_MATE1=$4
    shift 4

    # -0 to disable quality trimming
    do_run "${DST}" ${EXEC_FASTQ_MCF} -0 \
        "${ADAPTERS}" \
        "${FQ_MATE1}" \
        -o "${DST}/reads.trimmed" "$@"

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode SE \
        --dimers-are-discarded
}


function run_fastq_mcf_pe()
{
    DST=$1
    FQ_INFO=$2
    ADAPTERS=$3
    FQ_MATE1=$4
    FQ_MATE2=$5
    shift 5

    # -0 to disable quality trimming
    do_run "${DST}" ${EXEC_FASTQ_MCF} -0 \
        "${ADAPTERS}" \
        "${FQ_MATE1}" \
        "${FQ_MATE2}" \
        -o "${DST}/reads.trimmed.mate1" \
        -o "${DST}/reads.trimmed.mate2" \
        "$@"

    do_evaluate "${DST}" "${FQ_INFO}" --read-mode PE \
        --dimers-are-discarded
}


###############################################################################

function run_minion()
{
    DST=$1
    FQ_MATE1=$2
    FQ_MATE2=$3
    shift 3

    if [ ! -e "${DST}/DONE" ];
    then
        echo "Runing ${DST} ..." > /dev/stderr
        rm -rf "${DST:?}/"
        mkdir -p "${DST}/"

        ${EXEC_MINION} search-adapter -show 5 -i "${FQ_MATE1}" \
            > "${DST}/mate1.txt"
        ${EXEC_MINION} search-adapter -show 5 -i "${FQ_MATE2}" \
            > "${DST}/mate2.txt"

        touch "${DST}/DONE"
    fi
}


function run_adapterremoval_id()
{
    DST=$1
    FQ_MATE1=$2
    FQ_MATE2=$3
    shift 3

    if [ ! -e "${DST}/DONE" ];
    then
        echo "Runing ${DST} ..." > /dev/stderr
        rm -rf "${DST:?}/"
        mkdir -p "${DST}/"

        ${EXEC_ADAPTERREMOVAL2x} --identify-adapters \
            --file1 "${FQ_MATE1}" --file2 "${FQ_MATE2}" \
            > "${DST}/mates.txt"

        touch "${DST}/DONE"
    fi
}


###############################################################################

# Common replicates strings
REPLICATES=$(printf "%03i\n" $(seq 1 ${NUM_REPLICATES}))


function benchmark_se()
{
    simulate_reads

    echo "------------------------------------------------------------" > /dev/stderr
    echo "Running SE benchmarks ..."  > /dev/stderr
    echo ""  > /dev/stderr

    for readlen in ${READ_LENGTHS[*]};
    do
        for run_n in ${REPLICATES};
        do
            { # Shuffle each individual run
                SIMULATED_MATE1="${SIMULATED_PREFIX}_${run_n}_${readlen}_1.fq"
                SIMULATED_MATE2="${SIMULATED_PREFIX}_${run_n}_${readlen}_2.fq"
                SIMULATED_INFO="${SIMULATED_PREFIX}_${run_n}_${readlen}.read.info"
                RESULTS="results/se/${readlen}_${run_n}"

                # -mm 3 corresponds to AR 2.x defaults
                echo run_adapterremoval SE ${EXEC_ADAPTERREMOVAL1x} 1 \
                    "${RESULTS}/adapterremoval1x_mm3" \
                    "${SIMULATED_INFO}" \
                    --file1 "${SIMULATED_MATE1}" \
                    --mm 3

                # -mm 3 --minadapteroverlap 3 (test)
                echo run_adapterremoval SE ${EXEC_ADAPTERREMOVAL2x} 1 \
                    "${RESULTS}/adapterremoval2x_min3_mm3" \
                    "${SIMULATED_INFO}" \
                    --file1 "${SIMULATED_MATE1}" \
                    --minadapteroverlap 3 --mm 3

                # -mm 5 --minadapteroverlap 3 (test)
                echo run_adapterremoval SE ${EXEC_ADAPTERREMOVAL2x} 1 \
                    "${RESULTS}/adapterremoval2x_min3_mm5" \
                    "${SIMULATED_INFO}" \
                    --file1 "${SIMULATED_MATE1}" \
                    --minadapteroverlap 3 --mm 5

                for nthreads in $(seq 1 ${MAX_THREADS});
                do
                    AR_PREFIX="${RESULTS}/adapterremoval2x"
                    echo run_adapterremoval SE ${EXEC_ADAPTERREMOVAL2x} "${nthreads}" \
                        "${AR_PREFIX}_t${nthreads}" \
                        "${SIMULATED_INFO}" \
                        --file1 "${SIMULATED_MATE1}" \

                    echo run_skewer SE "${nthreads}" \
                        "${RESULTS}/skewer_t${nthreads}" \
                        "${SIMULATED_INFO}" \
                        "${SIMULATED_MATE1}"

                    echo run_flexbar SE "${nthreads}" \
                        "${RESULTS}/flexbar_t${nthreads}" \
                        "${SIMULATED_INFO}" \
                        -a "./adapters/adapter_1.fasta" \
                        -r "${SIMULATED_MATE1}"

                    echo run_trimmomatic_se SE "${nthreads}" \
                        "${RESULTS}/trimmomatic_t${nthreads}" \
                        "${SIMULATED_INFO}" \
                        "${SIMULATED_MATE1}" \
                        ".fastq"

                    echo run_peat_se "${nthreads}" \
                        "${RESULTS}/peat_t${nthreads}" \
                        "${SIMULATED_INFO}" \
                        "${SIMULATED_MATE1}"
                done

                echo run_cutadapt_se "${RESULTS}/cutadapt" \
                    "${SIMULATED_INFO}" \
                    "${SIMULATED_MATE1}" \
                    -a "$(cat adapters/adapter_1.txt)"

                echo run_scythe "${RESULTS}/scythe" \
                    "${SIMULATED_INFO}" \
                    "${SIMULATED_MATE1}"

                echo run_alientrimmer_se "${RESULTS}/alientrimmer_q00" \
                    "${SIMULATED_INFO}" \
                    "${SIMULATED_MATE1}" \
                    -c "adapters/adapter_1.txt" \
                    -q 0

                echo run_leeHom SE \
                    "${RESULTS}/leeHom" \
                    "${SIMULATED_INFO}" \
                    -fq1 "${SIMULATED_MATE1}"

                echo run_leeHom SE \
                    "${RESULTS}/leeHom_ancient" \
                    "${SIMULATED_INFO}" \
                    -fq1 "${SIMULATED_MATE1}" \
                    --ancientdna

                echo run_fastq_mcf_se "${RESULTS}/fastq_mcf" \
                    "${SIMULATED_INFO}" \
                    "./adapters/adapter_1.fasta" \
                    "${SIMULATED_MATE1}"
            } | shuffle_and_run
        done
    done

    ${SCRIPT_MERGE} results/se > results/se.table
}


function benchmark_pe()
{
    simulate_reads

    echo "------------------------------------------------------------" > /dev/stderr
    echo "Running PE benchmarks ..."  > /dev/stderr
    echo ""  > /dev/stderr

    for readlen in ${READ_LENGTHS[*]};
    do
        for run_n in ${REPLICATES};
        do
            { # Shuffle each individual run
                SIMULATED_MATE1="${SIMULATED_PREFIX}_${run_n}_${readlen}_1.fq"
                SIMULATED_MATE2="${SIMULATED_PREFIX}_${run_n}_${readlen}_2.fq"
                SIMULATED_INFO="${SIMULATED_PREFIX}_${run_n}_${readlen}.read.info"
                RESULTS="results/pe/${readlen}_${run_n}"

                DEFAULT_ARGS=("${SIMULATED_INFO}" "${SIMULATED_MATE1}" "${SIMULATED_MATE2}")

                # -mm 3 corresponds to AR 2.x defaults
                echo run_adapterremoval PE ${EXEC_ADAPTERREMOVAL1x} 1 \
                    "${RESULTS}/adapterremoval1x_mm3" \
                    "${SIMULATED_INFO}" \
                    --file1 "${SIMULATED_MATE1}" \
                    --file2 "${SIMULATED_MATE2}" \
                    --mm 3

                for nthreads in $(seq 1 ${MAX_THREADS});
                do
                    AR_PREFIX="${RESULTS}/adapterremoval2x"

                    echo run_adapterremoval PE ${EXEC_ADAPTERREMOVAL2x} "${nthreads}" \
                        "${AR_PREFIX}_t${nthreads}" \
                        "${SIMULATED_INFO}" \
                        --file1 "${SIMULATED_MATE1}" \
                        --file2 "${SIMULATED_MATE2}" \

                    echo run_skewer PE "${nthreads}" \
                        "${RESULTS}/skewer_t${nthreads}" \
                        "${DEFAULT_ARGS[*]}"

                    echo run_flexbar PE "${nthreads}" \
                        "${RESULTS}/flexbar_t${nthreads}" \
                        "${SIMULATED_INFO}" \
                        -a "./adapters/adapters.fasta" \
                        -r "${SIMULATED_MATE1}" \
                        -p "${SIMULATED_MATE2}"

                    echo run_trimmomatic_pe PE "${nthreads}" \
                        "${RESULTS}/trimmomatic_t${nthreads}" \
                        "${DEFAULT_ARGS[*]}" \
                        ".fastq"

                    echo run_peat_pe "${nthreads}" \
                        "${RESULTS}/peat_t${nthreads}" \
                        "${DEFAULT_ARGS[*]}"
                done

                echo run_cutadapt_pe "${RESULTS}/cutadapt" \
                    "${DEFAULT_ARGS[*]}" \
                    -a "$(cat adapters/adapter_1.txt)" \
                    -A "$(cat adapters/adapter_2.txt)"

                echo run_alientrimmer_pe "${RESULTS}/alientrimmer_q00" \
                    "${DEFAULT_ARGS[*]}" \
                    -cf "adapters/adapter_1.txt" \
                    -cr "adapters/adapter_2.txt" \
                    -q 0

                echo run_leeHom PE \
                    "${RESULTS}/leeHom" \
                    "${SIMULATED_INFO}" \
                    -fq1 "${SIMULATED_MATE1}" \
                    -fq2 "${SIMULATED_MATE2}"

                echo run_leeHom PE \
                    "${RESULTS}/leeHom_ancient" \
                    "${SIMULATED_INFO}" \
                    -fq1 "${SIMULATED_MATE1}" \
                    -fq2 "${SIMULATED_MATE2}" \
                    --ancientdna

                echo run_fastq_mcf_pe "${RESULTS}/fastq_mcf" \
                    "${SIMULATED_INFO}" \
                    "./adapters/adapters.fasta" \
                    "${SIMULATED_MATE1}" \
                    "${SIMULATED_MATE2}"
            } | shuffle_and_run
        done
    done

    ${SCRIPT_MERGE} results/pe > results/pe.table
}


function benchmark_collapse()
{
    simulate_reads

    echo "------------------------------------------------------------" > /dev/stderr
    echo "Running collapsing benchmarks ..."  > /dev/stderr
    echo ""  > /dev/stderr

    for readlen in ${READ_LENGTHS[*]};
    do
        for run_n in ${REPLICATES};
        do
            { # Shuffle each individual run
                SIMULATED_MATE1="${SIMULATED_PREFIX}_${run_n}_${readlen}_1.fq"
                SIMULATED_MATE2="${SIMULATED_PREFIX}_${run_n}_${readlen}_2.fq"
                SIMULATED_INFO="${SIMULATED_PREFIX}_${run_n}_${readlen}.read.info"
                RESULTS="results/collapse/${readlen}_${run_n}"

                DEFAULT_ARGS=("${SIMULATED_INFO}" "${SIMULATED_MATE1}" "${SIMULATED_MATE2}")

                # -mm 3 corresponds to AR 2.x defaults
                echo run_adapterremoval COLLAPSE ${EXEC_ADAPTERREMOVAL1x} 1 \
                    "${RESULTS}/adapterremoval1x_mm3" \
                    "${SIMULATED_INFO}" \
                    --file1 "${SIMULATED_MATE1}" \
                    --file2 "${SIMULATED_MATE2}" \
                    --mm 3 --collapse

                for nthreads in $(seq 1 ${MAX_THREADS});
                do
                    AR_PREFIX="${RESULTS}/adapterremoval2x"

                    echo run_adapterremoval COLLAPSE ${EXEC_ADAPTERREMOVAL2x} "${nthreads}" \
                        "${AR_PREFIX}_t${nthreads}" \
                        "${SIMULATED_INFO}" \
                        --file1 "${SIMULATED_MATE1}" \
                        --file2 "${SIMULATED_MATE2}" \
                        --collapse

                    echo run_pear "${nthreads}" \
                        "${RESULTS}/pear_t${nthreads}" \
                        "${SIMULATED_INFO}" \
                        "${SIMULATED_MATE1}" \
                        "${SIMULATED_MATE2}"
                done

                echo run_leeHom COLLAPSE \
                    "${RESULTS}/leeHom" \
                    "${SIMULATED_INFO}" \
                    -fq1 "${SIMULATED_MATE1}" \
                    -fq2 "${SIMULATED_MATE2}"

                echo run_leeHom COLLAPSE \
                    "${RESULTS}/leeHom_ancient" \
                    "${SIMULATED_INFO}" \
                    -fq1 "${SIMULATED_MATE1}" \
                    -fq2 "${SIMULATED_MATE2}" \
                    --ancientdna
            } | shuffle_and_run
        done
    done

    ${SCRIPT_MERGE} results/collapse > results/collapse.table
}



function benchmark_mixed_se
{
    simulate_mixed_reads

    echo "------------------------------------------------------------" > /dev/stderr
    echo "Running mixed se adapters benchmarks ..."  > /dev/stderr
    echo ""  > /dev/stderr

    for readlen in ${READ_LENGTHS[*]};
    do
	    for run_n in ${REPLICATES};
	    do
	        { # Shuffle each individual run
	            SIMULATED_MATE1="${SIMULATED_MIXED_PREFIX}_${run_n}_${readlen}_1.fq"
	            SIMULATED_MATE2="${SIMULATED_MIXED_PREFIX}_${run_n}_${readlen}_2.fq"
	            SIMULATED_INFO="${SIMULATED_MIXED_PREFIX}_${run_n}_${readlen}.read.info"
	            RESULTS="results/mixed_se/${readlen}_${run_n}"

	            DEFAULT_ARGS=("${SIMULATED_INFO}" "${SIMULATED_MATE1}" "${SIMULATED_MATE2}")
	            AR_PREFIX="${RESULTS}/adapterremoval2x"

	            # -mm 3 --minadapteroverlap 3 (test)
	            echo run_adapterremoval SE ${EXEC_ADAPTERREMOVAL2x} 1 \
	                "${AR_PREFIX}_min3_mm3" \
	                "${SIMULATED_INFO}" \
	                --file1 "${SIMULATED_MATE1}" \
	                --adapter-list ./adapters/mixed.table \
	                --minadapteroverlap 3 --mm 3

	            # -mm 5 --minadapteroverlap 3 (test)
	            echo run_adapterremoval SE ${EXEC_ADAPTERREMOVAL2x} 1 \
	                "${AR_PREFIX}_min3_mm5" \
	                "${SIMULATED_INFO}" \
	                --file1 "${SIMULATED_MATE1}" \
	                --adapter-list ./adapters/mixed.table \
	                --minadapteroverlap 3 --mm 5

	            for nthreads in $(seq 1 ${MAX_THREADS});
	            do
	                echo run_adapterremoval SE ${EXEC_ADAPTERREMOVAL2x} "${nthreads}" \
	                    "${AR_PREFIX}_t${nthreads}" \
	                    "${SIMULATED_INFO}" \
	                    --file1 "${SIMULATED_MATE1}" \
	                    --adapter-list ./adapters/mixed.table

	                echo run_trimmomatic_se MIXED "${nthreads}" \
	                    "${RESULTS}/trimmomatic_t${nthreads}" \
	                    "${SIMULATED_INFO}" \
	                    "${SIMULATED_MATE1}" \
	                    ".fastq"
	            done

	            echo run_alientrimmer_se "${RESULTS}/alientrimmer_q00" \
	                "${SIMULATED_INFO}" \
	                "${SIMULATED_MATE1}" \
	                -c "adapters/mixed_1.txt" \
	                -q 0

	            echo run_cutadapt_se "${RESULTS}/cutadapt" \
	                "${SIMULATED_INFO}" \
	                "${SIMULATED_MATE1}" \
	                $(for seq in $(cat adapters/mixed_1.txt); do echo "-a ${seq}";done)

	            echo run_fastq_mcf_se "${RESULTS}/fastq_mcf" \
	                "${SIMULATED_INFO}" \
	                "./adapters/mixed_1.fasta" \
	                "${SIMULATED_MATE1}"
	        } | shuffle_and_run
	    done
	done

    ${SCRIPT_MERGE} results/mixed_se > results/mixed_se.table
}


function benchmark_mixed_pe
{
    simulate_mixed_reads

    echo "------------------------------------------------------------" > /dev/stderr
    echo "Running mixed pe adapters benchmarks ..."  > /dev/stderr
    echo ""  > /dev/stderr

    for readlen in ${READ_LENGTHS[*]};
    do
	    for run_n in ${REPLICATES};
	    do
	        { # Shuffle each individual run
	            SIMULATED_MATE1="${SIMULATED_MIXED_PREFIX}_${run_n}_${readlen}_1.fq"
	            SIMULATED_MATE2="${SIMULATED_MIXED_PREFIX}_${run_n}_${readlen}_2.fq"
	            SIMULATED_INFO="${SIMULATED_MIXED_PREFIX}_${run_n}_${readlen}.read.info"
	            RESULTS="results/mixed_pe/${readlen}_${run_n}"

	            DEFAULT_ARGS=("${SIMULATED_INFO}" "${SIMULATED_MATE1}" "${SIMULATED_MATE2}")

	            for nthreads in $(seq 1 ${MAX_THREADS});
	            do
	                AR_PREFIX="${RESULTS}/adapterremoval2x"

	                echo run_adapterremoval PE ${EXEC_ADAPTERREMOVAL2x} "${nthreads}" \
	                    "${AR_PREFIX}_t${nthreads}" \
	                    "${SIMULATED_INFO}" \
	                    --file1 "${SIMULATED_MATE1}" \
	                    --file2 "${SIMULATED_MATE2}" \
	                    --adapter-list ./adapters/mixed.table

	                echo run_peat_pe "${nthreads}" \
	                    "${RESULTS}/peat_t${nthreads}" \
	                    "${DEFAULT_ARGS[*]}"

	                echo run_trimmomatic_pe MIXED "${nthreads}" \
	                    "${RESULTS}/trimmomatic_t${nthreads}" \
	                    "${DEFAULT_ARGS[*]}" \
	                    ".fastq"
	            done

	            echo run_alientrimmer_pe "${RESULTS}/alientrimmer_q00" \
	                "${DEFAULT_ARGS[*]}" \
	                -cf "adapters/mixed_1.txt" \
	                -cr "adapters/mixed_2.txt" \
	                -q 0

	            echo run_cutadapt_pe "${RESULTS}/cutadapt" \
	                "${DEFAULT_ARGS[*]}" \
	                $(for seq in $(cat adapters/mixed_1.txt); do echo "-a ${seq}";done) \
	                $(for seq in $(cat adapters/mixed_2.txt); do echo "-A ${seq}";done)

	            echo run_fastq_mcf_pe "${RESULTS}/fastq_mcf" \
	                "${SIMULATED_INFO}" \
	                "./adapters/mixed.fasta" \
	                "${SIMULATED_MATE1}" \
	                "${SIMULATED_MATE2}"
	        } | shuffle_and_run
	    done
	done

    ${SCRIPT_MERGE} results/mixed_pe > results/mixed_pe.table
}


function benchmark_adapter_id
{
    simulate_adapter_id_reads

    echo "------------------------------------------------------------" > /dev/stderr
    echo "Running adapter id benchmarks ..."  > /dev/stderr
    echo ""  > /dev/stderr

    readlen=100
    for INSERT_MEAN in ${ADAPTER_ID_INSERT_SIZES[*]};
    do
        for run_n in ${REPLICATES};
        do
            SIMULATED_MATE1="${SIMULATED_ADAPTER_ID_PREFIX}_${run_n}_${readlen}_${INSERT_MEAN}_1.fq"
            SIMULATED_MATE2="${SIMULATED_ADAPTER_ID_PREFIX}_${run_n}_${readlen}_${INSERT_MEAN}_2.fq"
            SIMULATED_INFO="${SIMULATED_ADAPTER_ID_PREFIX}_${run_n}_${readlen}_${INSERT_MEAN}.read.info"
            RESULTS="results/adapter_id/${readlen}_${INSERT_MEAN}_${run_n}"

            run_minion "${RESULTS}/minion" \
                "${SIMULATED_MATE1}" \
                "${SIMULATED_MATE2}"

            run_adapterremoval_id "${RESULTS}/adapterremovalv2" \
                "${SIMULATED_MATE1}" \
                "${SIMULATED_MATE2}"
        done
    done

    ${SCRIPT_EVALUATE_ID} results/adapter_id > results/adapter_id.table
}



cd "$(dirname "$0")"
mkdir -p results

fetch_reference

echo > /dev/stderr

if [ "$* " = "all " ];
then
    benchmark_se
    benchmark_pe
    benchmark_collapse
    benchmark_adapter_id
    benchmark_mixed_pe
    benchmark_mixed_se
elif [ "$* " = "se " ];
then
    benchmark_se
elif [ "$* " = "pe " ];
then
    benchmark_pe
elif [ "$* " = "collapse " ];
then
    benchmark_collapse
elif [ "$* " = "adapter_id " ];
then
    benchmark_adapter_id
elif [ "$* " = "mixed " ];
then
    benchmark_mixed_pe
    benchmark_mixed_se
else
    echo "Usage: benchmark.sh <command>" > /dev/stderr
    echo "Commands: all se pe collapse mixed adapter_id" > /dev/stderr
fi
