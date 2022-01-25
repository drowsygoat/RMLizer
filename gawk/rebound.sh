
#!/bin/bash
#.
#.                         |__  o\
#.                         | W    \O
#.      ###############    |       |\_           |\
#.      rebound.sh v2.0    |      /-\            \O
#.      ###############    |    /     \           |
#.                         |                     /|
#.                         |                    |  \
#.
#. Usage:
#.        rebound.sh [OPTIONS] [BAM FILES]
#.
#. Prerequisites:
#.      Samtools
#.      GAWK
#.
#. Required arguments:
#. -r     Fasta reference file. This should be the file that was used to generate
#.        the BAM file(s) and have .fa or .fasta extension. Alternatively use
#.        SNP-updated reference (in such case, [-f] flag must be used, see manual).
#. -g     Annotation (GTF/GFF) file with path. If using unstranded libraries, this file
#.        may include strandedness information appended to gene identifiers as double
#.        underscore followed by + or - sign, e.g., ENSMUSG00000079037__+.
#. BAM    BAM files to process (wildcards [*] accepted).
#.
#. Optional arguments:
#. -a     Path to rebound.awk script (./).
#. -o     Output directory (./).
#. -l     <STAR/BBmap> Aligner used to generate BAM files (STAR).
#. -f     Force addition/replacement of MD tags. BAM file(s) will be done automatically
#.        scanned and MD tags will be added it not found.
#. -s     <forward/reverse> Data strandedness (reverse). SLAM-seq requires
#.        stranded reads.
#. -i     <ban|1-30|disabled> Removes reads with (D)eletions egual or greater than this value.
#.        "Ban" removes all reads with (D)eletions, i.e., those containing D in CIGAR strings.
#. -q     <INT,INT,FLOAT,INT,INT> Comma-separated read filter settings:
#.        1. <0-255, INT> MAPQ/mapping quality filter (255 for STAR, 20 for BBmap). Keeps
#.        reads with MAPQ equal or higher than this value. Please note that 255 is STAR's
#.        default value assigned to unique mappers.
#.        2. <INT> Read filtering threshold by edit distance. (20 for STAR, 100000 for BBmap)
#.        Reads with at most this value are preserved. If i(N)trons are marked as as (D)eletions
#.        in CIGARs, e.g., as by BBmap, then setting this value to 100000 disables the filter.
#.        3. <0-1, FLOAT> Read filtering threshold by mismatch fraction (0.90). Reads with at
#.        least this value are kept.
#.        4. <1-41, INT> Ignore mismatches with called base quality below this number (27).
#.        5. <1-41, INT> Ignore TC mismatches with called base quality below this number (27).
#. -p     <0-1, FLOAT> Proportion of reads to be randomly sampled from BAM file(s), useful
#.        for shortening processing times for testing parameters.
#. -t     <INT> Integer specifying the number of CPU threads to use (1). This setting only
#.        affects samtools.
#. -u     Count only uniquely mapped reads.
#. -h     Displays this help while ignoring any remaining flags.

# -x     <STR> Per thread memory block size for samtools sort (1G).

NORM=$(tput sgr0)
BOLD=$(tput bold)
REV=$(tput smso)

help() {
grep "#\." $0 | sed 's/#\.//'
exit 0
}
REBOUND_COMMAND="rebound_lite.sh $@"

usage() {
   echo "Usage:
   ${BOLD}rebound.sh [options] [file1.bam file2.bam ...]${NORM}
   Type ${BOLD}rebound.sh -h${NORM} to get more help."
}
exit_err() {
  usage
  exit 1
}
[ $# -eq 0 ] && exit_err
# Set some defaults
UNIQ_ARG=1
THR_ARG=1
MD_ARG=0
PROCENT="1.00"
MAPPER_ARG="STAR"
MEM_ARG="1G"
MODE_ARG="normal"
STRANDEDNESS_ARG="reverse"
OUTDIR_ARG="./"
AWK_FILE="./rebound_lite.awk"
while [ $# -gt 0 ]
do
    unset OPTIND
    unset OPTARG
    while getopts ":a:o:l:r:fi:g:m:q:s:p:ux:t:h" opcje
    do
    case $opcje in
        a) # Specify a path to AWK script.
            AWK_FILE=${OPTARG}
            ;;
        o) # Specify a destination folder for the results.
            OUTDIR_ARG=${OPTARG}
            ;;
        l) # Specify a destination folder for the results.
            MAPPER_ARG=${OPTARG}
            ;;
        r) # Specify a valid reference file with full path.
            REF_ARG=${OPTARG}
            ;;
        f) # Force MD.
            MD_ARG=1
            ;;
        i) # Ban reads with indels.
            INDEL_BAN_ARG=${OPTARG}
            ;;
        g) # Specify a valid annotation (gtf) file with  path.
            GTF_ARG=${OPTARG}
            ;;
        m) # Mode of function.
            MODE_ARG=${OPTARG}
            ;;
        q) # Quality filtering settings (comma-separated values).
            QUALITY_ARG=${OPTARG}
            ;;
        s) # Strandedness.
            STRANDEDNESS_ARG=${OPTARG}
            ;;
        p) # Procent.
            PROCENT=${OPTARG}
            ;;
        u) # Include only unique reads.
            UNIQ_ARG=1
            ;;
        x) # Memory for samtools sort in G.
            MEM_ARG=${OPTARG}
            ;;
        t) # Number of threads for samtools.
            THR_ARG=${OPTARG}
            ;;
        h) # Display help.
            help
            ;;
        \?) #unrecognized option - show help
            echo -e \\n"Flag -${BOLD}$OPTARG${NORM} not allowed."
            usage
    esac
    done
    shift "$((OPTIND-1))"
    ARGS="${ARGS} $1"
    shift
done

BAM_LIST=$ARGS

##### check for software #####

check_soft() {
    if ! [ -f "$AWK_FILE" ] ; then
        echo "Fatal error: $AWK_FILE not found."
        exit 1
    else
        echo "Using script: $AWK_FILE"
    fi

    if ! command -v gawk &> /dev/null ; then
        echo "Fatal error. GAWK could not be found."
        exit 1
    else
        echo Found $(gawk -V | grep -E "wk.*\d*\.\d*" | cut -d, -f1)
    fi

    if ! command -v samtools &> /dev/null ; then
        echo "Fatal error. Samtools could not be found."
        exit 1
    else
        echo Found Samtools $(samtools --version | head -1 | cut -d' ' -f2)
        if [[ $MODE_ARG == "normal" ]] ; then
            if ! command -v htseq-count --version &> /dev/null ; then
                echo "Fatal error: HTseq could not be found."
                exit 1
            else
                echo "Found HTseq" $(htseq-count --version)
                # echo "Found HTseq version" $(htseq-count --help | grep -o " \d\..*")
            fi
        fi
    fi
}

check_soft

printf "CPU threads: ${THR_ARG}\n"

# printf "RAM per thread: ${MEM_ARG}\n"

##### parameter check #####
if ! [[ $MODE_ARG =~ ^.*nocount.*$|^.*normal.*$|^.*htseq_only.*$  ]] ; then
    echo "Fatal error: Mode must be set to \"normal\"" # or \"htseq_only\""
    exit_err

elif [[ $MODE_ARG =~ ^.*htseq_only.*$  ]] ; then
    echo "Mode: htseq_only"

elif [[ $MODE_ARG =~ ^.*norm.*$  ]] ; then
    echo "Mode: normal"

    if [[ $STRANDEDNESS_ARG =~ reverse|forward|unstranded ]] ; then # change if statement to case statement
        echo "Strandedness is set to $STRANDEDNESS_ARG."

        if [[ $STRANDEDNESS_ARG =~ forward ]] ; then
              STRANDEDNESS_ARG="yes"
        elif [[ $STRANDEDNESS_ARG =~ reverse ]] ; then
              STRANDEDNESS_ARG="reverse"
        elif [[ $STRANDEDNESS_ARG =~ unstranded ]] ; then
              STRANDEDNESS_ARG="no"
        fi

    else

        echo "Fatal error: Strandedness must be one of \"forward\", \"reverse\" or \"unstranded\"."
        exit_err

    fi
fi

if [[ $OUTDIR_ARG == "." ]] ; then
    echo "Output will be saved in the current directory"
else
    mkdir -p $OUTDIR_ARG
    mkdir -p $OUTDIR_ARG/stderr
fi

if [[ $MAPPER_ARG =~ ^.*BBmap.*$ ]] ; then
    echo "Aligner: BBmap"
elif [[ $MAPPER_ARG =~ ^.*STAR.*$ ]] ; then
    echo "Aligner: STAR"
fi

if [[ $MD_ARG -eq 1 ]] ; then
    echo "Forcing MD tags"
fi

if [[ $UNIQ_ARG == 1 ]] ; then
   UNIQ_REP=YES
else
   UNIQ_REP=NO
fi

samtools view -@ $((THR_ARG-1)) $(echo $BAM_LIST | cut -d" " -f1) | head -1000 | gawk '{if ($6 ~ /N/){exit 1}}'

if [[ $? != 1 && $INDEL_BAN_ARG =~ ^.*ban.*$ ]] ; then
    echo "${BOLD}Warning${NORM}: provided BAM file(s) seem not to discriminate i(N)trons from (D)eletions and [ -i ] flag is set to \"ban\". Current command will effectively remove all reads that span splice junctions. In case this is not expected behavior, set the [ -i ] parameter value in range of [1-30] or \"disabled\" instead."
fi

# indel filter

if ! [[ -z $INDEL_BAN_ARG ]] ; then
    if ! [[ ( $INDEL_BAN_ARG -gt 0 && $INDEL_BAN_ARG -lt 31 ) || $INDEL_BAN_ARG =~ ^.*disabled|ban.*$ ]] ; then

        echo "Indel threshold must be in range [1-30], \"disabled\" (default) or \"ban\"."
        exit_err

    else

        INDEL_BAN_VAL=$INDEL_BAN_ARG # for later use

        if [[ $INDEL_BAN_ARG =~ ban ]] ; then
            INDEL_BAN_ARG="[DI]"
        elif [[ $INDEL_BAN_ARG -lt 10 ]] ; then
            INDEL_BAN_ARG=$(echo "[1-$INDEL_BAN_ARG][DI]")
        elif [[ $INDEL_BAN_ARG -lt 20 ]] ; then
            INDEL_BAN_ARG=$(echo "[^1-9][1-9][DI]|^[1-9][DI]|[1][0-$((INDEL_BAN_ARG-10))][DI]")
        elif [[ $INDEL_BAN_ARG -lt 30 ]] ; then
            INDEL_BAN_ARG=$(echo "[^1-9][1-9][DI]|^[1-9][DI]|[1-2][0-$((INDEL_BAN_ARG-20))][DI]")
        else
            INDEL_BAN_ARG=$(echo "[^1-9][1-9][DI]|^[1-9][DI]|[1-2][0-9][DI]|30[DI]")
        fi

    fi
fi

if [[ $BAM_LIST == " " ]] ; then

     echo "No BAM file(s) provided!"
     exit_err

else

    for BAM in $BAM_LIST
    do
        if ! [[ $BAM =~ ^.+\.bam$ || $BAM =~ ^.+\.sam$ ]] ; then # add verification by samtools
            echo $BAM "is not a valid input file name!"
            exit_err
        fi

        if [[ $MD_ARG -eq 0 ]] ; then

            samtools view -@$((THR_ARG-1)) $BAM | head -1 | gawk '{if($0~/^.*MD:Z:.*$/){exit 0}}'

            if ! [ $? -eq 0 ]; then

                echo "MD tags not found. Will be added."
                MD_ARG=1

                if [[ $(samtools view -H $BAM | head -1) =~ ^.*oordinate.*$ ]]; then

                    echo $BAM "is not sorted by chromosomal coordinates. Please sort and strat from scratch."
                    exit_err

                fi
            fi
        fi
    done
fi

if ! [ -f "$REF_ARG" ] ; then
    echo "Fatal error: $REF_ARG not found."
    exit_err
elif ! [[ $REF_ARG =~ ^.+\.fasta$ || $REF_ARG =~ ^.+\.fa$ ]] ; then
    echo "Fatal error: $REF_ARG does not look like a valid reference file."
    exit_err
else
    echo "Genome reference: $REF_ARG"
fi

if ! [ -f "$GTF_ARG" ] ; then
    echo "Fatal error: $GTF_ARG not found."
    exit_err
elif ! [[ $GTF_ARG =~ ^.+\.gtf$ || $REF_ARG =~ ^.+\.gff$ ]] ; then
    echo "Fatal error: $GTF_ARG does not look like a valid annotation file."
    exit_err
else
    echo "Annotations: $GTF_ARG"
fi

if ! [[ $MODE_ARG =~ ^.*htseq_only.*$ ]] ; then
    if ! [[ -z $QUALITY_ARG ]] ; then
        MAPQ_FILTER=$(echo $QUALITY_ARG | cut -d"," -f1)
        NM_FILTER=$(echo $QUALITY_ARG | cut -d"," -f2)
        MISMATCH_FILTER=$(echo $QUALITY_ARG | cut -d"," -f3)
        MISMATCH_QUALITY=$(echo $QUALITY_ARG | cut -d"," -f4)
        MISMATCH_QUALITY_TC=$(echo $QUALITY_ARG | cut -d"," -f5)

    elif [[ $MAPPER_ARG =~ ^.*STAR|star.*$ ]] ; then # defaults STAR
        MAPQ_FILTER=254
        NM_FILTER=20
        MISMATCH_FILTER="0.90"
        MISMATCH_QUALITY=27
        MISMATCH_QUALITY_TC=27

    elif [[ $MAPPER_ARG =~ ^.*BBmap|bbmap.*$ ]] ; then # defaults BBmap
        MAPQ_FILTER=20
        NM_FILTER=100000 # == off
        MISMATCH_FILTER="0.90"
        MISMATCH_QUALITY=27
        MISMATCH_QUALITY_TC=27
    fi
    printf " --------------------\n Read filter settings\n --------------------\n Min mapping quality (MAPQ): ${BOLD}${MAPQ_FILTER}${NORM}\n Max edit distance (NM): ${BOLD}${NM_FILTER}${NORM}\n Min fraction of matching bases (XF): ${BOLD}${MISMATCH_FILTER}${NORM}\n Min mismatch call quality (QV): ${BOLD}${MISMATCH_QUALITY}${NORM}\n Min mismatch call quality TC (QVTC): ${BOLD}${MISMATCH_QUALITY_TC}${NORM}\n Counting only unique reads: ${BOLD}${UNIQ_REP}${NORM}\n"

    # echo $INDEL_BAN_VAL
    if ! [[ -z $INDEL_BAN_ARG ]] ; then
        if [[ $INDEL_BAN_VAL =~ disabled ]] ; then
            printf " ----------------------------------\n"
        elif [[ $INDEL_BAN_VAL =~ ban ]] ; then
            printf " Ignoring reads with indels: ${BOLD}YES${NORM}\n ------------------------------------\n"
        else
            printf " Ignoring reads with indels ${BOLD}<= ${INDEL_BAN_VAL}${NORM} nucleotides.\n ------------------------------------\n"
        fi
    else
        printf " ----------------------------------\n"
    fi

    XI_VALUE=$(echo ${MISMATCH_FILTER}*100|bc)
    REP1="Min mapping quality (MAPQ):\t${MAPQ_FILTER}"
    REP2="Max edit distance (NM):\t${NM_FILTER}"
    REP3="Min fraction of matching bases (XF):\t${XI_VALUE}%%"
    REP4="Min mismatch call quality (QV):\t${MISMATCH_QUALITY}"
    REP5="Min mismatch call quality TC (QVTC):\t${MISMATCH_QUALITY_TC}"
    REP6="Counting only unique reads:\t${UNIQ_REP}"
    REP7="Indel ban threshold:\t${INDEL_BAN_VAL}"
fi

if ! [[ $THR_ARG =~ ^[0-9]+$ ]] ; then
    echo "Threads number must be an integer!"
    exit_err
fi

##### main routine #####

START_TIME=$SECONDS

if [[ $MODE_ARG == "htseq_only" ]]; then # test mode
    for BAM in $BAM_LIST
    do
        BAM_NAME=$( echo $BAM | xargs basename | awk -F. '{print $1}' )
        BAM_FILE=$( echo $BAM | xargs basename)

        echo "Running only HTSseq. Unused parameters will be ignored."

        if [[ $MD_ARG -eq 1 ]] ; then

            samtools calmd -@ $((THR_ARG-1)) $BAM $REF_ARG 2> /dev/null |\

            samtools collate -@$((THR_ARG-1)) -O --output-fmt SAM - 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_samtools_stderr.log

        else
            samtools collate -@$((THR_ARG-1)) -O --output-fmt SAM $BAM 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_samtools_stderr.log
        fi |\

        htseq-count -o ${OUTDIR_ARG}/${BAM_NAME}_hts_strand.bam -p BAM -r name -a 0 -f sam -s $STRANDEDNESS_ARG --nonunique $([[ $UNIQ_ARG -eq 1 ]] && echo "none" || echo "all") -i gene_id - $GTF_ARG > ${OUTDIR_ARG}/${BAM_NAME}_counts.txt

            echo "Finished!"
    done

    echo "All completed. Exiting."
    exit 0

fi

for BAM in $BAM_LIST
do

    BAM_NAME=$( echo $BAM | xargs basename | awk -F. '{print $1}' )
    BAM_FILE=$( echo $BAM | xargs basename )

#####################################
    printf "Processing%s ${BOLD}$BAM_FILE${NORM}\n"
#####################################

    if ! [[ $PROCENT == "1.00" ]] ; then
        echo "Subsampling "$(echo "scale=1;$PROCENT*100" | bc)"% of all reads"
        samtools view -h -b -s $PROCENT $BAM > ${OUTDIR_ARG}/${BAM_NAME}_subsampled.bam
        BAM=${OUTDIR_ARG}/${BAM_NAME}_subsampled.bam
        BAM_NAME=$( echo $BAM | xargs basename | awk -F. '{print $1}' )
        BAM_FILE=$( echo $BAM | xargs basename )
        echo "Done!"
        echo "--> Indexing"
        samtools index $BAM
        echo "Done!"
        echo "New (subsampled) BAM: $BAM"
    fi

    if [[ $MODE_ARG == "normal" ]] ; then
        echo "--> Counting reads."

        if [[ $MD_ARG -eq 1 ]] ; then

            samtools calmd -@ $((THR_ARG-1)) $BAM $REF_ARG 2> /dev/null |\
            samtools collate -@$((THR_ARG-1)) -O --output-fmt SAM - 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_samtools_stderr.log

        else

            samtools collate -@$((THR_ARG-1)) -O --output-fmt SAM $BAM 2>> $OUTDIR_ARG/stderr/${BAM_NAME}_samtools_stderr.log

        fi |\

        htseq-count -o ${OUTDIR_ARG}/${BAM_NAME}_hts_strand.bam -p BAM -r name -a 0 -f sam -s $STRANDEDNESS_ARG --nonunique $([[ $UNIQ_ARG -eq 1 ]] && echo "none" || echo "all") -i gene_id - $GTF_ARG > ${OUTDIR_ARG}/${BAM_NAME}_counts.txt

    elif [[ $MODE_ARG == "nocount" ]] ; then # another test mode

        echo "Skipping counting with HTseq."

    fi

    echo "Processing reads"

    TOTAL_READS=$(samtools view -c -@$((THR_ARG-1)) ${OUTDIR_ARG}/${BAM_NAME}_hts_strand.bam)
    samtools view -@$((THR_ARG-1)) ${OUTDIR_ARG}/${BAM_NAME}_hts_strand.bam |\
    gawk -v TOTAL_READS=$TOTAL_READS -v MAPPER_ARG=$MAPPER_ARG -v MODE_ARG=$MODE_ARG -v MAPQ_FILTER=$MAPQ_FILTER -v NM_FILTER=$NM_FILTER -v MISMATCH_FILTER=$MISMATCH_FILTER -v MISMATCH_QUALITY_TC=$MISMATCH_QUALITY_TC -v MISMATCH_QUALITY=$MISMATCH_QUALITY -v INDEL_BAN_ARG=$INDEL_BAN_ARG -v REBOUND_COMMAND="$REBOUND_COMMAND" -v BAM_NAME="$BAM_NAME" -v OUTDIR_ARG="$OUTDIR_ARG" -v STRANDEDNESS_ARG="$STRANDEDNESS_ARG" -f $AWK_FILE 2> $OUTDIR_ARG/stderr/${BAM_NAME}_rebound_stderr.log

    REBOUND_EXIT=$?

    if [[ $REBOUND_EXIT -eq 0 ]]; then
        TOTAL_TIME=$(($SECONDS - $START_TIME))
        echo "Successfully finished processing" $(echo "$BAM_LIST" | wc -w | awk '{print $1}') "file(s) in" $(echo "scale=2;$TOTAL_TIME/60" | bc) "minute(s)."
        echo "Mismatch counts saved to ${OUTDIR_ARG}/${BAM_NAME}_slam_counts.txt"
        cat $OUTDIR_ARG/${BAM_NAME}_rebound_stats.log
        printf "${REP1}\n${REP2}\n${REP3}\n${REP4}\n${REP5}\n${REP6}\n${REP7}\n" >> ${OUTDIR_ARG}/${BAM_NAME}_rebound_stats.log
    else
        echo "${BOLD}Processing $BAM_FILE FAILED${NORM}"
    fi
done
