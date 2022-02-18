# exclude fragments mapping to different genes
# for positions, look at same read names (groups) for different mismatches at the same position; this is indicative of overlaps; use dplyr (count groups > 1 (2)) to make a list of such reads; this may possible be faster way to filter for overlaps; at first, just collapse (maybe with uniq shell function)
# add functionality for unstranded experiments
# enable strandedness for artifacts detection (strand-soecific artifacts in the data)
# make useful for all mismatches equally
BEGIN{
    start_time = systime()
    OFS="\t"
    FS="\t"
    SUBSEP="@"
    for (n=0;n<256;n++){
        phred_conv[sprintf("%c",n+33)]=n
    }
    to = 0 # total overlap
    pair = "TRUE"
    name_read_previous = "NULL" # starting read for mate comparison (no read)
    mm_read_previous = "NULL" # starting read for mate comparison (no read)
    len_MDN_previous = 0

    print "S(L)AM converter started @ " strftime() > "/dev/stderr"

    if (MODE_ARG == "normal"){
        print "gene_id", "read_name", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", "OVER", "DUP" > OUTDIR_ARG "/" BAM_NAME "_RMLizer_counts.txt"
    }
}

# --- FUNCTIONS --- #

function tag_finder(tag,    i){
    if (substr($0,1,1)!="@"){
        for (i=1; i<=NF+1; i++){
            if (i==NF+1) {
                print "Fatal error: unable to find " tag " at record " $0 > "/dev/stderr"
                exit 1
            }else if (substr($i,1,5)!=tag){
                continue
            }else{
                return $i
            }
        }
    }
}

function copy_array(orig, copy,     i){
    delete copy
    for (i in orig){
        copy[i] = orig[i]
    }
}

# returns matched + mismatched (M) portion of the reads and the corresponding quality string

function MM_extractor(read_, quali_, cig_, md_value_,     md_ext, cig_ext, r, q, i, j, k){
    delete cig_ext; delete md_ext; delete mmdel_md_ext; delete mmdel_quali; delete mmdel_read; delete q; delete r

    ### expand MD tag ###
    n = patsplit(md_value_, b, /[0-9]+/, seps)
    k = 1
    for (i=1; i<=n; i++){
        for(j=1; j<=b[i]; j++){
            md_ext[k] = "."
            k++
        }
        if (seps[i] ~ /[ACTGN]/){
            md_ext[k] = seps[i]
            k++
        }
    }

    ### expand CIGAR ###
    n = patsplit(cig_, b, /[^[:digit:]]/, seps)
    # loop through CIGAR "segments" to get their lengths (--> seps)
    k = 1
    for (i=0; i<=n-1; i++){
        # iteratively add letters CIGAR operations
        for (j=1; j<=seps[i]; j++){
            cig_ext[k] = b[i+1]
            k++
        }
    }

    split(quali_, q, "")
    split(read_, r, "")

    len_MDN = 0 # cumulative lenght of Ms, Ds and Ns
    j = 1
    k = 1

    for (i in cig_ext) {
        if (cig_ext[i] ~ /[IS]/){
            j++
        }else if (cig_ext[i] ~ /[DN]/){
            len_MDN++
        }else if (cig_ext[i] ~ /[M=X]/){
            len_MDN++
            mmdel_read[len_MDN] = r[j]
            mmdel_quali[len_MDN] = q[j]
            mmdel_md_ext[len_MDN] = md_ext[k]
            j++
            k++
        }else{
            print "Fatal error: Unexpected character(s) in CIGAR string of read" $0 > "/dev/stderr"
            exit 1
        }
    }
}

function MM_extractor_simple(read_, cig_,     cig_ext,r,i,j,k){
    delete cig_ext; delete mmdel_md_ext; delete mmdel_read; delete r

    ### expand CIGAR ###
    n = patsplit(cig_, b, /[^[:digit:]]/, seps)
    # loop through CIGAR "segments" to get their lengths (--> seps)
    k = 1
    for (i=0; i<=n-1; i++){
        # iteratively add letters CIGAR operations
        for (j=1; j<=seps[i]; j++){
            cig_ext[k]=b[i+1]
            k++
        }
    }

    split(read_, r, "")

    len_MDN=0
    j=1
    k=1

    for (i in cig_ext){
        if (cig_ext[i] ~ /[IS]/){
            j++
        }else if (cig_ext[i] ~ /[DN]/){
            len_MDN++
        }else if (cig_ext[i] ~ /[M=X]/){
            len_MDN++
            mmdel_read[len_MDN] = r[j]
            mmdel_md_ext[len_MDN] = md_ext[k]
            j++
            k++
        }else{
            print "Fatal error: Unexpected character(s) in CIGAR string of read" > "/dev/stderr"
            print $0 > "/dev/stderr"
            exit 1
        }
    }
}

function TCRA_producer_for(name_read_current_, name_read_previous_, mmdel_read_current_, mmdel_read_previous_, mmdel_quali_current_, mmdel_md_ext_, len_MDN_current_, len_MDN_previous_,    i,a){
    delete mmms; delete mmms_overlap; delete count_all; delete count_overlap; delete cnt
    os=0

    if (name_read_current_ != name_read_previous_){ # if previous read is NOT the mate read of the current read ...
        fr_start = $4
        fr_end = $4 + len_MDN_current_
        for (i in mmdel_md_ext_){
            if (mmdel_md_ext_[i] == "."){
                mmms[i] = mmdel_read_current_[i] mmdel_read_current_[i]
            }else if (mmdel_md_ext_[i] == ""){
                print "Mismatched array indexes (error in function \"TCRA_producer_for\"); read skipped" > "/dev/stderr"
                print $0 > "/dev/stderr"
                # exit 1
                next
            }else{
                if (phred_conv[mmdel_quali_current_[i]] < MISMATCH_QUALITY){
                    mmms[i] = mmdel_md_ext_[i] mmdel_md_ext_[i]
                }else{
                    mmms[i] = mmdel_md_ext_[i] mmdel_read_current_[i]
                }
            }
        }
        for (i in mmms){
            cnt["mmms"]++
            ###### A -> X
            if       (mmms[i] == "AG"){
                count_all["AG"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++ ;
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "AA"){
                count_all["AA"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                cnt["match"]++
            }else if (mmms[i] == "AC"){
                count_all["AC"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "AT"){
                count_all["AT"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "AN"){
                count_all["AN"]++
            ###### C -> X
            }else if (mmms[i] == "CA"){
                count_all["CA"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "CC"){
                count_all["CC"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                cnt["match"]++
            }else if (mmms[i] == "CG"){
                count_all["CG"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "CT"){
                count_all["CT"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "CN"){
                count_all["CN"]++
            ###### G -> X
            }else if (mmms[i] == "GA"){
                count_all["GA"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "GC"){
                count_all["GC"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "GG"){
                count_all["GG"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                cnt["match"]++
            }else if (mmms[i] == "GT"){
                count_all["GT"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "GN"){
                count_all["GN"]++
            ###### T -> X
            }else if (mmms[i] == "TA"){
                count_all["TA"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "TC"){
                if (phred_conv[mmdel_quali_current_[i]] >= MISMATCH_QUALITY_TC){
                    count_all["TC"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    if (phred_conv[mmdel_quali_current_[i]] > 41 || phred_conv[mmdel_quali_current_[i]] < 0){
                        print "Error. Quality scores out of range. Is it Illumina 1.8+ (Phred 33+) encoding?" > "/dev/stderr"
                        print $0 > "/dev/stderr"
                        exit 1
                    }
                }
            }else if (mmms[i] == "TG"){
                count_all["TG"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "TT"){
                count_all["TT"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                cnt["match"]++
            }else if (mmms[i] == "TN"){
                count_all["TN"]++
            ###### N -> X
            }else if (mmms[i] == "NA"){
                count_all["NA"]++
            }else if (mmms[i] == "NC"){
                count_all["NC"]++
            }else if (mmms[i] == "NG"){
                count_all["NG"]++
            }else if (mmms[i] == "NT"){
                count_all["NT"]++
            }else if (mmms[i] == "NN"){
                count_all["NN"]++
            }
        }
    }else{ # if mate found (read current == read previous)
        pos_shift = $8-$4
        if (pos_shift <= 0){
            fr_start = $8
            fr_end = $4 + len_MDN_current_
        }else if (pos_shift > 0){
            fr_start = $4
            fr_end = $8 + len_MDN_previous_
        }
        for (i in mmdel_md_ext_){
            if (mmdel_md_ext_[i] == "."){ # match
                mmms[i] = mmdel_read_current_[i] mmdel_read_current_[i]
                if (pos_shift <= 0){ # checks if the mate starts in the same position or downstream
                    if (i+0 <= pos_shift + len_MDN_previous_){
                        mmms_overlap[i] = mmdel_read_current_[i] mmdel_read_current_[i] mmdel_read_previous_[i-pos_shift]
                    }
                }else if (pos_shift > 0){ # check if the mate starts upstream
                    if (i+0 > pos_shift && i+0 <= len_MDN_current_){
                        mmms_overlap[i] = mmdel_read_current_[i] mmdel_read_current_[i] mmdel_read_previous_[i-pos_shift]
                    }
                }else{
                    print "Mate offset incorrect ("pos_shift"). Perhaps the data is single-ended? Record:" > "/dev/stderr"
                    print $0 > "/dev/stderr"
                    exit 1
                }
            }else if (mmdel_md_ext_[i] == ""){ # (D)eletion
                print "Mismatched array indexes (error in function \"TCRA_producer_for\"); read skipped" > "/dev/stderr"
                print $0 > "/dev/stderr"
                # exit 1
                next
            }else{ # mismatch
                if (phred_conv[mmdel_quali_current_[i]] < MISMATCH_QUALITY){
                    mmms[i] = mmdel_md_ext_[i] mmdel_md_ext_[i]
                    if (pos_shift <= 0){
                        if (i+0 <= pos_shift + len_MDN_previous_){
                            mmms_overlap[i] = mmdel_md_ext_[i] mmdel_md_ext_[i] mmdel_read_previous_[i-pos_shift]
                        }
                    }else if (pos_shift > 0){
                        if (i+0 > pos_shift && i+0 <= len_MDN_current_){
                            mmms_overlap[i] = mmdel_md_ext_[i] mmdel_md_ext_[i] mmdel_read_previous_[i-pos_shift]
                        }
                    }
                }else{ # quality filter TC passed
                    mmms[i] = mmdel_md_ext_[i] mmdel_read_current_[i]
                    if (pos_shift <= 0){
                        if (i+0 <= pos_shift + len_MDN_previous_){
                            mmms_overlap[i] = mmdel_md_ext_[i] mmdel_read_current_[i] mmdel_read_previous_[i-pos_shift]
                        }
                    }else if (pos_shift > 0){
                        if (i+0 > pos_shift && i+0 <= len_MDN_current_){
                            mmms_overlap[i] = mmdel_md_ext_[i] mmdel_read_current_[i] mmdel_read_previous_[i-pos_shift]
                        }
                    }else{
                        print "Mate offset incorrect ("pos_shift"). Perhaps the data is single-ended? Record:" > "/dev/stderr"
                        print $0 > "/dev/stderr"
                        exit 1
                    }
                }
            }
        }
        for (i in mmms){
            cnt["mmms"]++
            if (mmms_overlap[i] == ""){
                ###### A -> X
                if       (mmms[i] == "AG"){
                    count_all["AG"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "AA"){
                    count_all["AA"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms[i] == "AC"){
                    count_all["AC"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "AT"){
                    count_all["AT"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "AN"){
                    count_all["AN"]++
                ###### C -> X
                }else if (mmms[i] == "CA"){
                    count_all["CA"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "CC"){
                    count_all["CC"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms[i] == "CG"){
                    count_all["CG"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "CT"){
                    count_all["CT"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "CN"){
                    count_all["CN"]++
                ###### G -> X
                }else if (mmms[i] == "GA"){
                    count_all["GA"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "GC"){
                    count_all["GC"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "GG"){
                    count_all["GG"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms[i] == "GT"){
                    count_all["GT"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "GN"){
                    count_all["GN"]++
                ###### T -> X
                }else if (mmms[i] == "TA"){
                    count_all["TA"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "TC"){
                    if (phred_conv[mmdel_quali_current_[i]] >= MISMATCH_QUALITY_TC){
                        count_all["TC"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                        print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        if (phred_conv[mmdel_quali_current_[i]] > 41 || phred_conv[mmdel_quali_current_[i]] < 0){
                            print "Error. Quality scores out of range. Is it Illumina 1.8+ (Phred 33+) encoding?" > "/dev/stderr"
                            print $0 > "/dev/stderr"
                            exit 1
                        }
                    }
                }else if (mmms[i] == "TG"){
                    count_all["TG"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "TT"){
                    count_all["TT"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms[i] == "TN"){
                    count_all["TN"]++
                ###### N -> X
                }else if (mmms[i] == "NA"){
                    count_all["NA"]++
                }else if (mmms[i] == "NC"){
                    count_all["NC"]++
                }else if (mmms[i] == "NG"){
                    count_all["NG"]++
                }else if (mmms[i] == "NT"){
                    count_all["NT"]++
                }else if (mmms[i] == "NN"){
                    count_all["NN"]++
                }
            }else{
                os++ # overlap size counter
                to++ # total overlap size counter
                ###### A -> X
                if       (mmms_overlap[i] == "AGG"){
                    count_overlap["AG"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "AAA"){
                    count_overlap["AA"]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms_overlap[i] == "ACC"){
                    count_overlap["AC"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "ATT"){
                    count_overlap["AT"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "ANN"){
                ###### C -> X
                }else if (mmms_overlap[i] == "CAA"){
                    count_overlap["CA"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "CCC"){
                    count_overlap["CC"]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms_overlap[i] == "CGG"){
                    count_overlap["CG"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "CTT"){
                    count_overlap["CT"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "CNN"){
                    count_overlap["CN"]++
                ###### G -> X
                }else if (mmms_overlap[i] == "GAA"){
                    count_overlap["GA"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "GCC"){
                    count_overlap["GC"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "GGG"){
                    count_overlap["GG"]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms_overlap[i] == "GTT"){
                    count_overlap["GT"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "GNN"){
                    count_overlap["GN"]++
                ###### T -> X
                }else if (mmms_overlap[i] == "TAA"){
                    count_overlap["TA"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "TCC"){
                    if (phred_conv[mmdel_quali_current_[i]] >= MISMATCH_QUALITY_TC){
                        count_overlap["TC"]++
                        print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        if (phred_conv[mmdel_quali_current_[i]] > 41 || phred_conv[mmdel_quali_current_[i]] < 0){
                            print "Error. Quality scores out of range. Is it Illumina 1.8+ (Phred 33+) encoding?" > "/dev/stderr"
                            print $0 > "/dev/stderr"
                            exit 1
                        }
                    }
                }else if (mmms_overlap[i] == "TGG"){
                    count_overlap["TG"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "TTT"){
                    count_overlap["TT"]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms_overlap[i] == "TNN"){
                    count_overlap["TN"]++
                ###### N -> X
                }else if (mmms_overlap[i] == "NAA"){
                    count_overlap["NA"]++
                }else if (mmms_overlap[i] == "NCC"){
                    count_overlap["NC"]++
                }else if (mmms_overlap[i] == "NGG"){
                    count_overlap["NG"]++
                }else if (mmms_overlap[i] == "NTT"){
                    count_overlap["NT"]++
                }else if (mmms_overlap[i] == "NNN"){
                    count_overlap["NN"]++
                #############################################
                }else if (mmms_overlap[i] == "A[^G]G"){
                    count_all["AG"]--
                }else if (mmms_overlap[i] == "A[^A]A"){
                    count_all["AA"]--
                }else if (mmms_overlap[i] == "A[^C]C"){
                    count_all["AC"]--
                }else if (mmms_overlap[i] == "A[^T]T"){
                    count_all["AT"]--
                }else if (mmms_overlap[i] == "A[^N]N"){
                ###### C -> X
                }else if (mmms_overlap[i] == "C[^A]A"){
                    count_all["CA"]--
                }else if (mmms_overlap[i] == "C[^C]C"){
                    count_all["CC"]--
                }else if (mmms_overlap[i] == "C[^G]G"){
                    count_all["CG"]--
                }else if (mmms_overlap[i] == "C[^T]T"){
                    count_all["CT"]--
                }else if (mmms_overlap[i] == "C[^N]N"){
                    count_all["CN"]--
                ###### G -> X
                }else if (mmms_overlap[i] == "G[^A]A"){
                    count_all["GA"]--
                }else if (mmms_overlap[i] == "G[^C]C"){
                    count_all["GC"]--
                }else if (mmms_overlap[i] == "G[^G]G"){
                    count_all["GG"]--
                }else if (mmms_overlap[i] == "G[^T]T"){
                    count_all["GT"]--
                }else if (mmms_overlap[i] == "G[^N]N"){
                    count_all["GN"]--
                ###### T -> X
                }else if (mmms_overlap[i] == "T[^A]A"){
                    count_all["TA"]--
                }else if (mmms_overlap[i] == "T[^C]C"){
                    count_all["TC"]--
                }else if (mmms_overlap[i] == "T[^G]G"){
                    count_all["TG"]--
                }else if (mmms_overlap[i] == "T[^T]T"){
                    count_all["TT"]--
                }else if (mmms_overlap[i] == "T[^N]N"){
                    count_all["TN"]--
                ###### N -> X
                }else if (mmms_overlap[i] == "N[^A]A"){
                    count_all["NA"]--
                }else if (mmms_overlap[i] == "N[^C]C"){
                    count_all["NC"]--
                }else if (mmms_overlap[i] == "N[^G]G"){
                    count_all["NG"]--
                }else if (mmms_overlap[i] == "N[^T]T"){
                    count_all["NT"]--
                }else if (mmms_overlap[i] == "N[^N]N"){
                    count_all["NN"]--
                }
            }
        }
    }

    sum_count_all = cnt["mmms"]+0
    if (sum_count_all > 0){
        xi = (cnt["match"]+0) / sum_count_all
    }else{ # reads with zero counts will be excluded anyways
        xi = 0
    }
}

function TCRA_producer_rev(name_read_current_, name_read_previous_, mmdel_read_current_, mmdel_read_previous_, mmdel_quali_current_, mmdel_md_ext_, len_MDN_current_, len_MDN_previous_,      i,a){
    delete mmms; delete mmms_overlap; delete count_all; delete count_overlap; delete cnt
    os=0
    if (name_read_current_ != name_read_previous_){ # if previous read is NOT the mate read of the current
        fr_start = $4
        fr_end = $4 + len_MDN_current_
        for (i in mmdel_md_ext_){
            if (mmdel_md_ext_[i] == "."){
                mmms[i] = mmdel_read_current_[i] mmdel_read_current_[i]
            }else if (mmdel_md_ext_[i] == ""){
                print "Mismatched array indexes (error in function \"TCRA_producer_rev\"); read skipped" > "/dev/stderr"
                print $0 > "/dev/stderr"
                # exit 1
                next
            }else{
                if (phred_conv[mmdel_quali_current_[i]] < MISMATCH_QUALITY){
                    mmms[i] = mmdel_md_ext_[i] mmdel_md_ext_[i]
                }else{
                    mmms[i] = mmdel_md_ext_[i] mmdel_read_current_[i]
                }
            }
        }
        for (i in mmms){ # REVERSE
            cnt["mmms"]++
            ###### A -> X
            # print phred_conv[mmdel_quali_current_[i]],"---", i
            if (mmms[i] == "AG"){
                if (phred_conv[mmdel_quali_current_[i]] >= MISMATCH_QUALITY_TC){
                    count_all["TC"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    if (phred_conv[mmdel_quali_current_[i]] > 41 || phred_conv[mmdel_quali_current_[i]] < 0){
                        print "Error. Quality scores out of range. Is it Illumina 1.8+ (Phred 33+) encoding?" > "/dev/stderr"
                        print $0 > "/dev/stderr"
                        exit 1
                    }
                }
            }else if (mmms[i] == "AA"){
                count_all["TT"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                cnt["match"]++
            }else if (mmms[i] == "AC"){
                count_all["TG"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "AT"){
                count_all["TA"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "AN"){
                count_all["TN"]++
            ###### C -> X
            }else if (mmms[i] == "CA"){
                count_all["GT"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "CC"){
                count_all["GG"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                cnt["match"]++
            }else if (mmms[i] == "CG"){
                count_all["GC"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "CT"){
                count_all["GA"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "CN"){
                count_all["GN"]++
            ###### G -> X
            }else if (mmms[i] == "GA"){
                count_all["CT"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "GC"){
                count_all["CG"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "GG"){
                count_all["CC"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                cnt["match"]++
            }else if (mmms[i] == "GT"){
                count_all["CA"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "GN"){
                count_all["CN"]++
            ###### T -> X
            }else if (mmms[i] == "TA"){
                count_all["AT"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "TC"){
                count_all["AG"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "TG"){
                count_all["AC"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmms[i] == "TT"){
                count_all["AA"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                cnt["match"]++
            }else if (mmms[i] == "TN"){
                count_all["AN"]++
            ###### N -> X
            }else if (mmms[i] == "NA"){
                count_all["NT"]++
            }else if (mmms[i] == "NC"){
                count_all["NG"]++
            }else if (mmms[i] == "NG"){
                count_all["NC"]++
            }else if (mmms[i] == "NT"){
                count_all["NA"]++
            }else if (mmms[i] == "NN"){
                count_all["NN"]++
            }
        }
    }else{ # if mate found
        pos_shift = $8-$4
        if (pos_shift <= 0){
            fr_start = $8
            fr_end = $4 + len_MDN_current_
        }else if (pos_shift > 0){
            fr_start = $4
            fr_end = $8 + len_MDN_previous_
        }
        for (i in mmdel_md_ext_){
            if (mmdel_md_ext_[i] == "."){ # match
                mmms[i] = mmdel_read_current_[i] mmdel_read_current_[i]
                if (pos_shift <= 0){ # check if the mate starts in the same position or downstream
                    if (i+0 <= pos_shift + length(mmdel_read_previous_)){
                        mmms_overlap[i] = mmdel_read_current_[i] mmdel_read_current_[i] mmdel_read_previous_[i-pos_shift]
                    }
                }else if (pos_shift > 0){ # check if the mate starts upstream
                    if (i+0 > pos_shift && i+0 <= length(mmdel_read_current_)){
                        mmms_overlap[i] = mmdel_read_current_[i] mmdel_read_current_[i] mmdel_read_previous_[i-pos_shift]
                    }
                }else{
                    print "Mate offset incorrect ("pos_shift"). Perhaps the data is single-ended? Record:" > "/dev/stderr"
                    print $0 > "/dev/stderr"
                    exit 1
                }
            }else if (mmdel_md_ext_[i] == ""){ # (D)eletion
                print "Mismatched array indexes (error in function \"TCRA_producer_rev\"); read skipped" > "/dev/stderr"
                print $0 > "/dev/stderr"
                # exit 1
                next
            }else{ # mismatch
                if (phred_conv[mmdel_quali_current_[i]] < MISMATCH_QUALITY){
                    mmms[i] = mmdel_md_ext_[i] mmdel_md_ext_[i]
                    if (pos_shift <= 0){
                        if (i+0 <= pos_shift + len_MDN_previous_){
                            mmms_overlap[i] = mmdel_md_ext_[i] mmdel_md_ext_[i] mmdel_read_previous_[i-pos_shift]
                        }
                    }else if (pos_shift > 0){
                        if (i+0 > pos_shift && i+0 <= len_MDN_current_){
                            mmms_overlap[i] = mmdel_md_ext_[i] mmdel_md_ext_[i] mmdel_read_previous_[i-pos_shift]
                        }
                    }
                }else{
                    mmms[i] = mmdel_md_ext_[i] mmdel_read_current_[i]
                    if (pos_shift <= 0){
                        if (i+0 <= pos_shift + len_MDN_previous_){
                            mmms_overlap[i] = mmdel_md_ext_[i] mmdel_read_current_[i] mmdel_read_previous_[i-pos_shift]
                        }
                    }else if (pos_shift > 0){
                        if (i+0 > pos_shift && i+0 <= len_MDN_current_){
                            mmms_overlap[i] = mmdel_md_ext_[i] mmdel_read_current_[i] mmdel_read_previous_[i-pos_shift]
                        }
                    }else{
                        print "Mate offset incorrect ("pos_shift"). Perhaps the data is single-ended? Record:" > "/dev/stderr"
                        print $0 > "/dev/stderr"
                        exit 1
                    }
                }
            }
        }
        for (i in mmms){ # REVERSE
            cnt["mmms"]++
            # print phred_conv[mmdel_quali_current_[i]],"---", i, mmms[i], mmms_overlap[i]
            if (mmms_overlap[i] == ""){
                ###### A -> X
                if (mmms[i] == "AG"){
                    if (phred_conv[mmdel_quali_current_[i]] >= MISMATCH_QUALITY_TC){
                        count_all["TC"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                        print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        if (phred_conv[mmdel_quali_current_[i]] > 41 || phred_conv[mmdel_quali_current_[i]] < 0){
                            print "Error. Quality scores out of range. Is it Illumina 1.8+ (Phred 33+) encoding?" > "/dev/stderr"
                            print $0 > "/dev/stderr"
                            exit 1
                        }
                    }
                }else if (mmms[i] == "AA"){
                    count_all["TT"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms[i] == "AC"){
                    count_all["TG"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "AT"){
                    count_all["TA"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "AN"){
                    count_all["TN"]++
                ###### C -> X
                }else if (mmms[i] == "CA"){
                    count_all["GT"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "CC"){
                    count_all["GG"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms[i] == "CG"){
                    count_all["GC"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "CT"){
                    count_all["GA"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "CN"){
                    count_all["GN"]++
                ###### G -> X
                }else if (mmms[i] == "GA"){
                    count_all["CT"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "GC"){
                    count_all["CG"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "GG"){
                    count_all["CC"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms[i] == "GT"){
                    count_all["CA"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "GN"){
                    count_all["CN"]++
                ###### T -> X
                }else if (mmms[i] == "TA"){
                    count_all["AT"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "TC"){
                    count_all["AG"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "TG"){
                    count_all["AC"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms[i] == "TT"){
                    count_all["AA"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms[i] == "TN"){
                    count_all["AN"]++
                ###### N -> X
                }else if (mmms[i] == "NA"){
                    count_all["NT"]++
                }else if (mmms[i] == "NC"){
                    count_all["NG"]++
                }else if (mmms[i] == "NG"){
                    count_all["NC"]++
                }else if (mmms[i] == "NT"){
                    count_all["NA"]++
                }else if (mmms[i] == "NN"){
                    count_all["NN"]++
                }
            }else{
                os++ # overlap size counter
                to++ # total overlap size counter
                ###### A -> X
                if (mmms_overlap[i] == "AGG"){ # REVERSE!!!!!
                    if (phred_conv[mmdel_quali_current_[i]] >= MISMATCH_QUALITY_TC){
                        count_overlap["TC"]++
                        print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        if (phred_conv[mmdel_quali_current_[i]] > 41 || phred_conv[mmdel_quali_current_[i]] < 0){
                            print "Error. Quality scores out of range. Is it Illumina 1.8+ (Phred 33+) encoding?" > "/dev/stderr"
                            print $0 > "/dev/stderr"
                            exit 1
                        }
                    }
                }else if (mmms_overlap[i] == "AAA"){
                    count_overlap["TT"]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms_overlap[i] == "ACC"){
                    count_overlap["TG"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "ATT"){
                    count_overlap["TA"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "ANN"){
                    count_overlap["TN"]++
                ###### C -> X
                }else if (mmms_overlap[i] == "CAA"){
                    count_overlap["GT"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "CCC"){
                    count_overlap["GG"]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms_overlap[i] == "CGG"){
                    count_overlap["GC"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "CTT"){
                    count_overlap["GA"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "CNN"){
                    count_overlap["GN"]++
                ###### G -> X
                }else if (mmms_overlap[i] == "GAA"){
                    count_overlap["CT"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "GCC"){
                    count_overlap["CG"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "GGG"){
                    count_overlap["CC"]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms_overlap[i] == "GTT"){
                    count_overlap["CA"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "GNN"){
                    count_overlap["CN"]++
                ###### T -> X
                }else if (mmms_overlap[i] == "TAA"){
                    count_overlap["AT"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "TCC"){
                    count_overlap["AG"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "TGG"){
                    count_overlap["AC"]++
                    print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                }else if (mmms_overlap[i] == "TTT"){
                    count_overlap["AA"]++
                    # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    cnt["match"]++
                }else if (mmms_overlap[i] == "TNN"){
                    count_overlap["AN"]++
                ###### N -> X
                }else if (mmms_overlap[i] == "NAA"){
                    count_overlap["NT"]++
                }else if (mmms_overlap[i] == "NCC"){
                    count_overlap["NG"]++
                }else if (mmms_overlap[i] == "NGG"){
                    count_overlap["NC"]++
                }else if (mmms_overlap[i] == "NTT"){
                    count_overlap["NA"]++
                }else if (mmms_overlap[i] == "NNN"){
                    count_overlap["NN"]++
                ####################################
                }else if (mmms_overlap[i] == "A[^G]G"){ # REVERSE!!!!!
                    count_all["TC"]--
                }else if (mmms_overlap[i] == "A[^A]A"){
                    count_all["TT"]--
                }else if (mmms_overlap[i] == "A[^C]C"){
                    count_all["TG"]--
                }else if (mmms_overlap[i] == "A[^T]T"){
                    count_all["TA"]--
                }else if (mmms_overlap[i] == "A[^N]N"){
                    count_all["TN"]--
                ###### C -> X
                }else if (mmms_overlap[i] == "C[^A]A"){
                    count_all["GT"]--
                }else if (mmms_overlap[i] == "C[^C]C"){
                    count_all["GG"]--
                }else if (mmms_overlap[i] == "C[^G]G"){
                    count_all["GC"]--
                }else if (mmms_overlap[i] == "C[^T]T"){
                    count_all["GA"]--
                }else if (mmms_overlap[i] == "C[^N]N"){
                    count_all["GN"]--
                ###### G -> X
                }else if (mmms_overlap[i] == "G[^A]A"){
                    count_all["CT"]--
                }else if (mmms_overlap[i] == "G[^C]C"){
                    count_all["CG"]--
                }else if (mmms_overlap[i] == "G[^G]G"){
                    count_all["CC"]--
                }else if (mmms_overlap[i] == "G[^T]T"){
                    count_all["CA"]--
                }else if (mmms_overlap[i] == "G[^N]N"){
                    count_all["CN"]--
                ###### T -> X
                }else if (mmms_overlap[i] == "T[^A]A"){
                    count_all["AT"]--
                }else if (mmms_overlap[i] == "T[^C]C"){
                    count_all["AG"]--
                }else if (mmms_overlap[i] == "T[^G]G"){
                    count_all["AC"]--
                }else if (mmms_overlap[i] == "T[^T]T"){
                    count_all["AA"]--
                }else if (mmms_overlap[i] == "T[^N]N"){
                    count_all["AN"]--
                ###### N -> X
                }else if (mmms_overlap[i] == "N[^A]A"){
                    count_all["NT"]--
                }else if (mmms_overlap[i] == "N[^C]C"){
                    count_all["NG"]--
                }else if (mmms_overlap[i] == "N[^G]G"){
                    count_all["NC"]--
                }else if (mmms_overlap[i] == "N[^T]T"){
                    count_all["NA"]--
                }else if (mmms_overlap[i] == "N[^N]N"){
                    count_all["NN"]--
                }
            }
        }
    }

    sum_count_all = cnt["mmms"]+0
    if (sum_count_all > 0){
        xi = (cnt["match"]+0) / sum_count_all
    }else{
        xi = 0
    }
}

function TCRA_producer_simple_for(name_read_current_, name_read_previous_, mmdel_read_current_, mmdel_read_previous_, len_MDN_current_, len_MDN_previous_,      i,j){
    delete count_overlap; delete count_all
    os = 0 # overlap size counter
    if (name_read_current_ != name_read_previous_){ # previous read is NOT the mate read of the current read
        fr_start = $4
        fr_end = $4 + len_MDN_current_
        for (i in mmdel_read_current_){
            if (mmdel_read_current_[i] == "A"){
                count_all["AA"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmdel_read_current_[i] == "C"){
                count_all["CC"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmdel_read_current_[i] == "G"){
                count_all["GG"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmdel_read_current_[i] == "T"){
                count_all["TT"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmdel_read_current_[i] == "N"){
                count_all["NN"]++
            }
        }
    }else{ # previous read is the mate read of the current read
        pos_shift = $8-$4
        if (pos_shift <= 0){
            fr_start = $8
            fr_end = $4 + len_MDN_current_
            for (i in mmdel_read_current_){
                if (i+0 > pos_shift + len_MDN_previous_){
                    if (mmdel_read_current_[i] == "A"){
                        count_all["AA"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "C"){
                        count_all["CC"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "G"){
                        count_all["GG"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "T"){
                        count_all["TT"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "N"){
                        count_all["NN"]++
                    }
                }
                if (i+0 <= len_MDN_previous_ + pos_shift){
                    if (mmdel_read_current_[i] == "A"){
                        count_overlap["AA"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "C"){
                        count_overlap["CC"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "G"){
                        count_overlap["GG"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "T"){
                        count_overlap["TT"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "N"){
                        count_overlap["NN"]++
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }
                }
            }
        }else if (pos_shift > 0){ # check if the mate starts in the same position or downstream
            fr_start = $4
            fr_end = $8 + len_MDN_previous_
            for (i in mmdel_read_current_){
                if (i+0 <= pos_shift){
                    if (mmdel_read_current_[i] == "A"){
                        count_all["AA"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "C"){
                        count_all["CC"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "G"){
                        count_all["GG"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "T"){
                        count_all["TT"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "N"){
                        count_all["NN"]++
                    }
                }
                if (i+0 > pos_shift && i+0 <= len_MDN_current_){
                    if (mmdel_read_current_[i] == "A"){
                        count_overlap["AA"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "C"){
                        count_overlap["CC"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "G"){
                        count_overlap["GG"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "T"){
                        count_overlap["TT"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "N"){
                        count_overlap["NN"]++
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }
                }
            }
        }else{
            print "Mate offset incorrect ("pos_shift"). Perhaps the data is single-ended? Record:" > "/dev/stderr"
            print $0 > "/dev/stderr"
            exit 1
        }
    }
}

function TCRA_producer_simple_rev(name_read_current_, name_read_previous_, mmdel_read_current_, mmdel_read_previous_, len_MDN_current_, len_MDN_previous_,      i,j){
    delete count_overlap; delete count_all

    os = 0 # overlap size counter reset

    if (name_read_current_ != name_read_previous_){ # previous read is NOT the mate read of the current read
        # delete count_previous
        fr_start = $4
        fr_end = $4 + len_MDN_current_
        for (i in mmdel_read_current_){ # REVERSE !!!
            if (mmdel_read_current_[i] == "A"){
                count_all["TT"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmdel_read_current_[i] == "C"){
                count_all["GG"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmdel_read_current_[i] == "G"){
                count_all["CC"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmdel_read_current_[i] == "T"){
                count_all["AA"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
            }else if (mmdel_read_current_[i] == "N"){
                count_all["NN"]++
            }
        }
    }else{ # previous read is the mate read of the current read
        pos_shift = $8-$4
        if (pos_shift <= 0){
            fr_start = $8
            fr_end = $4 + len_MDN_current_
            for (i in mmdel_read_current_){
                if (i+0 > pos_shift + len_MDN_previous_){ # REVERSE !!!
                    if (mmdel_read_current_[i] == "A"){
                        count_all["TT"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "C"){
                        count_all["GG"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "G"){
                        count_all["CC"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "T"){
                        count_all["AA"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "N"){
                        count_all["NN"]++
                    }
                }
                if (i+0 <= len_MDN_previous_ + pos_shift){ # REVERSE !!!
                    if (mmdel_read_current_[i] == "A"){
                        count_overlap["TT"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "C"){
                        count_overlap["GG"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "G"){
                        count_overlap["CC"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "T"){
                        count_overlap["AA"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "N"){
                        count_overlap["NN"]++
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }
                }
            }
        }else if (pos_shift > 0){ # check if the mate starts in the same position or downstream
            fr_start = $4
            fr_end = $8 + len_MDN_previous_
            for (i in mmdel_read_current_){ # REVERSE !!!
                if (i+0 <= pos_shift){ # REVERSE !!!
                    if (mmdel_read_current_[i] == "A"){
                        count_all["TT"]++ ; depth[$3, $4+i-1, "T", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "T", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "C"){
                        count_all["GG"]++ ; depth[$3, $4+i-1, "G", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "G", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "G"){
                        count_all["CC"]++ ; depth[$3, $4+i-1, "C", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "C", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "T"){
                        count_all["AA"]++ ; depth[$3, $4+i-1, "A", gene_id, strand3, is_duplicate_current]++ ; breadth[$3, "A", gene_id, strand3, is_duplicate_current]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                    }else if (mmdel_read_current_[i] == "N"){
                        count_all["NN"]++
                    }
                }
                if (i+0 > pos_shift && i+0 <= len_MDN_current_){ # REVERSE !!!
                    if (mmdel_read_current_[i] == "A"){
                        count_overlap["TT"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "TT", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "C"){
                        count_overlap["GG"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "GG", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "G"){
                        count_overlap["CC"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "CC", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "T"){
                        count_overlap["AA"]++
                        # print $3, $4+i-2, $4+i-1, gene_id, strand3, name_read_current, "AA", is_duplicate_current > OUTDIR_ARG "/" BAM_NAME "_mismatch_positions_" $3 ".txt"
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }else if (mmdel_read_current_[i] == "N"){
                        count_overlap["NN"]++
                        os++ # overlap size counter
                        to++ # total overlap size counter
                    }
                }
            }
        }else{
            print "Mate offset incorrect ("pos_shift"). Perhaps the data is single-ended? Record:" > "/dev/stderr"
            print $0 > "/dev/stderr"
            exit 1
        }
    }
}

# ---- main routine ---- #

{
    # chunk size might be provided as a parameter
    # consider moving in case many record skipped, thoug 1000000 record unlikely to take less than 1s
    if ((NR/1000000)%1==0){

        timestamp=systime()

        for (i in depth){
            print i, depth[i] > OUTDIR_ARG "/" BAM_NAME "_depth_per_position_" timestamp ".txt"
        }
        delete depth

        for (i in breadth){
            print i, breadth[i] > OUTDIR_ARG "/" BAM_NAME "_breadth_per_gene_" timestamp ".txt"
        }
        delete breadth

    }

    if ($1 ~ /^@/){ # write header lines; replace this by adding the header separately via concatenating the SAM file instead
        header=1
        print $0 > OUTDIR_ARG "/" BAM_NAME "_RMLizer.sam"
        header_lines=NR
        next
    }

    if (NR < 100000){
        if ((NR/(10000-header_lines))%1==0){
            printf "--> %s out of %s records (%.2f\%) processed\n", NR, TOTAL_READS, (NR/TOTAL_READS)*100
        }
    }else{
        if ((NR/(100000-header_lines))%1==0){
            printf "--> %s out of %s records (%.2f\%) processed\n", NR, TOTAL_READS, (NR/TOTAL_READS)*100
        }
    }

    ############################## Save read name #########################################

    name_read_current = $1
    # split($1,a,":")
    # name_read_current=a[5]":"a[6]":"a[7]
    # delete a

    if (and($2, 0x400)){ # skip duplicates

        if (NODUPS_ARG == 1){
            next
        }else{
            is_duplicate_current="T"
        }

    }else{

        is_duplicate_current="F"
    }

    stats["total_reads_after_duplicate_removal"]++

    if (header == 1){ # first non-header record; change it so that header is extracted separately instead of header end being separately extracted for on every record

        header = 0
        is_duplicate_previous = is_duplicate_current
        print "@PG","ID:RMLizer.sh","VN:1.0","CL:"RMLizer_COMMAND > OUTDIR_ARG "/" BAM_NAME "_RMLizer.sam"

    }

    if ($6 ~ INDEL_BAN_ARG){

        next
    }

    stats["indel_filter_pass"]++

    nm_value = gensub(/NM:i:/, "", "g", tag_finder("NM:i:"))
    if ($5 < MAPQ_FILTER && nm_value > NM_FILTER){
        next
    }
    stats["reads_passing_quality_filters"]++

    if ($0 ~ /^.*XF:Z:.*$/){ # get XF tag field if present
        xf = tag_finder("XF:Z:")
        xf_value = gensub(/XF:Z:/, "", "g", xf)
        if (xf_value ~ /.+__.+/){ # only if gene IDs are appended with strandedness
            split(xf_value, xf_val, "__")
            gene_id = xf_val[1]; strand = xf_val[2]
            delete xf_val
        }else{
            gene_id = xf_value
        }
    }
    # implement counting of non-unique mapping!
    if (xf ~ /_ambiguous|_no_feature|_too_low_aQual|_not_aligned|_alignment_not_unique/){
        next
    }
    stats["reads_uniquely_mapped_to_features"]++

    ############################## Determine strand ####################################
    if ((and($2, 0x20) && and($2, 0x80)) || (and($2, 0x10) && and($2, 0x40)) || $2 == 0){
        strand2 = "forward"
        # stats["reads_with_reverse_strandedness"]++
    }else if ((and($2, 0x10) && and($2, 0x80)) || (and($2, 0x20) && and($2, 0x40)) || $2 == 16){
        strand2 = "reverse"
        # stats["reads_with_forward_strandedness"]++
    }else{
        strand2 = "undetermined"
        stats["reads_with_undetermined_strandedness"]++
        next
    }

    if ((strand == "+" && strand2 == "reverse") || (strand == "-" && strand2 == "forward")){ # make a conditional to print only first record like this
        stats["reads_with_undetermined_strandedness"]++
        printf "Warning: conflicting strandedness information between BAM and GTF files. Read skipped." > "/dev/stderr"
        print $0 > "/dev/stderr"
        next
    }

    if ($6 !~ /[M=X]/){ # only include reads with matches to the template; overlap clipping not used any longer, comment left as legacy
        next
    }

    md_value = gensub(/MD:Z:/, "", "g", tag_finder("MD:Z:"))

    if (md_value ~ /^[0-9]*(\^[ACTGN]+[0-9]+)*$/){
        stats["reads_without_mismatches"]++
        if ($6 ~ /^[^SDNI]+$/){ # CIGAR cases without I,N,D,S
            split($10, mmdel_read, "") # the variable naming with "del" infix preserved for consistency
            len_MDN = length($10)
        }else{
            MM_extractor_simple($10, $6)
        }
        ########## Counter function call for reads without mismatches ##############################
        if (strand == "+" || strand2 == "forward"){
            strand3 = "+"
            TCRA_producer_simple_for(name_read_current, name_read_previous, mmdel_read, mmdel_read_previous, len_MDN, len_MDN_previous)
        }else if (strand == "-" || strand2 == "reverse"){
            strand3 = "-"
            TCRA_producer_simple_rev(name_read_current, name_read_previous, mmdel_read, mmdel_read_previous, len_MDN, len_MDN_previous)
        }
        ############################################################################################

        print $0, "TC:i:0" > OUTDIR_ARG "/" BAM_NAME "_RMLizer.sam" # add data to SAM file

        if ($7 == "*"){

            pair = "FALSE"

        }else if (name_read_current != name_read_previous){

            if (pair == "FALSE"){ # for the first read in BAM pair == TRUE

                print chr_previous, fr_start_previous-1, fr_end_previous-1, strand3_previous, name_read_previous, gene_id_previous, count_previous["TC"]+0, count_previous["TA"] + count_previous["TC"] + count_previous["TG"] + count_previous["TT"] + 0, 0, 0, is_duplicate_previous > OUTDIR_ARG "/" BAM_NAME "_RMLizer.bed"

                print gene_id_previous,name_read_previous,count_previous["AA"]+0, count_previous["AC"]+0, count_previous["AG"]+0, count_previous["AT"]+0, count_previous["CA"]+0, count_previous["CC"]+0, count_previous["CG"]+0, count_previous["CT"]+0, count_previous["GA"]+0, count_previous["GC"]+0, count_previous["GG"]+0, count_previous["GT"]+0, count_previous["TA"]+0, count_previous["TC"]+0, count_previous["TG"]+0, count_previous["TT"]+0, "a", is_duplicate_previous > OUTDIR_ARG "/" BAM_NAME "_RMLizer_counts.txt"
            }

            pair = "FALSE"

        }else if (name_read_current == name_read_previous){

            # abs_span = $9
            # sub(/^-/,"",abs_span)
            # print $3, fr_start, fr_start + abs_span, strand3, $1, $2 > OUTDIR_ARG "/" BAM_NAME "_RMLizer.bed"

            print $3, fr_start-1, fr_end-1, strand3, $1, gene_id, count_previous["TC"]+0, count_previous["TA"] + count_previous["TC"] + count_previous["TG"] + count_all["TT"] + count_previous["TT"]+0, 0, count_overlap["TT"]+0, is_duplicate_previous > OUTDIR_ARG "/" BAM_NAME "_RMLizer.bed" # no mismatches in the current read, possibly mismatches in previous if mate

            print gene_id,$1,count_overlap["AA"]+0, 0, 0, 0, 0, count_overlap["CC"]+0, 0, 0, 0, 0, count_overlap["GG"]+0, 0, 0, 0, 0, count_overlap["TT"]+0, "o", is_duplicate_previous > OUTDIR_ARG "/" BAM_NAME "_RMLizer_counts.txt"

            print gene_id, $1, count_all["AA"] + count_previous["AA"]+0, count_previous["AC"]+0, count_previous["AG"]+0, count_previous["AT"]+0, count_previous["CA"]+0, count_all["CC"] + count_previous["CC"]+0, count_previous["CG"]+0, count_previous["CT"]+0, count_previous["GA"]+0, count_previous["GC"]+0, count_all["GG"] + count_previous["GG"]+0, count_previous["GT"]+0, count_previous["TA"]+0, count_previous["TC"]+0, count_previous["TG"]+0, count_all["TT"] + count_previous["TT"]+0, "a", is_duplicate_previous > OUTDIR_ARG "/" BAM_NAME "_RMLizer_counts.txt"

            pair = "TRUE"
        }

    }else{

        stats["reads_with_mismatches"]++

        MM_extractor($10, $11, $6, md_value)
        #######################################################
        if (strand == "+" || strand2 == "forward"){
            strand3 = "+"
            TCRA_producer_for(name_read_current, name_read_previous, mmdel_read, mmdel_read_previous, mmdel_quali, mmdel_md_ext, len_MDN, len_MDN_previous)
        }else if (strand == "-" || strand2 == "reverse"){
            strand3 = "-"
            TCRA_producer_rev(name_read_current, name_read_previous, mmdel_read, mmdel_read_previous, mmdel_quali, mmdel_md_ext, len_MDN, len_MDN_previous)
        }
        #######################################################
        if (xi < MISMATCH_FILTER){
            next
        }
        stats["reads_with_mismatches_passing_identity_filter"]++

        print $0,"TC:i:"count_all["TC"]+0 > OUTDIR_ARG "/" BAM_NAME "_RMLizer.sam" # add data to sam file

        if ($7 == "*"){

            pair = "FALSE"

        }else if (name_read_current != name_read_previous){

            if (pair == "FALSE"){

                print chr_previous, fr_start_previous-1, fr_end_previous-1, strand3_previous, name_read_previous, gene_id_previous, count_previous["TC"]+0, count_previous["TA"] + count_previous["TC"] + count_previous["TG"] + count_previous["TT"]+0, 0, 0, is_duplicate_previous > OUTDIR_ARG "/" BAM_NAME "_RMLizer.bed" # overlap not applicable

                print gene_id_previous, name_read_previous, count_previous["AA"]+0, count_previous["AC"]+0, count_previous["AG"]+0, count_previous["AT"]+0, count_previous["CA"]+0, count_previous["CC"]+0, count_previous["CG"]+0, count_previous["CT"]+0, count_previous["GA"]+0, count_previous["GC"]+0, count_previous["GG"]+0, count_previous["GT"]+0, count_previous["TA"]+0, count_previous["TC"]+0, count_previous["TG"]+0, count_previous["TT"]+0, "a", is_duplicate_previous > OUTDIR_ARG "/" BAM_NAME "_RMLizer_counts.txt"
            }

            pair = "FALSE"

        }else if (name_read_current == name_read_previous){

            # abs_span = $9
            # sub(/^-/,"",abs_span)
            # print $3, fr_start, fr_start + abs_span, strand3, $1, gene_id > OUTDIR_ARG "/" BAM_NAME "_RMLizer.bed"

            print $3, fr_start-1, fr_end-1, strand3, $1, gene_id, count_all["TC"] + count_previous["TC"]+0, count_all["TA"] + count_previous["TA"] + count_all["TC"] + count_previous["TC"] + count_all["TG"] + count_previous["TG"] + count_all["TT"] + count_previous["TT"]+0, count_overlap["TC"]+0, count_overlap["TA"] + count_overlap["TC"] + count_overlap["TG"] + count_overlap["TT"]+0, is_duplicate_previous > OUTDIR_ARG "/" BAM_NAME "_RMLizer.bed"

            print gene_id,$1,count_overlap["AA"]+0, count_overlap["AC"]+0, count_overlap["AG"]+0, count_overlap["AT"]+0, count_overlap["CA"]+0, count_overlap["CC"]+0, count_overlap["CG"]+0, count_overlap["CT"]+0, count_overlap["GA"]+0, count_overlap["GC"]+0, count_overlap["GG"]+0, count_overlap["GT"]+0, count_overlap["TA"]+0, count_overlap["TC"]+0, count_overlap["TG"]+0, count_overlap["TT"]+0, "o", is_duplicate_previous > OUTDIR_ARG "/" BAM_NAME "_RMLizer_counts.txt"

            print gene_id,$1,count_all["AA"]+count_previous["AA"]+0, count_all["AC"]+count_previous["AC"]+0, count_all["AG"]+count_previous["AG"]+0, count_all["AT"]+count_previous["AT"]+0, count_all["CA"]+count_previous["CA"]+0, count_all["CC"]+count_previous["CC"]+0, count_all["CG"]+count_previous["CG"]+0, count_all["CT"]+count_previous["CT"]+0, count_all["GA"]+count_previous["GA"]+0, count_all["GC"]+count_previous["GC"]+0, count_all["GG"]+count_previous["GG"]+0, count_all["GT"]+count_previous["GT"]+0, count_all["TA"]+count_previous["TA"]+0, count_all["TC"]+count_previous["TC"]+0, count_all["TG"]+count_previous["TG"]+0, count_all["TT"]+count_previous["TT"]+0, "a", is_duplicate_previous > OUTDIR_ARG "/" BAM_NAME "_RMLizer_counts.txt"

            pair = "TRUE"
        }
    }
    ################
    if (os >= 1){

        stats["pairs_with_overlaps"]++

    }
    name_read_previous = name_read_current
    is_duplicate_previous = is_duplicate_current
    gene_id_previous = gene_id
    strand3_previous = strand3
    fr_start_previous = fr_start
    fr_end_previous = fr_end
    chr_previous = $3
    len_MDN_previous = len_MDN
    # flag_previous = $2
    copy_array(count_all, count_previous)
    copy_array(mmdel_read, mmdel_read_previous)

} # end of main routine

END{

    if (pair == "FALSE"){

        print gene_id_previous, name_read_previous, count_previous["AA"]+0, count_previous["AC"]+0, count_previous["AG"]+0, count_previous["AT"]+0, count_previous["CA"]+0, count_previous["CC"]+0, count_previous["CG"]+0, count_previous["CT"]+0, count_previous["GA"]+0, count_previous["GC"]+0, count_previous["GG"]+0, count_previous["GT"]+0, count_previous["TA"]+0, count_previous["TC"]+0, count_previous["TG"]+0, count_previous["TT"]+0, "a", is_duplicate_previous > OUTDIR_ARG "/" BAM_NAME "_RMLizer_counts.txt"
    }

    timestamp=systime()

    for (i in depth){
        print i, depth[i] > OUTDIR_ARG "/" BAM_NAME "_depth_per_position_" timestamp ".txt"
    }
    delete depth

    for (i in breadth){
        print i, breadth[i] > OUTDIR_ARG "/" BAM_NAME "_breadth_per_gene_" timestamp ".txt"
    }
    delete breadth

    printf "Total reads: %s\nReads after deduplication: %s\nReads passing indel filter: %s\nReads passing quality filters: %s\nReads uniquely mapped to features: %s\nSkipped reads due to undetermined strandedness: %s\nReads without mismatches: %s\nReads with mismatches: %s\nReads passing alignment identity filter: %s\nNumber of read pairs with overlap of at least 1 base: %s\nCumulative sum of overlapping bases: %s\n", TOTAL_READS+0, stats["total_reads_after_duplicate_removal"]+0, stats["indel_filter_pass"]+0, stats["reads_passing_quality_filters"]+0, stats["reads_uniquely_mapped_to_features"]+0, stats["reads_with_undetermined_strandedness"]+0, stats["reads_without_mismatches"]+0, stats["reads_with_mismatches"]+0, stats["reads_with_mismatches_passing_identity_filter"]+0, stats["pairs_with_overlaps"]+0, to+0 >> OUTDIR_ARG "/" BAM_NAME "_RMLizer_stats.log"

    printf "--> %s out of %s records (%.2f\%) processed\n", NR-header_lines, TOTAL_READS, ((NR-header_lines)/TOTAL_READS)*100

    finish_time = systime()

    printf "S(L)AM converter ended @ %s\nRun time: %.0f hrs which equals %.0f min\n", strftime(), (finish_time - start_time)/3600, (finish_time - start_time)/60 > "/dev/stderr"

}
