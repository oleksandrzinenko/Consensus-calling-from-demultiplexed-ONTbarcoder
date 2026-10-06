#!/bin/bash

# ============================================================
# ONTbarcoder 2.0 -> VSEARCH clustering -> centroid + consensus
# ============================================================

IDENTITY=0.80
MIN_READS=3

# Output directories
mkdir -p centroids consensus msa cleaned

# Combined output files
ALL_CENTROIDS="all_centroids.fasta"
ALL_CONSENSUS="all_consensus.fasta"
SUMMARY="cluster_summary.tsv"

# Start fresh
rm -f "$ALL_CENTROIDS" "$ALL_CONSENSUS" "$SUMMARY"

echo -e "Sample\tInput_records\tDNA_reads\tNonDNA_reads\tClusters\tRetained_clusters\tReads_in_retained" \
    > "$SUMMARY"


for f in *_all.fa; do

    sample="${f%_all.fa}"

    clean="cleaned/${sample}_clean.fa"
    centroid="centroids/${sample}_centroids.fa"
    consensus="consensus/${sample}_consensus.fa"
    msa="msa/${sample}_msa.fa"

    echo
    echo "============================================================"
    echo "Processing: $f"
    echo "============================================================"

    # --------------------------------------------------------
    # STEP 1
    # Extract the sequence field from each ONTbarcoder record.
    #
    # Keep only sequences consisting exclusively of A,C,G,T,N.
    # This removes EPFP... records but does NOT alter genuine DNA.
    # --------------------------------------------------------

    awk '
    {
        # The sequence is the final whitespace-separated field
        seq=$NF

        if (seq ~ /^[ACGTNacgtn]+$/) {
            print ">" ++n
            print seq
            dna++
        }
        else {
            nondna++
        }
    }

    END {
        print n > "/tmp/dna_count.txt"
        print nondna > "/tmp/nondna_count.txt"
    }
    ' "$f" > "$clean"

    dna_reads=$(cat /tmp/dna_count.txt)
    nondna_reads=$(cat /tmp/nondna_count.txt)

    input_records=$((dna_reads + nondna_reads))

    echo "Input records:       $input_records"
    echo "DNA reads:           $dna_reads"
    echo "Non-DNA reads:       $nondna_reads"

    # --------------------------------------------------------
    # STEP 2
    # Cluster DNA reads
    # --------------------------------------------------------

    vsearch --cluster_fast "$clean" \
        --id "$IDENTITY" \
        --strand both \
        --centroids "$centroid" \
        --sizeout \
        --consout "$consensus" \
        --msaout "$msa"

    # --------------------------------------------------------
    # STEP 3
    # Count clusters
    # --------------------------------------------------------

    clusters=$(grep -c '^>' "$centroid")

    echo "Clusters:             $clusters"

    # --------------------------------------------------------
    # STEP 4
    # Extract clusters >= MIN_READS
    #
    # Do this separately for centroid and consensus.
    # --------------------------------------------------------

    retained=0
    retained_reads=0

    awk \
        -v sample="$sample" \
        -v minreads="$MIN_READS" \
        -v output="$ALL_CENTROIDS" '

    BEGIN {
        RS=">"
    }

    NR > 1 {

        n=split($0,line,"\n")
        header=line[1]

        reads=0

        if (match(header,/size=[0-9]+/)) {
            x=substr(header,RSTART,RLENGTH)
            sub(/^size=/,"",x)
            reads=x+0
        }

        if (reads >= minreads) {

            kept++

            printf ">%s_cluster_%d_reads_%d\n",
                   sample, kept, reads >> output

            for(i=2;i<=n;i++)
                if(line[i]!="")
                    print line[i] >> output

            print "" >> output

            total_reads+=reads
        }
    }

    END {
        print kept > "/tmp/kept.txt"
        print total_reads > "/tmp/kept_reads.txt"
    }

    ' "$centroid"

    retained=$(cat /tmp/kept.txt)
    retained_reads=$(cat /tmp/kept_reads.txt)

    # --------------------------------------------------------
    # STEP 5
    # Process consensus FASTA.
    #
    # Consensus headers normally contain seqs=N.
    # --------------------------------------------------------

    awk \
        -v sample="$sample" \
        -v minreads="$MIN_READS" \
        -v output="$ALL_CONSENSUS" '

    BEGIN {
        RS=">"
    }

    NR > 1 {

        n=split($0,line,"\n")
        header=line[1]

        reads=0

        # VSEARCH consensus: seqs=N
        if (match(header,/seqs=[0-9]+/)) {
            x=substr(header,RSTART,RLENGTH)
            sub(/^seqs=/,"",x)
            reads=x+0
        }

        # Some VSEARCH versions may use size=N
        if (reads==0 && match(header,/size=[0-9]+/)) {
            x=substr(header,RSTART,RLENGTH)
            sub(/^size=/,"",x)
            reads=x+0
        }

        if (reads >= minreads) {

            consensus_n++

            printf ">%s_cluster_%d_reads_%d\n",
                   sample, consensus_n, reads >> output

            for(i=2;i<=n;i++)
                if(line[i]!="")
                    print line[i] >> output

            print "" >> output
        }
    }
    ' "$consensus"

    # --------------------------------------------------------
    # STEP 6
    # Summary
    # --------------------------------------------------------

    echo -e "$sample\t$input_records\t$dna_reads\t$nondna_reads\t$clusters\t$retained\t$retained_reads" \
        >> "$SUMMARY"

    echo "Retained >= $MIN_READS: $retained"
    echo "Reads in retained:     $retained_reads"

done


rm -f /tmp/dna_count.txt \
      /tmp/nondna_count.txt \
      /tmp/kept.txt \
      /tmp/kept_reads.txt

echo
echo "============================================================"
echo "DONE"
echo "============================================================"

echo "Centroids:"
echo "  $ALL_CENTROIDS"

echo "Consensus:"
echo "  $ALL_CONSENSUS"

echo "MSA files:"
echo "  msa/"

echo "Summary:"
echo "  $SUMMARY"

echo
echo "Centroids in combined file:"
grep -c '^>' "$ALL_CENTROIDS"

echo "Consensus sequences in combined file:"
grep -c '^>' "$ALL_CONSENSUS"
