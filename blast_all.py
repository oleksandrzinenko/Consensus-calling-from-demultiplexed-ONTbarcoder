
#!/usr/bin/env python3

import os
import csv
import time
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# ============================================================
# SETTINGS
# ============================================================

INPUT_FASTA = "all_consensus.fasta"
OUTPUT_CSV = "blast_results.csv"

PROGRAM = "blastn"
DATABASE = "nt"

MAX_HITS = 10

# NCBI identification
# CHANGE THIS TO YOUR EMAIL ADDRESS
EMAIL = "oleksandrzinenko@gmail.com"
TOOL = "ONT_consensus_BLAST"

# Number of attempts if network/NCBI fails
MAX_RETRIES = 3

# Wait between different BLAST submissions
# NCBI recommends at least 10 seconds
SUBMISSION_DELAY = 12


# ============================================================
# SET NCBI IDENTIFICATION
# ============================================================

NCBIWWW.email = EMAIL
NCBIWWW.tool = TOOL


# ============================================================
# FIND ALREADY COMPLETED QUERIES
# ============================================================

completed = set()

if os.path.exists(OUTPUT_CSV):

    with open(OUTPUT_CSV, newline="") as f:

        reader = csv.DictReader(f)

        for row in reader:
            completed.add(row["query_id"])

    print(f"Existing results found: {len(completed)} queries completed")


# ============================================================
# CREATE CSV IF NEEDED
# ============================================================

new_file = not os.path.exists(OUTPUT_CSV)

csvfile = open(OUTPUT_CSV, "a", newline="")
writer = csv.writer(csvfile)

if new_file:

    writer.writerow([
        "query_id",
        "query_length",
        "rank",
        "accession",
        "title",
        "identity_percent",
        "alignment_length",
        "query_coverage_percent",
        "mismatches",
        "gaps",
        "evalue",
        "bitscore"
    ])

    csvfile.flush()


# ============================================================
# READ FASTA
# ============================================================

records = list(SeqIO.parse(INPUT_FASTA, "fasta"))

total = len(records)

print()
print("================================================")
print("NCBI BLAST batch search")
print("================================================")
print(f"Input file:        {INPUT_FASTA}")
print(f"Total sequences:   {total}")
print(f"Already completed: {len(completed)}")
print(f"Remaining:         {total - len(completed)}")
print()


# ============================================================
# BLAST EACH SEQUENCE
# ============================================================

for number, record in enumerate(records, start=1):

    query_id = record.id
    sequence = str(record.seq)

    # --------------------------------------------------------
    # SKIP COMPLETED QUERIES
    # --------------------------------------------------------

    if query_id in completed:

        print(
            f"[{number}/{total}] "
            f"{query_id} -> already completed, skipping"
        )

        continue


    print()
    print("------------------------------------------------")
    print(f"[{number}/{total}] BLASTing:")
    print(f"Query:  {query_id}")
    print(f"Length: {len(sequence)} bp")
    print("------------------------------------------------")


    # --------------------------------------------------------
    # SUBMIT WITH RETRIES
    # --------------------------------------------------------

    success = False

    for attempt in range(1, MAX_RETRIES + 1):

        try:

            print(
                f"Submitting to NCBI "
                f"(attempt {attempt}/{MAX_RETRIES})..."
            )

            result_handle = NCBIWWW.qblast(
                program=PROGRAM,
                database=DATABASE,
                sequence=sequence,
                hitlist_size=MAX_HITS,
                format_type="XML"
            )

            xml_data = result_handle.read()

            result_handle.close()

            success = True

            print("BLAST completed.")

            break


        except Exception as e:

            print()
            print(f"ERROR: {e}")
            print()

            if attempt < MAX_RETRIES:

                wait = 60 * attempt

                print(
                    f"Waiting {wait} seconds before retry..."
                )

                time.sleep(wait)


    # --------------------------------------------------------
    # FAILED AFTER ALL RETRIES
    # --------------------------------------------------------

    if not success:

        print(
            f"FAILED: {query_id}"
        )

        writer.writerow([
            query_id,
            len(sequence),
            "ERROR",
            "",
            "BLAST_FAILED",
            "",
            "",
            "",
            "",
            "",
            "",
            ""
        ])

        csvfile.flush()

        continue


    # --------------------------------------------------------
    # SAVE INDIVIDUAL XML
    # --------------------------------------------------------

    xml_filename = f"blast_xml_{number}.xml"

    with open(xml_filename, "w") as f:

        f.write(xml_data)


    # --------------------------------------------------------
    # PARSE XML
    # --------------------------------------------------------

    try:

        with open(xml_filename) as result_handle:

            blast_record = NCBIXML.read(result_handle)


    except Exception as e:

        print(f"ERROR parsing XML: {e}")

        writer.writerow([
            query_id,
            len(sequence),
            "ERROR",
            "",
            "XML_PARSE_FAILED",
            "",
            "",
            "",
            "",
            "",
            "",
            ""
        ])

        csvfile.flush()

        continue


    # --------------------------------------------------------
    # EXTRACT TOP HITS
    # --------------------------------------------------------

    rank = 0

    if len(blast_record.alignments) == 0:

        print("No BLAST hits found.")

        writer.writerow([
            query_id,
            len(sequence),
            0,
            "",
            "NO_HIT",
            "",
            "",
            "",
            "",
            "",
            "",
            ""
        ])

        csvfile.flush()

    else:

        for alignment in blast_record.alignments:

            for hsp in alignment.hsps:

                rank += 1

                if rank > MAX_HITS:
                    break


                identity_percent = (
                    100 * hsp.identities /
                    hsp.align_length
                )


                # Approximate query coverage
                query_coverage = (
                    100 * hsp.align_length /
                    len(sequence)
                )


                mismatches = (
                    hsp.align_length -
                    hsp.identities -
                    hsp.gaps
                )


                writer.writerow([

                    query_id,

                    len(sequence),

                    rank,

                    alignment.accession,

                    alignment.title,

                    round(identity_percent, 3),

                    hsp.align_length,

                    round(query_coverage, 2),

                    mismatches,

                    hsp.gaps,

                    hsp.expect,

                    hsp.bits

                ])

                csvfile.flush()


            if rank >= MAX_HITS:

                break


        print(
            f"Saved {rank} BLAST hits."
        )


    # --------------------------------------------------------
    # WAIT BEFORE NEXT QUERY
    # --------------------------------------------------------

    print(
        f"Waiting {SUBMISSION_DELAY} seconds "
        f"before next query..."
    )

    time.sleep(SUBMISSION_DELAY)


# ============================================================
# FINISH
# ============================================================

csvfile.close()

print()
print("================================================")
print("ALL DONE")
print("================================================")
print()
print(f"Results saved to:")
print(f"  {OUTPUT_CSV}")

