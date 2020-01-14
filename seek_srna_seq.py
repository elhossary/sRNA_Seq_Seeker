# Author: Muhammad Elhossary | elhossary@zbmed.de
import sys
from numpy import genfromtxt
import argparse
import glob


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tss_in", required=True, type=str, help="")
    parser.add_argument("--term_in", required=True, type=str, help="")
    parser.add_argument("--gff_out", required=True, type=str, help="")
    parser.add_argument("--max_len", required=True, type=int, help="")
    parser.add_argument("--min_len", required=True, type=int, help="")
    args = parser.parse_args()
    tss_arr = build_arr_form_gff(glob.glob(args.tss_in)[0])
    term_arr = build_arr_form_gff(glob.glob(args.term_in)[0])
    output_file_path = glob.glob(args.gff_out)[0]
    srna_max_length = int(args.max_len)
    srna_min_length = int(args.min_len)
    print("\n\n--- sRNA Seq Seeker ---\n\n")
    print(f"Seeking for possible sRNA at sequences terminated at maximum length of {srna_max_length}"
          f" and length between tss and terminator not shorter that {srna_min_length}")
    srna_gff_str = find_possible_sRNA(srna_max_length, tss_arr, term_arr, srna_min_length)
    print("\nWriting output to file")
    outfile = open(output_file_path, "w")
    outfile.write(f"###gff-version 3\n{srna_gff_str}###")
    outfile.close()
    print("DONE")


def find_possible_sRNA(srna_max_length, tss_arr, term_arr, srna_min_length):
    # Number of Column names
    # 0 = accession
    # 1 = source
    # 2 = type
    # 3 = start
    # 4 = end
    # 5 = dot1
    # 6 = strand
    # 7 = dot2
    # 8 = attributes

    r_srna_gff_str = ""
    srna_count = 0
    tss_arr_len = len(tss_arr)
    for tss_index, tss_row in enumerate(tss_arr):
        sys.stdout.flush()
        sys.stdout.write("\r" + f"Progress: {round(tss_index / tss_arr_len * 100, 2)}% | " +
                         f"{srna_count} possible sRNAs found      ...")
        for term_index, term_row in enumerate(term_arr):
            if tss_row[0] == term_row[0]:
                if tss_row[6] == term_row[6] == "+":
                    if tss_row[4] < term_row[3] and \
                            (term_row[4] - tss_row[3]) <= srna_max_length and \
                            srna_min_length <= (term_row[3] - tss_row[3]):
                        srna_count += 1
                        r_srna_gff_str += \
                            f"{tss_row[0]}\t" + \
                            f"sRNA_Seq_Seeker\t" + \
                            f"possible_sRNA_seq\t" + \
                            f"{tss_row[3]}\t" + \
                            f"{term_row[4]}\t" + \
                            f".\t" + \
                            f"{term_row[6]}\t" + \
                            f".\t" + \
                            f"id=possible_srna{srna_count};" + \
                            f"name=possible_srna{srna_count};" + \
                            f"seq_len={term_row[4] - tss_row[3]};" + \
                            f"matched_tss={parse_attributes(tss_row[8])['id']};" + \
                            f"matched_terminator={parse_attributes(term_row[8])['id']}\n"

                if tss_row[6] == term_row[6] == "-":
                    if term_row[4] < tss_row[3] and \
                            (tss_row[4] - term_row[3]) <= srna_max_length and \
                            srna_min_length <= (tss_row[3] - term_row[3]):
                        srna_count += 1
                        r_srna_gff_str += \
                            f"{tss_row[0]}\t" + \
                            f"sRNA_Seq_Seeker\t" + \
                            f"possible_sRNA_seq\t" + \
                            f"{term_row[3]}\t" + \
                            f"{tss_row[4]}\t" + \
                            f".\t" + \
                            f"{term_row[6]}\t" + \
                            f".\t" + \
                            f"id=possible_srna{srna_count};" + \
                            f"name=possible_srna{srna_count};" + \
                            f"seq_len={tss_row[4] - term_row[3]};" + \
                            f"matched_tss={parse_attributes(tss_row[8])['id']};" + \
                            f"matched_terminator={parse_attributes(term_row[8])['id']}\n"

    sys.stdout.write("\r" + f"Progress 100% with total {srna_count} possible sRNAs could be found")
    return r_srna_gff_str


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


def build_arr_form_gff(path):
    data_arr = genfromtxt(path, delimiter="\t", comments="#", dtype=None, encoding=None)
    return data_arr





main()
