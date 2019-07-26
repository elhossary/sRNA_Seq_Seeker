# Author: Muhammad Elhossary | elhossary@zbmed.de
import sys
from numpy import genfromtxt


def parse_attributes(attr_str):
    return dict(item.split("=") for item in attr_str.split(";"))


def build_arr_form_gff(path):
    data_arr = genfromtxt(path, delimiter="\t", comments="#", dtype=None, encoding=None)
    return data_arr

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


def find_possible_sRNA(srna_max_length, tss_arr, term_arr):
    r_srna_gff_str = ""
    srna_count = 0
    tss_arr_len = len(tss_arr)
    for tss_index, tss_row in enumerate(tss_arr):
        sys.stdout.flush()
        sys.stdout.write("\r" + f"Progress: {round(tss_index / tss_arr_len * 100, 2)}% | " +
                         f"{srna_count } possible sRNAs found      ...")
        for term_index, term_row in enumerate(term_arr):
            if tss_row[0] == term_row[0]:
                if tss_row[6] == term_row[6] == "+":
                    if tss_row[4] < term_row[3] and \
                            (term_row[4] - tss_row[3]) <= srna_max_length:
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
                            f"matched_tss={parse_attributes(term_row[8])['Name']};" + \
                            f"matched_terminator={parse_attributes(term_row[8])['Name']}\n"
                if tss_row[6] == term_row[6] == "-":
                    if term_row[4] < tss_row[3] and \
                            (tss_row[4] - term_row[3]) <= srna_max_length:
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
                            f"matched_tss={parse_attributes(term_row[8])['Name']};" + \
                            f"matched_terminator={parse_attributes(term_row[8])['Name']}\n"
    sys.stdout.write("\r" + f"Progress 100% with total {srna_count} possible sRNAs could be found")
    return r_srna_gff_str


if sys.argv[1] is None or sys.argv[2] is None or sys.argv[3] is None or int(sys.argv[4]) is None:
    print("check your inputs")
    exit()
else:
    tss_arr = build_arr_form_gff(sys.argv[1])
    term_arr = build_arr_form_gff(sys.argv[2])
    output_file_path = sys.argv[3]
    srna_max_length = int(sys.argv[4])
    print("\n\n--- sRNA Seq Seeker ---\n\n")
    print(f"Seeking for possible sRNA at sequence length of {srna_max_length} nucleotides")
    srna_gff_str = find_possible_sRNA(srna_max_length, tss_arr, term_arr)
    print("\nWriting output to file")
    outfile = open(output_file_path, "w")
    outfile.write(f"###gff-version 3\n{srna_gff_str}\n###")
    outfile.close()
    print("DONE")
