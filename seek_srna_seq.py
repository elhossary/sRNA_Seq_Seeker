# Author: Muhammad Elhossary | elhossary@zbmed.de
import pandas as pd
import sys
import time

def parse_attributes(attr_str):
    return dict(item.split("=") for item in attr_str.split(";"))

def build_df_form_gff(path):
    in_file = open(path, "r")
    r_df = pd.read_csv(in_file, sep="\t", comment="#", header=None)
    in_file.close()
    r_df[3] = pd.to_numeric(r_df[3], downcast='integer')
    r_df[4] = pd.to_numeric(r_df[4], downcast='integer')
    r_df.sort_values(by=[3])
    return r_df

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

def find_possible_sRNA(srna_max_length, tss_df, term_df):
    r_srna_gff_str = ""
    srna_count = 0
    tss_df_len = len(tss_df.index)
    term_df_len = len(term_df.index)
    for tss_row_index in range(0, tss_df_len, 1):
        for term_row_index in range(0, term_df_len, 1):
            if tss_df.iloc[tss_row_index, 0] == term_df.iloc[term_row_index, 0] and \
                    tss_df.iloc[tss_row_index, 6] == term_df.iloc[term_row_index, 6] and \
                    term_df.iloc[term_row_index, 3] > tss_df.iloc[tss_row_index, 4] and \
                    (term_df.iloc[term_row_index, 4] - tss_df.iloc[tss_row_index, 3]) <= srna_max_length:
                srna_count += 1
                sys.stdout.write('\r' + f"Progress: {round(tss_row_index / tss_df_len * 100, 2)}%, " +
                                 f"Possible sRNAs found: {srna_count}")
                time.sleep(1)
                sys.stdout.flush()
                r_srna_gff_str += f"{tss_df.iloc[tss_row_index, 0]}\t" + \
                                  f"sRNA_Seq_Seeker\t" + \
                                  f"possible_sRNA_seq\t" + \
                                  f"{tss_df.iloc[tss_row_index, 4]}\t" + \
                                  f"{term_df.iloc[term_row_index, 4]}\t" + \
                                  f".\t" + \
                                  f"{term_df.iloc[term_row_index, 6]}\t" + \
                                  f".\t" + \
                                  f"id=possible_srna{srna_count};" + \
                                  f"name=possible_srna{srna_count};" + \
                                  f"seq_len={term_df.iloc[term_row_index, 4] - tss_df.iloc[tss_row_index, 3]};" + \
                                  f"matched_tss={parse_attributes(tss_df.iloc[term_row_index, 8])['Name']};" + \
                                  f"matched_terminator={parse_attributes(term_df.iloc[term_row_index, 8])['Name']}\n"
    return r_srna_gff_str


if sys.argv[1] is None or sys.argv[2] is None or sys.argv[3] is None or int(sys.argv[4]) is None:
    print("check your inputs")
    exit()
else:
    tss_df = build_df_form_gff(sys.argv[1])
    term_df = build_df_form_gff(sys.argv[2])
    output_file_path = sys.argv[3]
    srna_max_length = int(sys.argv[4])

    print(f"Seeking for possible sRNA at sequence length of {srna_max_length}")
    srna_gff_str = find_possible_sRNA(srna_max_length, tss_df, term_df)
    print("\nWriting output to file")
    outfile = open(output_file_path, "w")
    outfile.write(f"###gff-version 3\n{srna_gff_str}")
    outfile.close()
    print("DONE")
