# Author: Muhammad Elhossary | elhossary@zbmed.de
import sys
from numpy import genfromtxt
import argparse
import glob
from io import StringIO
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tss_in", required=True, type=str, help="")
    parser.add_argument("--term_in", required=True, type=str, help="")
    parser.add_argument("--gff_out", required=True, type=str, help="")
    parser.add_argument("--max_len", required=True, type=int, help="")
    parser.add_argument("--min_len", required=True, type=int, help="")
    parser.add_argument("--merge_overlaps", action='store_true', help="")
    args = parser.parse_args()
    tss_arr = build_arr_form_gff(glob.glob(args.tss_in)[0])
    term_arr = build_arr_form_gff(glob.glob(args.term_in)[0])

    print("\n\n--- sRNA Seq Seeker ---\n\n")
    print(f"Seeking for possible sRNA at sequences at length between {args.min_len} and {args.max_len}")
    srna_gff_str = find_possible_sRNA(args.max_len, tss_arr, term_arr, args.min_len)
    if args.merge_overlaps:
        srna_gff_str = merge_overlaps(srna_gff_str)
    print("\nWriting output to file")
    outfile = open(args.gff_out, "w")
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


def merge_overlaps(srna_gff_str):
    col_names = ["accession", "source", "type", "start", "end", "dot1", "strand", "dot2", "attributes"]
    ret_srna_gff_str = ""
    sran_gff_df = pd.read_csv(StringIO(srna_gff_str), names=col_names, sep="\t", comment="#")
    accession_list = list(sran_gff_df.accession.unique())
    df_dict = {}
    for acc in accession_list:
        df_dict[f"{acc}_f"] = \
            merge_interval_lists(sran_gff_df[(sran_gff_df['accession'] == acc) & (sran_gff_df['strand'] == "+")]
                                 .loc[:, ['start', 'end']].sort_values(by=['start']).values.tolist())
        df_dict[f"{acc}_r"] = \
            merge_interval_lists(sran_gff_df[(sran_gff_df['accession'] == acc) & (sran_gff_df['strand'] == "-")]
                                 .loc[:, ['start', 'end']].sort_values(by=['start']).values.tolist())
    strand_func = lambda x: "+" if "_f" in x else "-"

    for acc in accession_list:
        srna_count = 0
        for dict_key in df_dict.keys():
            if dict_key == f"{acc}_f" or dict_key == f"{acc}_r":
                for loc in df_dict[dict_key]:
                    srna_count += 1
                    ret_srna_gff_str += \
                        f"{acc}\t" + \
                        f"sRNA_Seq_Seeker\t" + \
                        f"merged_possible_sRNA_seq\t" + \
                        f"{loc[0]}\t" + \
                        f"{loc[1]}\t" + \
                        f".\t" + \
                        f"{strand_func(dict_key)}\t" + \
                        f".\t" + \
                        f"id={dict_key}_possible_srna_{srna_count};" + \
                        f"name={dict_key}_possible_srna_{srna_count};" + \
                        f"seq_len={loc[1] - loc[0]}\n"
    ret_sran_gff_df = pd.read_csv(StringIO(ret_srna_gff_str), names=col_names, sep="\t", comment="#")
    ret_sran_gff_df = ret_sran_gff_df.sort_values(by=['accession', 'start'])
    ret_srna_gff_str = ret_sran_gff_df.to_csv(sep="\t", index=False, header=False)
    print(f"Total sRNAs after merge: {ret_sran_gff_df.shape[0]}")
    return ret_srna_gff_str


def merge_interval_lists(list_in, merge_range=0):
    list_out = []
    for loc in list_in:
        if len(list_out) == 0:
            list_out.append(loc)
        else:
            if loc[0] in range(list_out[-1][0], list_out[-1][-1] + merge_range):
                list_out[-1][-1] = loc[-1]
            else:
                list_out.append(loc)
    return list_out


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


def build_arr_form_gff(path):
    data_arr = genfromtxt(path, delimiter="\t", comments="#", dtype=None, encoding=None)
    return data_arr





main()
