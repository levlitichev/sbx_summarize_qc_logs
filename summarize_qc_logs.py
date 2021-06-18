"""
summarize_qc_logs.py

This function finds and summarizes log files made during QC steps of Sunbeam. Output is a single text file.

Currently only works for paired-end data.

"""

import pandas as pd
import argparse
import os
import sys
import glob
import re

__author__ = "Lev Litichevskiy"
__email__ = "levlitichev@gmail.com"

CUTADAPT_COLUMN_NAMES = ["filename", "cutadapt_read_pairs_processed", "cutadapt_read_pairs_written",
                         "cutadapt_bps_processed", "cutadapt_bps_written"]
TRIMMOMATIC_COLUMN_NAMES = ["filename", "trimmomatic_read_pairs_processed", "trimmomatic_read_pairs_written"]
KOMPLEXITY_COLUMN_NAMES = ["filename", "komplexity_read_pairs_removed"]
DECONTAM_COLUMN_NAMES = ["filename", "decontam_read_pairs_written", "decontam_read_pairs_removed"]


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--qc_log_dir", "-i", required=True, type=str,
                        help="path to sunbeam_output/qc/log directory")
    parser.add_argument("--is_single_end", "-s", action="store_true", default=False,
                        help="indicate whether data is single-end, default is paired-end")
    parser.add_argument("--out_path", "-o", type=str, default="qc_summary.txt",
                        help="path to output file")

    return parser


def get_cutadapt_paths(qc_log_dir):
    cutadapt_wildcard = os.path.join(qc_log_dir, "cutadapt", "*.log")
    cutadapt_paths = glob.glob(os.path.expanduser(cutadapt_wildcard))

    assert len(cutadapt_paths) > 0, f"No cutadapt files found with the following wildcard: {cutadapt_wildcard}"

    return cutadapt_paths


def cutadapt_summary(cutadapt_paths, is_single_end=False):

    list_of_lists = [""] * len(cutadapt_paths)

    for ii, this_path in enumerate(cutadapt_paths):
        this_filename = os.path.splitext(os.path.basename(this_path))[0]
        this_file_handle = open(this_path, "rt")

        for line in this_file_handle:
            if re.search("Total read pairs processed", line):
                this_read_pairs_processed = int(line.split()[-1].replace(",", ""))
            elif re.search("Pairs written", line):
                this_read_pairs_written = int(line.split()[-2].replace(",", ""))
            elif re.search("Total basepairs processed", line):
                this_bps_processed = int(line.split()[-2].replace(",", ""))
            elif re.search("Total written", line):
                this_bps_written = int(line.split()[-3].replace(",", ""))

                # no need to read the rest of the file
                break
            else:
                pass
        this_file_handle.close()

        list_of_lists[ii] = [this_filename,
                             this_read_pairs_processed, this_read_pairs_written,
                             this_bps_processed, this_bps_written]

    out_df = pd.DataFrame(list_of_lists, columns=CUTADAPT_COLUMN_NAMES)
    return out_df


def get_trimmomatic_paths(qc_log_dir):
    trimmomatic_wildcard = os.path.join(qc_log_dir, "trimmomatic", "*.out")
    trimmomatic_paths = glob.glob(os.path.expanduser(trimmomatic_wildcard))

    assert len(trimmomatic_paths) > 0, f"No trimmomatic files found with the following wildcard: {trimmomatic_wildcard}"

    return trimmomatic_paths


def trimmomatic_summary(trimmomatic_paths, is_single_end=False):

    list_of_lists = [""] * len(trimmomatic_paths)

    for ii, this_path in enumerate(trimmomatic_paths):
        this_filename = os.path.splitext(os.path.basename(this_path))[0]
        this_file_handle = open(this_path, "rt")

        for line in this_file_handle:
            if re.search("Input Read Pairs", line):
                line_split_up = line.split()
                this_read_pairs_processed = int(line_split_up[3])
                this_read_pairs_written = int(line_split_up[6])
                break
        this_file_handle.close()

        list_of_lists[ii] = [this_filename, this_read_pairs_processed, this_read_pairs_written]

    out_df = pd.DataFrame(list_of_lists, columns=TRIMMOMATIC_COLUMN_NAMES)
    return out_df


def get_komplexity_paths(qc_log_dir):
    komplexity_wildcard = os.path.join(qc_log_dir, "komplexity", "*.filtered_ids")
    komplexity_paths = glob.glob(os.path.expanduser(komplexity_wildcard))

    assert len(komplexity_paths) > 0, f"No komplexity files found with the following wildcard: {komplexity_wildcard}"

    return komplexity_paths


def komplexity_summary(komplexity_paths, is_single_end=False):

    list_of_lists = [""] * len(komplexity_paths)

    for ii, this_path in enumerate(komplexity_paths):
        this_filename = os.path.splitext(os.path.basename(this_path))[0]
        this_file_handle = open(this_path, "rt")
        # this_read_pairs_removed = len(this_file_handle.readlines())
        this_read_pairs_removed = 0

        for _ in this_file_handle:
            this_read_pairs_removed += 1

        this_file_handle.close()

        list_of_lists[ii] = [this_filename, this_read_pairs_removed]

    out_df = pd.DataFrame(list_of_lists, columns=KOMPLEXITY_COLUMN_NAMES)
    return out_df


def get_decontam_paths(qc_log_dir):
    decontam_wildcard = os.path.join(qc_log_dir, "decontam", "*.txt")
    decontam_paths = glob.glob(os.path.expanduser(decontam_wildcard))

    assert len(decontam_paths) > 0, f"No decontam files found with the following wildcard: {decontam_wildcard}"

    return decontam_paths


def decontam_summary(decontam_paths, is_single_end=False):

    list_of_lists = [""] * len(decontam_paths)

    for ii, this_path in enumerate(decontam_paths):
        this_filename = os.path.splitext(os.path.basename(this_path))[0]

        in_df = pd.read_table(this_path)

        this_read_pairs_written = int(in_df.loc[:, "nonhost"][0])
        this_read_pairs_removed = int(in_df.loc[:, "host"][0])

        list_of_lists[ii] = [this_filename,
                             this_read_pairs_written, this_read_pairs_removed]

    tmp_out_df = pd.DataFrame(list_of_lists, columns=DECONTAM_COLUMN_NAMES)

    # if paired-end, need to sum the two files corresponding to the same sample
    if not is_single_end:
        out_df = tmp_out_df
        out_df["filename"] = tmp_out_df["filename"].str[:-2]
        out_df = tmp_out_df.groupby(["filename"]).sum().reset_index()

    return out_df


def combine_into_one(cutadapt_df, trimmomatic_df, komplexity_df, decontam_df, is_single_end=False):
    out_df = cutadapt_df.merge(
        trimmomatic_df, on="filename", how="outer").merge(
        komplexity_df, on="filename", how="outer").merge(
        decontam_df, on="filename", how="outer")
    return out_df


def write_to_file(out_df, out_path):
    out_df.to_csv(out_path, sep="\t", index=False)


def summarize_qc_logs_main(args):

    # cutadapt
    cutadapt_paths = get_cutadapt_paths(args.qc_log_dir)
    cutadapt_df = cutadapt_summary(cutadapt_paths, args.is_single_end)

    # trimmomatic
    trimmomatic_paths = get_trimmomatic_paths(args.qc_log_dir)
    trimmomatic_df = trimmomatic_summary(trimmomatic_paths, args.is_single_end)

    # komplexity
    komplexity_paths = get_komplexity_paths(args.qc_log_dir)
    komplexity_df = komplexity_summary(komplexity_paths, args.is_single_end)

    # decontam
    decontam_paths = get_decontam_paths(args.qc_log_dir)
    decontam_df = decontam_summary(decontam_paths, args.is_single_end)

    # combine all results
    out_df = combine_into_one(cutadapt_df, trimmomatic_df, komplexity_df, decontam_df, args.is_single_end)

    # write output
    write_to_file(out_df, args.out_path)


def main():
    args = build_parser().parse_args(sys.argv[1:])
    summarize_qc_logs_main(args)


if __name__ == "__main__":
    main()