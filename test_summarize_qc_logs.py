import os
import unittest
import numpy as np
import pandas as pd
import summarize_qc_logs

ASSETS_DIR = "test_summarize_qc_logs/"


class TestSummarizeQcLogs(unittest.TestCase):

    def test_get_cutadapt_paths(self):

        expected_paths = [os.path.join(ASSETS_DIR, "cutadapt", "sample_1.log"),
                          os.path.join(ASSETS_DIR, "cutadapt", "sample_2.log"),
                          os.path.join(ASSETS_DIR, "cutadapt", "sample_3.log")]
        out_paths = summarize_qc_logs.get_cutadapt_paths(ASSETS_DIR)
        self.assertCountEqual(expected_paths, out_paths)

        with self.assertRaises(AssertionError) as e:
            summarize_qc_logs.get_cutadapt_paths("bad_path")
        self.assertIn("No cutadapt files found", str(e.exception))

    def test_cutadapt_paired(self):
        in_paths = [os.path.join(ASSETS_DIR, "cutadapt", "sample_1.log"),
                          os.path.join(ASSETS_DIR, "cutadapt", "sample_2.log"),
                          os.path.join(ASSETS_DIR, "cutadapt", "sample_3.log")]
        expected_df = pd.DataFrame([["sample_1", 4521113, 4521113, 1356333900, 1236192329],
                                    ["sample_2", 7301623, 7301623, 2190486900, 2076703029],
                                    ["sample_3", 456123, 555555, 200000, 180000]],
                                   columns=["filename",
                                            "cutadapt_read_pairs_processed",
                                            "cutadapt_read_pairs_written",
                                            "cutadapt_bps_processed",
                                            "cutadapt_bps_written"])
        out_df = summarize_qc_logs.cutadapt_summary(in_paths, is_single_end=False)
        pd.testing.assert_frame_equal(expected_df, out_df)

    def test_get_trimmomatic_paths(self):

        expected_paths = [os.path.join(ASSETS_DIR, "trimmomatic", "sample_1.out"),
                          os.path.join(ASSETS_DIR, "trimmomatic", "sample_2.out"),
                          os.path.join(ASSETS_DIR, "trimmomatic", "sample_3.out")]
        out_paths = summarize_qc_logs.get_trimmomatic_paths(ASSETS_DIR)
        self.assertCountEqual(expected_paths, out_paths)

        with self.assertRaises(AssertionError) as e:
            summarize_qc_logs.get_trimmomatic_paths("bad_path")
        self.assertIn("No trimmomatic files found", str(e.exception))


    def test_trimmomatic_paired(self):
        in_paths = [os.path.join(ASSETS_DIR, "trimmomatic", "sample_1.out"),
                    os.path.join(ASSETS_DIR, "trimmomatic", "sample_2.out"),
                    os.path.join(ASSETS_DIR, "trimmomatic", "sample_3.out")]
        expected_df = pd.DataFrame([["sample_1", 4521113, 4309073],
                                    ["sample_2", 7301623, 7121205],
                                    ["sample_3", 100000, 50]],
                                   columns=["filename",
                                            "trimmomatic_read_pairs_processed",
                                            "trimmomatic_read_pairs_written"])
        out_df = summarize_qc_logs.trimmomatic_summary(in_paths, is_single_end=False)
        pd.testing.assert_frame_equal(expected_df, out_df)


    def test_get_komplexity_paths(self):

        expected_paths = [os.path.join(ASSETS_DIR, "komplexity", "sample_1.filtered_ids"),
                          os.path.join(ASSETS_DIR, "komplexity", "sample_2.filtered_ids"),
                          os.path.join(ASSETS_DIR, "komplexity", "sample_3.filtered_ids")]
        out_paths = summarize_qc_logs.get_komplexity_paths(ASSETS_DIR)
        self.assertCountEqual(expected_paths, out_paths)

        with self.assertRaises(AssertionError) as e:
            summarize_qc_logs.get_komplexity_paths("bad_path")
        self.assertIn("No komplexity files found", str(e.exception))


    def test_komplexity_paired(self):
        in_paths = [os.path.join(ASSETS_DIR, "komplexity", "sample_1.filtered_ids"),
                    os.path.join(ASSETS_DIR, "komplexity", "sample_2.filtered_ids"),
                    os.path.join(ASSETS_DIR, "komplexity", "sample_3.filtered_ids")]
        expected_df = pd.DataFrame([["sample_1", 3],
                                    ["sample_2", 20],
                                    ["sample_3", 0]],
                                   columns=["filename",
                                            "komplexity_read_pairs_removed"])
        out_df = summarize_qc_logs.komplexity_summary(in_paths, is_single_end=False)
        pd.testing.assert_frame_equal(expected_df, out_df)

    def test_get_decontam_paths(self):
        expected_paths = [os.path.join(ASSETS_DIR, "decontam", "sample_1_1.txt"),
                          os.path.join(ASSETS_DIR, "decontam", "sample_2_1.txt"),
                          os.path.join(ASSETS_DIR, "decontam", "sample_3_1.txt"),
                          os.path.join(ASSETS_DIR, "decontam", "sample_1_2.txt"),
                          os.path.join(ASSETS_DIR, "decontam", "sample_2_2.txt"),
                          os.path.join(ASSETS_DIR, "decontam", "sample_3_2.txt")]
        out_paths = summarize_qc_logs.get_decontam_paths(ASSETS_DIR)
        self.assertCountEqual(expected_paths, out_paths)

        with self.assertRaises(AssertionError) as e:
            summarize_qc_logs.get_decontam_paths("bad_path")
        self.assertIn("No decontam files found", str(e.exception))

    def test_decontam_paired(self):
        in_paths = [os.path.join(ASSETS_DIR, "decontam", "sample_1_1.txt"),
                    os.path.join(ASSETS_DIR, "decontam", "sample_2_1.txt"),
                    os.path.join(ASSETS_DIR, "decontam", "sample_3_1.txt"),
                    os.path.join(ASSETS_DIR, "decontam", "sample_1_2.txt"),
                    os.path.join(ASSETS_DIR, "decontam", "sample_2_2.txt"),
                    os.path.join(ASSETS_DIR, "decontam", "sample_3_2.txt")
                    ]
        expected_df = pd.DataFrame([["sample_1", 80000, 32000],
                                    ["sample_2", 310, 205],
                                    ["sample_3", 3027, 31]],
                                   columns=["filename",
                                            "decontam_read_pairs_written",
                                            "decontam_read_pairs_removed"])
        out_df = summarize_qc_logs.decontam_summary(in_paths, is_single_end=False)
        pd.testing.assert_frame_equal(expected_df, out_df)

    def test_combine_into_one_paired(self):
        cutadapt_df = pd.DataFrame([["sample_1", 1, 2, 3, 4],
                                    ["sample_2", 5, 7, 8, 6]],
                                   columns=["filename",
                                            "cutadapt_read_pairs_processed",
                                            "cutadapt_read_pairs_written",
                                            "cutadapt_bps_processed",
                                            "cutadapt_bps_written"])
        trimmomatic_df = pd.DataFrame([["sample_1", 8, 10]],
                                   columns=["filename",
                                            "trimmomatic_read_pairs_processed",
                                            "trimmomatic_read_pairs_written"])
        komplexity_df = pd.DataFrame([["sample_1", 3], ["sample_2", 20]],
                                     columns=["filename", "komplexity_read_pairs_removed"])
        decontam_df = pd.DataFrame([["sample_1", 20, 30],
                                    ["sample_2", 33, 31]],
                                   columns=["filename",
                                            "decontam_read_pairs_written",
                                            "decontam_read_pairs_removed"])

        expected_df = pd.DataFrame([["sample_1", 1, 2, 3, 4, 8, 10, 3, 20, 30],
                                    ["sample_2", 5, 7, 8, 6, np.nan, np.nan, 20, 33, 31]],
                                   columns=["filename",
                                            "cutadapt_read_pairs_processed", "cutadapt_read_pairs_written",
                                            "cutadapt_bps_processed", "cutadapt_bps_written",
                                            "trimmomatic_read_pairs_processed", "trimmomatic_read_pairs_written",
                                            "komplexity_read_pairs_removed",
                                            "decontam_read_pairs_written", "decontam_read_pairs_removed"])

        out_df = summarize_qc_logs.combine_into_one(cutadapt_df, trimmomatic_df, komplexity_df, decontam_df, is_single_end=False)
        pd.testing.assert_frame_equal(expected_df, out_df)

    def test_write_to_file(self):
        out_df = pd.DataFrame([[1, 2], [3, 4]], columns=["filename", "cutadapt_something"])
        out_path = "test_write_qc_summary.txt"

        # file should not already exist
        self.assertTrue(not(os.path.exists(out_path)), f"{out_path} should not yet exist")

        # write file
        summarize_qc_logs.write_to_file(out_df, out_path)

        # make sure it now exists
        self.assertTrue(os.path.exists(out_path))

        # clean up
        os.remove(out_path)

    def test_summarize_qc_logs_main_paired(self):
        out_path = "test_summarize_qc_logs_main_paired.txt"

        # file should not already exist
        self.assertTrue(not (os.path.exists(out_path)), f"{out_path} should not yet exist")

        args_string = f"-i {ASSETS_DIR} -o {out_path}"
        args = summarize_qc_logs.build_parser().parse_args(args_string.split())
        summarize_qc_logs.summarize_qc_logs_main(args)

        # make sure output exists
        self.assertTrue(os.path.exists(out_path))

        # clean up
        os.remove(out_path)


if __name__ == "__main__":
    unittest.main()