import unittest
import sys
import os 
sys.path.append("/Users/kobiekirven/Desktop/scims/SCIMS/src")

from determine_sex import *
from errors import WrongFile
import numpy as np

class TestDetermineSex(unittest.TestCase):

    def test_get_alignment_handle(self):
        type_file = "<class 'pysam.libcalignmentfile.AlignmentFile'>"
        self.assertEqual(str(get_alignment_handle("../tests/test_data/example.sam").__class__), type_file)
        self.assertEqual(str(get_alignment_handle("../tests/test_data/example.bam").__class__), type_file)
        with self.assertRaises(WrongFile):
            get_alignment_handle("../tests/test_data/not_a_bam.txt")

        with self.assertRaises(DoesNotExist):
            get_alignment_handle("../tests/test_data/not_a_bam2.txt")


    def test_chrom_len_from_sam(self):
        chrom_dict = {"NC_000913.3":4641652}
        sam_handle = get_alignment_handle("../tests/test_data/example.sam")
        bam_handle = get_alignment_handle("../tests/test_data/example.bam")
        self.assertEqual(chrom_len_from_sam(sam_handle), chrom_dict)
        self.assertEqual(chrom_len_from_sam(bam_handle), chrom_dict)

    def test_build_chrom_coverage_array(self):
        test_array = np.asarray([0,0,0,0,0])
        self.assertTrue((build_chrom_coverage_array(50, window_size=10) == test_array).all())
        self.assertTrue((build_chrom_coverage_array(53, window_size=10) == test_array).all())
        
        with self.assertRaises(ValueError):
            build_chrom_coverage_array(53, window_size="10")
            build_chrom_coverage_array(53, window_size=-10)
            build_chrom_coverage_array("4", window_size=10)
    
    def test_build_chrom_coverage_dict(self):
        test_dict = {"chrom1":np.asarray([0.,0.,0.]), "chrom2":np.asarray([0.,0.,0.,0.])}
        method_dict = build_chrom_coverage_dict({"chrom1":35,"chrom2":43},10)
        self.assertEqual(str(method_dict), str(test_dict))

    
    def test_add_to_coverage_dict(self):
        cov_array = np.asarray([0., 0., 0., 0., 0.])
        res = add_to_coverage_dict(cov_array, 0, 5, 10)
        print(res)
        self.assertTrue((res == np.asarray([5., 0., 0., 0., 0.])).all())
    
        res = add_to_coverage_dict(cov_array, 98, 11, 100)
        print(res)
        self.assertTrue((res == np.asarray([2., 9., 0., 0., 0.])).all())

        res = add_to_coverage_dict(cov_array, 299, 3, 100)
        print(res)
        self.assertTrue((res == np.asarray([0., 0., 1., 2., 0.])).all())

        res = add_to_coverage_dict(cov_array, 120, 82, 100)
        print(res)
        self.assertTrue((res == np.asarray([0., 80., 2., 0., 0.])).all())

        res = add_to_coverage_dict(cov_array, 320, 82, 100)
        print(res)
        self.assertTrue((res == np.asarray([0., 0., 0., 80., 2.])).all())

        with self.assertRaises(ValueError):
            add_to_coverage_dict(cov_array, 420, 82, 10)

    def test_decompose_sam_flag(self):
        self.assertEqual(decompose_sam_flag(25), ["PAIRED", "MUNMAP", "REVERSE"])
        self.assertEqual(decompose_sam_flag(144), ["REVERSE", "READ2"])

    
    def test_check_alignment(self):
        handle = pysam.AlignmentFile("../tests/test_data/one_sam_good.sam", "r", require_index=False)
        recs = [h for h in handle]
        align = check_alignment(recs[0])
        self.assertEqual(align, ("NC_000913.3", 191, 100))

        handle = pysam.AlignmentFile("../tests/test_data/one_bad_sam.sam", "r", require_index=False)
        recs = [h for h in handle]
        align = check_alignment(recs[0])
        self.assertEqual(align, None)

    def test_get_ks_stat(self):
        stat = get_ks_stat([1,5,3,4,5,3,3,4], [20,35,23,34,35,33,23,24])
        self.assertEqual(stat[0], 1.0)
        self.assertEqual(stat[1], 0.00015540015540015537)

    def test_get_chrom_windows_coverage(self):
        handle = pysam.AlignmentFile("../tests/test_data/one_sam_good.sam", "r", require_index=False)
        zeros = [0.]*4641
        test_dict  = build_chrom_coverage_dict({"NC_000913.3":4641652}, 1000)
        final = {"NC_000913.3":zeros}
        final["NC_000913.3"][0] = 100.
        self.assertEqual(list(get_chrom_windows_coverage(handle,test_dict, 1000)["NC_000913.3"]), final["NC_000913.3"])

        handle = pysam.AlignmentFile("../tests/test_data/sam_2.sam", "r", require_index=False)
        zeros = [0.]*4641
        test_dict  = build_chrom_coverage_dict({"NC_000913.3":4641652}, 1000)
        final = {"NC_000913.3":zeros}
        final["NC_000913.3"][1] = 100.
        final["NC_000913.3"][0] = 100.
        self.assertEqual(list(get_chrom_windows_coverage(handle,test_dict, 1000)["NC_000913.3"]), final["NC_000913.3"])

    def test_get_hom_het_lists(self):
        final = {"Chrom1":[0,4,5,6], "Chrom2":[7,8,9]}
        test = get_hom_het_lists(final, ["Chrom2"])
        self.assertEqual(list(test[0]), [0,4,5,6])
        self.assertEqual(list(test[1]), [7,8,9] )

        with self.assertRaises(ValueError):
            get_hom_het_lists(final, ["Not"])




if __name__ == "__main__":
    unittest.main()
