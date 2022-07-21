import unittest
import sys
import os 
sys.path.append("/Users/kobiekirven/Desktop/scims/src")

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
        method_dict = build_chrom_coverage_dict(("chrom1","chrom2"), (35, 43),10)
        self.assertEqual(str(method_dict), str(test_dict))
        
        with self.assertRaises(ValueError):
            method_dict = build_chrom_coverage_dict(("chrom2"), (35, 43),10)

    
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


if __name__ == "__main__":
    unittest.main()
