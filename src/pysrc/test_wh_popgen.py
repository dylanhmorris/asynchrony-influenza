#!/usr/bin/env python3

import wh_popgen as whp
import unittest
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal

class TestFM(unittest.TestCase):
        
    def test_no_emerge(self):
        """
        check that freq is 
        0 whenever t_final < t_emerge
        """
        x1 = whp.f_m(
            t_final = 5,
            t_emerge = 10,
            delta = 5,
            f_0 = 0.52,
            t_M = 0,
            t_peak = 1,
            R0 = 5,
            d_v = 4,
            mu = 0.33e-5)
        
        x2 = whp.f_m(
            t_final = 20,
            t_emerge = 23,
            delta = 0,
            f_0 = 0.28,
            t_M = 0.5,
            t_peak = 200,
            R0 = 5,
            d_v = 4,
            mu = 0.33e-5)
        
        self.assertAlmostEqual(x1, 0)
        self.assertAlmostEqual(x2, 0)

    def test_bounds_check(self):
        with self.assertRaises(ValueError):
            whp.freq_mut_ongoing(2, 1.1, 5, 0, 0)
        with self.assertRaises(ValueError):
            whp.freq_mut_ongoing(2, 0, 5, 0, 0)
        with self.assertRaises(ValueError):
            whp.freq_mut_ongoing(2, -1, 5, 0, 0)
        pass


    def test_t_star(self):
        """
        Check that the t-stars 
        are good approximations
        """
        f_target = 0.9
        t = 3
        R0 = 5
        d_v = 4
        bottleneck = 1
        t_M = 0.01525
        k = 10
        mu = .33e-5
        c_w = 1

        for c_m, tol in zip([0, 0.25, 0.5, 0.75],
                            [4, 4, 3, 2]):
            delta = k * (c_w - c_m)

            self.assertAlmostEqual(
                whp.invert_f_m(
                    f_target,
                    t,
                    t_M,
                    R0,
                    d_v,
                    bottleneck,
                    k,
                    c_w,
                    c_m,
                    ongoing_mutate = False,
                    mu = mu),
                whp.t_star(f_target,
                           t,
                           t_M,
                           R0,
                           d_v,
                           bottleneck,
                           k,
                           c_w,
                           c_m),
                places = tol)




class TestPTransmit(unittest.TestCase):

    def test_no_selection(self):
        """
        Check that all no selection situations
        are equivalent, as they should be
        """
        k_a, k_b, k_c = 0, 1000, 0.25
        x1 = whp.p_transmit(k_a, 0, 0, .33e-5, 4, 5, 1,
                            2, 2, 0.99, 0.98, 10, 2)
        x2 = whp.p_transmit(k_b, 0, 0, .33e-5, 4, 5, 1,
                            2, 2, 0.99, 0.98, 10, 2)
        self.assertAlmostEqual(x1, x2, places = 10)

        x3 = whp.p_transmit(k_b, 0.5, 0.5, .33e-5, 4, 5, 1,
                            5, 5, 0.25, 0.999, 1000, 2)
        x4 = whp.p_transmit(k_c, 0.5, 0.5, .33e-5, 4, 5, 1,
                            5, 5, 0.25, 0.999, 1000, 2)
        self.assertAlmostEqual(x3, x4, places = 10)

            

if __name__ == "__main__":
    unittest.main()
