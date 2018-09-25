#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the `buildup` package."""

import numpy as np

import buildup.api as bu


def test_interval_div():
    assert bu._interval_div(1, 2, 3, 4) == (1 / 4, 2 / 3)
    assert bu._interval_div(3, 4, 1, 2) == (3 / 2, 4)


def test__element_wise_max():
    assert np.array_equal(bu._element_wise_max(np.asarray([1, 2, 3]), np.asarray([1, 1, 4])), np.asarray([1, 2, 4]))
    assert np.array_equal(bu._element_wise_max(np.asarray([[1, 2, 3]]), np.asarray([[1, 1, 4]])),
                          np.asarray([[1, 2, 4]]))
    assert np.array_equal(bu._element_wise_max(np.asarray([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
                                               np.asarray([[1, 1, 4], [5, 5, 5], [7, 7, 7]])),
                          np.asarray([[1, 2, 4], [5, 5, 6], [7, 8, 9]]))
