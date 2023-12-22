#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module for internal utility functions.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import ctypes
import logging
import multiprocessing as mp

import numpy as np

logger = logging.getLogger(__name__)

__all__ = [
    "as_ndarray",
    "as_dtype",
    "as_float32",
    "as_int32",
    "as_uint8",
    "as_uint16",
    "as_c_float_p",
    "as_c_bool_p",
    "as_c_uint8_p",
    "as_c_uint16_p",
    "as_c_int",
    "as_c_int_p",
    "as_c_float",
    "as_c_char_p",
    "as_c_void_p",
    "as_c_size_t",
]


def as_ndarray(arr, dtype=None, copy=False):
    if not isinstance(arr, np.ndarray):
        arr = np.array(arr, dtype=dtype, copy=copy)
    return arr


def as_dtype(arr, dtype, copy=False):
    if not arr.dtype == dtype:
        arr = np.array(arr, dtype=dtype, copy=copy)
    return arr


def as_float32(arr):
    arr = as_ndarray(arr, np.float32)
    return as_dtype(arr, np.float32)


def as_int32(arr):
    arr = as_ndarray(arr, np.int32)
    return as_dtype(arr, np.int32)


def as_uint16(arr):
    arr = as_ndarray(arr, np.uint16)
    return as_dtype(arr, np.uint16)


def as_uint8(arr):
    arr = as_ndarray(arr, np.uint8)
    return as_dtype(arr, np.uint8)


def as_c_float_p(arr):
    c_float_p = ctypes.POINTER(ctypes.c_float)
    return arr.ctypes.data_as(c_float_p)


def as_c_bool_p(arr):
    c_bool_p = ctypes.POINTER(ctypes.c_bool)
    return arr.ctypes.data_as(c_bool_p)


def as_c_uint8_p(arr):
    c_uint8_p = ctypes.POINTER(ctypes.c_uint8)
    return arr.ctypes.data_as(c_uint8_p)


def as_c_uint16_p(arr):
    c_uint16_p = ctypes.POINTER(ctypes.c_uint16)
    return arr.ctypes.data_as(c_uint16_p)


def as_c_int(arr):
    return ctypes.c_int(arr)


def as_c_long(arr):
    return ctypes.c_long(arr)


def as_c_int_p(arr):
    arr = arr.astype(np.intc, copy=False)
    c_int_p = ctypes.POINTER(ctypes.c_int)
    return arr.ctypes.data_as(c_int_p)


def as_c_float(arr):
    return ctypes.c_float(arr)


def as_c_char_p(arr):
    return ctypes.c_char_p(arr.encode())


def as_c_void_p():
    return ctypes.POINTER(ctypes.c_void_p)


def as_c_size_t(arr):
    return ctypes.c_size_t(arr)


def to_numpy_array(obj, dtype, shape):
    return np.frombuffer(obj, dtype=dtype).reshape(shape)


def is_contiguous(arr):
    return arr.flags.c_contiguous
