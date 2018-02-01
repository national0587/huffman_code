# -*- coding: utf-8 -*-
from distutils.core import setup
from Cython.Build import cythonize
import numpy as np

setup(
    name = "optimal_reduction",
    include_dirs=[np.get_include()],
    ext_modules = cythonize("optimal_reduction.pyx")
)