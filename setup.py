# -*- coding: utf-8 -*-
# from distutils.core import setup
# from Cython.Build import cythonize
# import numpy as np

from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension("huffman_tree",
                sources=["huffman_tree.pyx"],
                extra_compile_args=["-std=c++1y"],
                language="c++")

setup(name="RNG",
      ext_modules=cythonize(ext))


# from distutils.core import setup
# from Cython.Build import cythonize
#
# setup(
#     ext_modules=cythonize(
#         "huffman_tree.pyx",                 # Cython ソース
#         sources=["huffman_tree_.cpp"],  # その他のソースファイル
#         language="c++",             # C++ コードを生成させる
#     )
# )