from distutils.core import setup, Extension
from distutils.sysconfig import *
from distutils.util import *
import os
import os.path
import numpy
from distutils.command.build_py import build_py
from Cython.Distutils import build_ext
from Cython.Build import cythonize

math_lib = ['m']

extra_compile_args = ['-g', '-Wall', '-O2']

py_inc = [get_python_inc()]

np_lib = os.path.dirname(numpy.__file__)
np_inc = [os.path.join(np_lib, 'core/include')]
cmdclass = {'build_py': build_py}

cmdclass.update({'build_ext': build_ext})
packages=['test']

extensions = [        Extension("dtw",  ["src/dtw/dtw.pyx", "src/dtw/cdtw.c"],
                      extra_compile_args=extra_compile_args,
                      libraries=math_lib,
                      include_dirs=py_inc + np_inc,
                      language = 'c')
            ]

setup(name = 'test',
      version='0.0.2',
      requires=['numpy (>=1.3.0)'],
      description='abea',
      cmdclass=cmdclass,
      ext_modules=cythonize(extensions),
      )