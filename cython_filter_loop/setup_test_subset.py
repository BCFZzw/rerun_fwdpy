from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np
import sys
sys.path.insert(1, '/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/cython_filter_loop')

extensions = [
    Extension(
        name="test_subset",
        sources=["test_subset.pyx"],
        define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_25_API_VERSION')]
    )
]

setup(
    ext_modules=cythonize("test_subset.pyx"),
    include_dirs=[np.get_include()]   
)