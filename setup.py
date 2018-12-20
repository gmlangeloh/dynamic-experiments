from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

modules = \
[
    "buchberger.pyx",
    "unrestricted.pyx",
    "caboara_perry.pyx",
    "polynomials.pyx",
    "heuristics.pyx",
    "stats.pyx"
]

my_include_dirs = \
[
    "$SAGE_ROOT/local/include/singular"
]

my_libraries = \
[
    "m",
    "readline",
    "Singular",
    "givaro",
    "gmpxx",
    "gmp"
]

extensions = \
[
    Extension(
        "Dynamic", modules,
        include_dirs = my_include_dirs,
        libraries = my_libraries
    )
]

setup(
    name = "Dynamic Buchberger",
    ext_modules = cythonize(extensions)
)
