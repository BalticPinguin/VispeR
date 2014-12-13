#!/usr/bin/env python
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(name="smallscript",
      py_modules=["smallscript", "OPA" ,"broaden", "functions_smsc","Dusch_unrest", "Btree"],
      ext_modules=cythonize([Extension("Dusch_unrest", ["Dusch_unrest.py"]),
         Extension("Btree.py" ,["Btree.py"]),
         Extension("functions_smsc.py", ["functions_smsc.py"]),
         Extension("OPA.py", ["OPA.py"])]),
      author='Tobias Moehle ',
      version='0.4',
      description='program package to calculate vibronic spectra using g09-log files. Besides the common FC-picture, also several Duschinsky-rotated approaches are implemented.',
      author_email='tobias.moehle@uni-rostock.de',
      url='https://github.com/BalticPinguin/smallscript',
     )
