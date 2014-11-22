#!/usr/bin/env python
from distutils.core import setup
setup(name="smallscript",
      py_modules=["smallscript", "OPA" ,"broaden", "functions_smsc","Dusch_unrest", "Btree"],
      author='Tobias Moehle ',
      version='0.3',
      description='program package to calculate vibronic spectra using g09-log files. 
		   Besides the common FC-picture, also several Duschinsky-rotated approaches are implemented.',
      author_email='tobias.moehle@uni-rostock.de',
      url='https://github.com/BalticPinguin/smallscript',
     )
