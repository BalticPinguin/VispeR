#!/usr/bin/env python
from distutils.core import setup
setup(name="smallscript",
      py_modules=["smallscript", "OPA" ,"broaden", "functions_smsc"],
      author='Tobias Moehle',
      version='0.2',
      description='program package to calculate vibronic spectra using g09-log files.',
      author_email='tobias.moehle@uni-rostock.de',
      url='https://github.com/BalticPinguin/smallscript',
     )
