#!/usr/bin/env python
from distutils.core import setup
from distutils.extension import Extension

setup(name="smallscript",
      py_modules=["smallscript", "functions_smsc", "broaden"],
      author='Tobias Moehle ',
      version='0.7',
      description='program package to calculate vibronic spectra using g09-log files. Besides the common FC-picture, also several Duschinsky-rotated approaches are implemented.',
      author_email='tobias.moehle@uni-rostock.de',
      url='https://github.com/BalticPinguin/smallscript',
     )
