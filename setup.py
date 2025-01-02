from setuptools import setup, find_packages

VERSION = '0.0.3'
DESCRIPTION = 'A package that allows to derive force field.'

setup(
      version=VERSION,
      packages = find_packages(),
      package_data={'hessfit': ['json_files/*']},
      entry_points = {
        'console_scripts': ['hessfit=hessfit.hessfit:main',]
      },
     )
