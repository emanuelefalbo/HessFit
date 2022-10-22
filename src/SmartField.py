#!/usr/bin/env python3

import subprocess
import os

scripts_to_run = ['build_4Smart.py'] #, 'Smart_harmonic.py']

for s in scripts_to_run:
    subprocess.call([os.path.join(os.getcwd(), s)])