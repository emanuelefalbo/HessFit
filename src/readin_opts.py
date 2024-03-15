#!/usr/bin/env python3

import json
import os
import argparse

def commandline_parser1():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('optfile', nargs='?', help='option file in json')
    # parser.add_argument('-m', '--mode', choices=['modsem', 'ric'],
                        # default='ric', help='method to compute harmonic factors')
    return parser

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def commandline_parser2():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    # requiredNamed = parser.add_argument_group('mandatory arguments')
    parser.add_argument('optfile', nargs='?', help='option file in json')
    # requiredNamed.add_argument('-f1','--log_file', help='Gaussian QM log file ')
    # requiredNamed.add_argument('-f2','--fchk_file', help='Gaussain QM fchk file')
    parser.add_argument('-m', '--mode', choices=['all', 'mean'],
                        default='mean', help='averaging across same types; default = mean')
    txt = """path/to/amber.prm in Gaussain root directory
default = current directory """
    parser.add_argument('-path', nargs='?', help=txt, default=os.getcwd(), type=dir_path)
    return parser


def commandline_parser3():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('optfile', nargs='?', help='option file in json')
    parser.add_argument('--path', default='$g09root', help='path to Gaussian09/16 program')
    return parser

def read_optfile(fname):
    """
    Read the option file. expected in json style
    """
    if not os.path.exists(str(fname)):
        raise FileNotFoundError('Missing JSON file')
    with open(fname, 'r')  as fopen:
        data = json.load(fopen)
    for i in ["log_qm_file", "fchk_qm_file", "atype_file", "fchk_mm_file", "fchk_nb_file"]:
        if not os.path.exists(data['files'][i]):
            raise FileNotFoundError(f'Missing {i} file')
    return data

def read_optfile_2(fname):
    """
    Read the option file. expected in json style
    """
    if not os.path.exists(str(fname)):
        raise FileNotFoundError('Missing JSON file')
    with open(fname, 'r')  as fopen:
        data = json.load(fopen)
    for i in ["atom2type", "force_file", "file_xyz", "topol"]:
        if not os.path.exists(data['files'][i]):
            raise FileNotFoundError(f'Missing {i} file')
    return data

def read_optfile_3(fname):
    """
    Read the option file. expected in json style
    """
    if not os.path.exists(str(fname)):
        raise FileNotFoundError('Missing JSON file')
    with open(fname, 'r')  as fopen:
        data = json.load(fopen)
    for i in ["log_qm_file", "fchk_qm_file"]:
        if not os.path.exists(data['files'][i]):
            raise FileNotFoundError(f'Missing {i} file')
    return data