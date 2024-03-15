#!/usr/bin/env python3

import re
import argparse
import pandas as pd

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def build_parser():
    par = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    txt = 'MM log files'
    par.add_argument('filelist', type=str, help=txt, nargs='+')
    par.add_argument('-o', type=str, help='name output file', default='output.csv', nargs='?')
    par.add_argument('-t', type=str, default='mm', help='type of log file', nargs='?')
    return par

if __name__ == '__main__':
    PAR = build_parser()
    OPTS = PAR.parse_args()
    file_list = OPTS.filelist

    match_mm = " Energy="
    match_qm = " SCF Done:"
    file_list.sort(key=natural_keys)
    res = []
    for file in file_list:
        # print(file)
        with open(file, 'r') as fname:
             for line in fname:
                 if OPTS.t == 'mm':
                    if line[:8] == match_mm:
                       res.append(line[11:28]) 
                 elif OPTS.t == 'qm':
                     if line[:10] == match_qm:
                       res.append(line[23:40]) 
    df =pd.DataFrame(res)
    fout = OPTS.o
    df.to_csv(fout)
                       


