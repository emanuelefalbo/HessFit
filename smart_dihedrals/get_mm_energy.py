#!/usr/bin/env python3

import re
import argparse

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
    return par

if __name__ == '__main__':
    PAR = build_parser()
    OPTS = PAR.parse_args()
    file_list = OPTS.filelist

    match = " Energy="
    file_list.sort(key=natural_keys)
    for file in file_list:
        # print(file)
        with open(file, 'r') as fname:
             for line in fname:
                 if line[:8] == match:
                    print(line[11:28]) 


