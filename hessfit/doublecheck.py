#!/usr/bin/env python3

import subprocess
import argparse

def remove_duplicates(file_path):
    seen = set()
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    with open(file_path, 'w') as file:
        for line in lines:
            if line not in seen:
                file.write(line)
                seen.add(line)

def main():
    parser = argparse.ArgumentParser(description="Remove duplicate lines from a file.")
    parser.add_argument("file_path", help="Path to the file to process")
    args = parser.parse_args()
    
    remove_duplicates(args.file_path)

if __name__ == "__main__":
    main()
