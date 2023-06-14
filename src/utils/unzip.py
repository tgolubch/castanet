#!/usr/bin/env python3
import sys
from zipfile import PyZipFile

if __name__ == "__main__":
    for zip_file in sys.argv[1:]:
        pzf = PyZipFile(zip_file)
        pzf.extractall()