#!/usr/bin/env python3
#-*- coding: utf-8 -*-
## Created on 2024-01-26 at 16:09:23 CET by David Gaspard <gaspard.dave@gmail.com>
## Python module providing the compile_tikz() function to compile TikZ files.
## This module is used by many template scripts in the folder ./tikz/
## This file is also callable as a script. When called, it compiles the given TikZ files.
## USAGE : ./tikz/compile_tikz.py FILE_1.tikz [ FILE_2.tikz FILE_3.tikz ]
## FILENAME_TIKZ = Path of the TikZ file to be compiled by the present script.
import sys
import os
import shutil

LATEX_COMPILER = "pdflatex"    ## Compiler used to compile the LaTeX file.
GLOBAL_PREAMBLE_TEX = "tikz/preamble.tex"  ## LaTeX file containing the preamble (which does not hold in this Python script).

def compile_tikz(filename_tikz):
    """
    Compile the TikZ file "filename_tikz" in LaTeX using the global template file "tikz/preamble.tex".
    This function also checks for possible errors. Returns 1 on error, and 0 otherwise. 
    This function assumes that "pdflatex" and "grep" are available.
    """
    ## 1. First check for possible errors:
    if (shutil.which(LATEX_COMPILER) == None):
        print("[ERROR] LaTeX compiler not found: '" + LATEX_COMPILER + "'...")
        return 1
    if (not os.path.isfile(GLOBAL_PREAMBLE_TEX)):
        print("[ERROR] LaTeX preamble not found: '" + GLOBAL_PREAMBLE_TEX + "'...")
        return 1
    if (not os.path.isfile(filename_tikz)):
        print("[ERROR] TikZ file not found: '" + filename_tikz + "'...")
        return 1
    
    ## 2. Prepare the substitution dictionary:
    jobname = os.path.splitext(filename_tikz)[0]  ## The LaTeX jobname is simply the filename without extension.
    ##search_path = os.path.dirname(filename_tikz)  ## The search path used by PGFPlots to locate data files.
    
    dic = {## Substitution dictionary:
        "compiler": LATEX_COMPILER,
        "preamble": GLOBAL_PREAMBLE_TEX,
        "jobname": jobname,
        "filename_tikz": filename_tikz
        ##"search_path": search_path
    }
    
    ## \\pgfplotsset{table/search path={%(search_path)s}}
    
    ## 3. Prepare the UNIX commands to be executed:
    cmd = "%(compiler)s -jobname '%(jobname)s' '\\documentclass[12pt]{article}\\input{%(preamble)s}\\pagestyle{empty}\\begin{document}\\noindent\\input{\detokenize{%(filename_tikz)s}}\\end{document}' | grep -C 1 -wi --color=auto '^!\\|^l\\|error\\|undefined\\|unknown\\|missing\\|runaway\\|misplaced\\|multiply\\|exceeded\\|too\\|ended\\|extra\\|double\\|forget\\|forgotten\\|unbounded\\|overfull\\|underfull' " % dic
    
    ##print("COMMAND: " + cmd)
    
    ## 4. Execute the commands:
    print("[INFO] Compiling TikZ file: '" + filename_tikz + "'...")
    os.system(cmd)
    
    ## 5. Removes LaTeX's auxiliary files:
    flist = [jobname + ".aux", jobname + ".log", jobname + ".out"]
    for f in flist:
        if (os.path.isfile(f)):
            ##print("[INFO] Removing file: '" + f + "'...")
            os.remove(f)
    
    return 0

def answer_is_yes(msg):
    """
    Prompts the user with a binary choice.
    Returns True if the answer is "yes", False otherwise.
    """
    while(True):
        ans = input(msg)
        if (ans.lower() in ["y", "yes"]):
            return True
        elif (ans.lower() in ["n", "no"]):
            return False
        print("Please answer yes (y) or no (n)...")

def can_overwrite(filename):
    """
    Prompts the user if the file "filename" can be overwritten.
    Returns True if the file can be overwritten, False otherwise.
    """
    if (os.path.isfile(filename)):
        print("[WARN] File already exists: '" + filename + "'...")
        if (not answer_is_yes("Overwrite ? (y/n): ")):
            print("[INFO] OK, keeping file...")
            return False
    return True

def get_header_comment(filename, comment_char):
    """
    Returns the first commented lines of a file assuming the given comment character.
    """
    res = ""
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(comment_char):
                res += line
    return res

def main(args):
    """
    Main function of the script. This function is called when the script is prompted. 
    """
    if (len(args) == 1):
        print("[ERROR] No input file, doing nothing...")
        print("[USAGE] " + args[0] + " file_1.tikz [file_2.tikz file_3.tikz ...] ")
        return 1
    
    ## Compile all the given files:
    for f in args[1:]:
        compile_tikz(f)
    
    return 0

if (__name__ == '__main__'):
    exit(main(sys.argv))
