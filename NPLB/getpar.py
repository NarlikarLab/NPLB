
##################### NPLB #####################

#    No Promoter Left Behind (NPLB) is a tool to 
#    find the different promoter architectures within a set of promoter
#    sequences. More information can be found in the README file.
#    Copyright (C) 2015  Sneha Mitra and Leelavati Narlikar

#    NPLB is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    NPLB is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

################################################

import re
import sys
import os
from config import *

def printHelp1():    # Print usage of promoterLearn
    print("usage:\n./promoterLearn  [options]")
    print("\t-f fasta file with equal number of features in all sequences(compulsory)")
    print("\t-o output prefix")
    print("\t-a Pseudo count. Default 1")
    print("\t-t maximum number of iterations. Default: It keeps on iterating until the maximum likelihood stops increasing for n consecutive iterations where n is the number of sequences in the input file")
    print("\t-i 0 or 1. 1 to save image matrix and 0 to ignore. Default 1")
    print("\t-v 0 or 1. 1 to save likelihood plots of the learnt models. Default 0")
    print("\t-kfold k. Value of k in K-fold cross validation. Default 5")
    print("\t-lambda l. Non negative real number. Default: Varying lambda to find the best one")
    print("\t-lcount lc. Number of models to be learnt while training. The best model is considered. Default 5")
    print("\t-minarch m. Minimum number of architectures possible for the given data set. Default 1")
    print("\t-maxarch m. Maximum number of architectures possible for the given data set. Default 20")
    print("\t-proc maximum number of processors to use for computation. Default is the number of processors the system has")
    print("\t-plotExtra filename. Plots data from given file as pie charts or boxplots(depending on the type of data) for each architecture of the best model. Note: -pCol must also be set.")
    print("\t-pCol Natural number. Plots the data in the given column of the given file in the form of pie charts or boxplots with respect to the architectures of the best model. Note: -plotExtra must specify a valid filename") 
    print("\t-sortBy Natural number. Sort architectures in increasing order of median values calculated from values in column -sortBy of file -plotExtra. Note: -plotExtra must specify a valid filename and -pCol must be set")
    exit(2)

def printHelp2():    # Print usage of promoterClassify
    print("usage:\n./promoterClassify  [options]")
    print("\t-f fasta file with equal number of features in all sequences(compulsory)")
    print("\t-m file containing the model(compulsory)")
    print("\t-o output prefix")
    print("\t-i 0 or 1. 1 to save image matrix and 0 to ignore. Default 0. R needs to be installed for this option")
    print("\t-plotExtra filename. Plots data from given file as pie charts or boxplots(depending on the type of data) for each architecture of the best model. Note: -pCol must also be set")
    print("\t-pCol Natural number. Plots the data in the given column of the given file in the form of pie charts or boxplots with respect to the architectures of the best model. Note: -plotExtra must specify a valid filename") 
    print("\t-sortBy Natural number. Sort architectures in increasing order of median values calculated from values in column -sortBy of file -plotExtra. Note: -plotExtra must specify a valid filename and -pCol must be set")
    exit(2)

def printHelp(fv):
    if fv == 1: printHelp1()
    else: printHelp2()

def validNum(s, opt, fv):    # Check is string is a valid positive number
    s1 = re.search('^[0-9]+?(\.[0-9]+)\Z', s)
    if s1 is None: 
        print "ERROR: Invalid option: " +opt + " " + s + ". Must be a positive real number"
        printHelp(fv)
    s = s1
    s = float(s.group(0))
    if s == 0: 
        printHelp(fv)
    return s

def validNum1(s, opt, fv):    # Check is string is a valid non negative number
    s1 = re.search('^[0-9]+(\.[0-9]+)?\Z', s)
    if s1 is None: 
        print "ERROR: Invalid option: " +opt + " " + s + ". Must be a non negative real number"
        printHelp(fv)
    s = s1
    s = float(s.group(0))
    return s

def validFD(s, opt, fv):    # Check is string is a valid filename
    s1 = re.search('^.*\Z', s)
    if s1 is None: 
        print "ERROR: Invalid filename " + s
        printHelp(fv)
    s = s1
    return s.group(0)

def validInt(s, opt, fv):    # Check if string is a valid positive integer
    s1 = re.search('^[1-9][0-9]*\Z', s)
    if s1 is None: 
        print "ERROR: Invalid option: " + opt + " " + s + ". Must be positive integer"
        printHelp(fv)
    s = s1
    return int(s.group(0))

def validR(s, opt, fv):    # Check if string is 0 or 1
    s1 = re.search('^(0|1)\Z', s)
    if s1 is None: 
        print "ERROR: Invalid option: " + opt + " " + s + ". Must be 0 or 1"
        printHelp(fv)
    s = s1
    s = int(s.group(0))
    return s

def validDir(s):    # Check if string is a valid directory
    defVal = 0
    if s == "": 
        s = defaultDir
        defVal = 1
    elif s[-1] == "/":
        s = s[:-1]
    old = s
    d = os.path.dirname(s)
    if d != "" and not os.path.isdir(d):
        print("Invalid directory: " + d)
        exit(2)
    if os.path.isdir(s):
        b = 1
        while True:
            if os.path.isdir(s + "_" + str(int(b))): b = b + 1
            else:
                new = s + "_" + str(int(b))
                break
    else: new = s
    return [old, new, defVal]

def validFile(s):    # Check if string is a valid file
    if s == "": printHelp(1)
    if not os.path.isfile(s):
        print("Could not open file: " + s)
        exit(2)
    f = open(s, "r")
    lines = f.readlines()
    f.close
    return len(lines), lines[0]

def getValues():    # Read settings
    d = {'-a': defaultPseudoCount, '-o': '', '-f': '', '-t': 0, '-i': 1, '-kfold': defaultKFold, '-m': '', '-lcount': defaultLearnCount, '-proc': 0, '-maxarch': defaultMaxArch, '-minarch': defaultMinArch, '-lambda': -1, '-tss': 0, '-v': 0, '-plotExtra': '', '-pCol': 0, '-sortBy': 0}
    dF = {'-a': validNum, '-o': validFD, '-f': validFD, '-t': validInt, '-i': validR, '-kfold': validInt, '-m': validFD, '-lcount': validInt, '-proc': validInt, '-maxarch': validInt, '-minarch': validInt, '-lambda': validNum1, '-tss': validInt, '-v': validR, '-plotExtra': validFD, '-pCol': validInt, '-sortBy': validInt}
    lst = []
    fv = int((sys.argv)[1])
    sysargv = (sys.argv)[2:]
    if sysargv[:-1] == []: 
        printHelp(fv)
    for i in range(len(sysargv)):
        if i % 2 == 0: tmp = sysargv[i]
        else: lst = lst + [tmp + " " + sysargv[i]]
    for i in lst:
        s = re.search('^\-[a-zA-Z]* ', i)
        if s is None: 
            print "ERROR: Invalid option", i.split(" ")[0]
            printHelp(fv)
        s = s.group(0).strip()
        if fv == 0 and (s == '-kfold' or s == '-a' or s == '-t' or s == '-lcount' or s == '-proc' or s == '-maxarch' or s == '-minarch' or s == '-lambda' or s == '-v'):
            print "ERROR: Invalid option", s
            exit(2)
        if fv == 1 and (s == '-m'):
            print "ERROR: Invalid option", s
            exit(2)
        if s in d: d[s] = dF[s](i.split()[1], s, fv)
        else:
            print "ERROR: Invalid option", s
            printHelp(fv)
    if(d['-maxarch'] < d['-minarch']):
        print "ERROR: -maxarch cannot be less than -minarch"
        exit(2)
    validFile(d['-f'])
    if d['-m'] != '': validFile(d['-m'])
    if d['-i'] != 0: d['-i'] = sysargv[-1] + makeImageFile
    if d['-v'] != 0: d['-v'] = sysargv[-1] + makePlotFile
    if d['-plotExtra'] != '' and d['-pCol'] != 0: d['-plotExtra'] = (d['-plotExtra'], (sysargv[-1] + boxPlotFile), (sysargv[-1] + pieChartFile))
    if d['-plotExtra'] != '' and d['-pCol'] == 0:
        print "ERROR: -pCol must be set to a positive integer"
        printHelp(fv)
        exit(2)
    if d['-plotExtra'] == '' and d['-pCol'] != 0:
        print "ERROR: -plotExtra must be set to a valid filename"
        printHelp(fv)
        exit(2)
    if d['-sortBy'] != 0 and (d['-plotExtra'] == '' or d['-pCol'] == 0):
        print "ERROR: -plotExtra and -pCol must be set if -sortBy is non zero"
        printHelp(fv)
        exit(2)

    d['-o'] = validDir(d['-o'])
    try:
        os.mkdir(d['-o'][1])
    except:
        print("ERROR: Cannot create directory " + d['-o'][1])
        exit(2)
    return d, fv
