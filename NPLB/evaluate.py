
##################### NPLB #####################

#    No Promoter Left Behind (NPLB) is a tool to 
#    find the different promoter architectures within a set of promoter
#    sequences. More information can be found in the README file.
#    Copyright (C) 2015  Sneha Mita and Leelavati Narlikar

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

import multiprocessing as mp
import multiprocessing.sharedctypes as mpc
from cstructures import *
from config import *
import pickle
import os
import savefiles as sf
import numpy as np
import time
import sys

trainSets = []
testSets = []
ds = 0
v = mp.Value('i', 0)
lockL = mp.Lock()
lcount = 0
d = []
toArr = []
cvals = []

# Train data during final model learning phase

def singleEval1(args):
    try:
        sds = args[3].value
        f = args[1]
        lock = args[0]
        args = [ds] + args[2:] 
        tout = libctest.callTrainData(*args)
        to = getTrainOut(tout)
        libctest.freeTo(tout)
        lockL.acquire()
        with lock:
            t = pickle.dump(to, open(f, "wb"))
            v.value = 1
    except (KeyboardInterrupt, SystemExit):
        exit(1)

# Train data during cross validation phase

def singleEval(args):
    try:
        sds = args[4].value
        k = args[2]
        f = args[1]
        lock = args[0]
        args = [trainSets[k-1]] + args[3:] 
        tout = libctest.callTrainData(*args)
        to = getTrainOut(tout)
        libctest.freeTo(tout)
        lockL.acquire()
        with lock:
            t = pickle.dump([[k, to, sds]], open(f, "wb"))
            v.value = 1
    except (KeyboardInterrupt, SystemExit):
        exit(1)

def nameTmpFile(a, i, s):
    if d['-v'] == 0: return "0"
    return (d['-o'][1] + "/." + str(a) + "_" + str(i) + "_" + str(s))

# Perform cross validation in a parallel manner

def multiEval(picklefile, count, archC):
    global v, best
    pos = [[] for i in range(d['-maxarch'])]
    countsArr = [0 for i in range(d['-maxarch'])]     # Array to count how many model learned with given number of architectures.
    # Array to check for occupancy of processors.

    procArr = [0 for i in range(d['-proc'])]    # Array to check for occupancy of processors. This is just to enumerate the processes in execution.

    checkList = np.zeros(3)
    chooseFrom = [[-1 for i in range(d['-kfold'])] for i in range(d['-maxarch'])]    # Array to check which model is chosen as best for given fold and architecture.

    cArr = [[0 for i in range(d['-kfold'])] for i in range(d['-maxarch'])]    # Count number of models learned for given fold and architecture

    rs = sum([range(1, (d['-lcount']+1)) for i in range(d['-kfold'])], [])   # Array of seeds for each fold
    rk = range(d['-lcount']*d['-kfold'])    # Array of number of models to learn per architecture

    # Cross validation likelihood results should ideally increase,
    # reach a maximum and decrease again. 
    # Here it attempts to find the maximum in a binary search fashion,
    # in parallel execution.
 
    low = d['-minarch'] - 1
    high = d['-maxarch'] - 1
    while low <= high:
        p = [-5 for i in range(d['-proc'])]
        mid = int((low+high)/2)
        if low == high: 
            if best[0] == 0: best = (low+1, 0)
            break
        trs = []
        tes = []
        seeds = []
        archs = []
        idx = 0
        c = 0
        pc = 0
        lock = mp.Lock()

        # Check values at mid-1, mid and mid + 1 positions and
        # count number of models to learn.

        if mid > d['-minarch'] - 1 and not(cvals[mid - 1]):
            seeds = seeds + rs
            archs = archs + [mid for i in rk]
            c = c + d['-lcount']*d['-kfold']
        if not(cvals[mid]):
            seeds = seeds + rs
            archs = archs + [mid+1 for i in rk]
            c = c + d['-lcount']*d['-kfold']
        if mid != d['-maxarch']-1 and not(cvals[mid+1]):
            seeds = seeds + rs
            archs = archs + [mid+2 for i in rk]
            c = c + d['-lcount']*d['-kfold']

        # Set variable to the first among the three (or less)
        # positions which are to be learned.

        checkIndex = -1
        if mid > d['-minarch'] - 1 and cvals[mid-1] == 0: checkIndex = mid-1
        elif cvals[mid] == 0: checkIndex = mid

        # Set the maximum number of models that can be
        # learned in parallel.

        if c < d['-proc']: pc = c
        else: pc = d['-proc']

        # Check which all models are to be learned.

        if checkIndex > -1: 
            r = map(lambda x: str(x), range(checkIndex + 1, (checkIndex + (c/(d['-lcount']*d['-kfold'])) + 1)))
            if r[0] != r[-1]: 
                r = r[:-1] + [" and ", r[-1]]
                if r[1] != " and ": r = [r[0], ", "] + r[1:]
            print "**** Learning models with", "".join(r), "architectures ****\n"    # Print a message about the models to be learned.
        flag = 0
        try:
            for i in range(pc):
                pr = mp.Process(target = singleEval, args = ([lock, picklefile, idx + 1, c_int(d['-t']), c_uint(seeds[i]), c_int(archs[i]), c_float(d['-a']), c_double(d['-lambda']), nameTmpFile(archs[i], idx + 1, seeds[i])],))    # Start the maximum number of processes that can be learned in parallel.

                if seeds[i] == d['-lcount']: idx = (idx + 1)%d['-kfold']        # Switch to next fold if all models of current fold are processing/learned.

                p[i] = pr
                p[i].start()        # Start process
                procArr[i] = archs[i]-1        # Set the architecture which is being learned.
            i = pc
            k = 0

            joined = 0
            while joined < c:            # Join processes which are complete and start the pending ones in a loop till all are over.
                time.sleep(0.5)
                for j in range(pc):
                    if procArr[j] != -1:
                        if not p[j].is_alive():    # Check if complete
                            p[j].join()
                            countsArr[procArr[j]] = countsArr[procArr[j]] + 1    # Increment the number of models learned for given architecture.
                            procArr[j] = -1
                            joined = joined + 1
                            if i < c:
                                pr = mp.Process(target = singleEval, args = ([lock, picklefile, idx + 1, c_int(d['-t']), c_uint(seeds[i]), c_int(archs[i]), c_float(d['-a']), c_double(d['-lambda']), nameTmpFile(archs[i], idx + 1, seeds[i])],))    # Start pending operations
                                if seeds[i] == d['-lcount']: idx = (idx + 1)%d['-kfold']
                                p[j] = pr
                                p[j].start()
                                procArr[j] = archs[i]-1
                                i = i + 1
                if os.path.isfile(picklefile):
                    with lock:
                        with open(picklefile, "rb") as f:    # Read output from binary file and save the best model per fold and architecture.
                            tmpL = pickle.load(f)
                            for tmp in tmpL:
                                if toArr[tmp[1]['m']['arch'] - 1][tmp[0] - 1] == 0 or tmp[1]['l'] > toArr[tmp[1]['m']['arch'] - 1][tmp[0] - 1]['l']:    # Update count and selected model accordingly.
                                    chooseFrom[tmp[1]['m']['arch'] - 1][tmp[0] - 1] = tmp[2]
                                    toArr[tmp[1]['m']['arch']-1][tmp[0] - 1] = tmp[1]
                                cArr[tmp[1]['m']['arch'] - 1][tmp[0] - 1] = cArr[tmp[1]['m']['arch'] - 1][tmp[0] - 1] + 1
                                if cArr[tmp[1]['m']['arch'] - 1][tmp[0] - 1] == d['-lcount']:    # Print message if all models of given fold and architecture have been learned.
                                    print "> Completed training fold", tmp[0], "of model with", tmp[1]['m']['arch'], "architectures"
                        if v.value == 1:
                            lockL.release()
                            v.value = 0
                tcheckIndex = checkIndex
                cnt = 0
                if checkIndex == mid - 1: cnt = 3
                elif checkIndex == mid: cnt = 2
                elif checkIndex == mid + 1: cnt = 1
                
                for ic in range(cnt):
                    if (checkIndex+ic) < d['-maxarch'] and countsArr[checkIndex+ic] == d['-kfold']*d['-lcount'] and cvals[checkIndex+ic] == 0:    # Compute cross validation likelihood if models of all folds of given architecture have been learned.


                        cvals[checkIndex+ic] = sum([libctest.cvLikelihood(testSets[x], byref(createModel(toArr[checkIndex+ic][x]['m']))) for x in range(d['-kfold'])])/d['-kfold']
                        saveInfo(checkIndex + ic)

                        if best[0] == 0 or cvals[checkIndex+ic] > best[1]:
                            best = (checkIndex+ic+1, cvals[checkIndex+ic])
                        print "\n>>Cross validation likelihood computed for model with", (checkIndex + ic + 1), "architectures\n"
                        print "Best model till now has", best[0], "architectures\n"
                        archC = archC + 1
                        print "PROGRESS:", str(int(round(float(archC)/count*100))) + "% complete\n"
                if checkIndex < d['-maxarch'] and countsArr[checkIndex] == d['-kfold']*d['-lcount']:        # Search for highest cross validation likelihood and update low and high accordingly.
                    if checkIndex == mid - 1: checkIndex = checkIndex + 1
                    if checkIndex == mid and mid > d['-minarch'] - 1 and cvals[checkIndex] != 0:
                        checkIndex = checkIndex + 1
                        if cvals[mid-1] <= cvals[mid]:
                            if mid == d['-maxarch']-1:
                                low = mid
                                high = mid
                                flag = 1
                                break
                        else:
                            high = mid - 1
                            flag = 1
                            break
                    elif checkIndex == mid and cvals[mid] != 0: checkIndex = checkIndex + 1
                    if checkIndex == mid+1 and cvals[checkIndex] != 0:
                        if cvals[mid] != 0 and cvals[mid] > cvals[mid+1] and ((mid > d['-minarch'] - 1 and cvals[mid-1] < cvals[mid]) or (mid <= d['-minarch'] - 1)):
                            low = mid
                            high = mid
                            flag = 1
                            break
                        else:
                            low = mid + 1
                            flag = 1
                            break
                if flag == 1:
                    break
            if flag == 1:
                for i in p:
                    if type(i) is int: continue
                    i.terminate()
                    i.join()
                pickle.dump([], open(picklefile, "wb"))
                continue
            if checkIndex < d['-maxarch'] and countsArr[checkIndex] == d['-kfold']*d['-lcount']:
                if checkIndex == mid - 1: checkIndex = checkIndex + 1
                if checkIndex == mid and mid > d['-minarch'] - 1 and cvals[checkIndex] != 0:
                    checkIndex = checkIndex + 1
                    if cvals[mid-1] <= cvals[mid]:
                        if mid == d['-maxarch']-1:
                            low = mid
                            high = mid
                            flag = 1
                            continue
                    else:
                        high = mid - 1
                        flag = 1
                        continue
                elif checkIndex == mid: checkIndex = checkIndex + 1
                if checkIndex == mid+1 and cvals[checkIndex] != 0:
                    if cvals[mid] > cvals[mid+1] and ((mid > d['-minarch'] - 1 and cvals[mid-1] < cvals[mid]) or (mid <= d['-minarch'] - 1)):
                        low = mid
                        high = mid
                        flag = 1
                        continue
                    else:
                        low = mid + 1
                        flag = 1
                        continue
        except (KeyboardInterrupt, SystemExit): 
           for i in p: 
                if type(i) is int: continue
                i.terminate()
                i.join()
           exit(1)

        if checkIndex == -1:    # If all models are learned just search for the best.
            if mid > 0 and cvals[mid-1] <= cvals[mid]:
                if mid == (d['-maxarch'] - 1):
                    low = mid
                    high = mid
                elif cvals[mid] > cvals[mid + 1]:
                    low = mid
                    high = mid
                else: low = mid + 1
            else: high = mid - 1

    return count, archC

# Learn best model after cross validation

def learnModel(arch, picklefile):
    pc = min(d['-proc'], d['-lcount'])
    seeds = range(1, (d['-lcount'] + 1))
    p = [-5 for i in range(d['-proc'])]
    procArr = [-1 for i in range(d['-proc'])]
    maxL = 0
    lock = mp.Lock()
    pickle.dump([], open(picklefile, "wb"))
    try:
        for i in range(pc):   # Learn in parallel
            pr = mp.Process(target = singleEval1, args = ([lock, picklefile, c_int(d['-t']), c_uint(seeds[i]), c_int(arch), c_float(d['-a']), c_double(d['-lambda']), "0"],))
            p[i] = pr
            p[i].start()
            procArr[i] = seeds[i]
        i = pc
        k = 0
        joined = 0
        while joined < d['-lcount']:
            time.sleep(1)
            for j in range(pc):
                if procArr[j] != -1:
                    if not p[j].is_alive():
                        p[j].join()
                        joined = joined + 1
                        procArr[j] = -1
                        progLearn(float(joined)/d['-lcount'])
                        if i < d['-lcount']:
                            pr = mp.Process(target = singleEval1, args = ([lock, picklefile, c_int(d['-t']), c_uint(seeds[i]), c_int(arch), c_float(d['-a']), c_double(d['-lambda']), "0"],))
                            p[j] = pr
                            p[j].start()
                            procArr[j] = seeds[i]
                            i = i + 1
            if os.path.isfile(picklefile):
                with lock:
                    with open(picklefile, "rb") as f:
                        tmpL = pickle.load(f)
                        if tmpL != []:        # Selct model with highest likelihood.
                            if maxL == 0 or tmpL['l'] > maxL['l']: 
                                maxL = tmpL
                        if v.value == 1:
                            lockL.release()
                            v.value = 0
    except (KeyboardInterrupt, SystemExit):
        for i in p:
            if type(i) is int: continue
            i.terminate()
            i.join()
        exit(1)

    return maxL


def getData(ds1):    # Copy data to shared memory for parallel processing
    ds = mpc.copy(ds1)
    return ds

def saveInfo(c):
    tDir = d['-o'][1] + "/" + "Arch_" + str(c + 1)
    try:
        # Save HTML file for best model of each fold an architecture

        os.mkdir(tDir)
        for fold in range(d['-kfold']):
            td = tDir + "/Fold_" + str(fold + 1)
            os.mkdir(td)
            sf.saveDetails(toArr[c][fold], td, d['-i'], [], d['-tss'], 1, d['-plotExtra'], d['-pCol'], d['-sortBy'])
    except OSError:
        print "ERROR: Cannot create directory in", d['-o'][1]
        exit(2)
    if c == 0 or d['-v'] == 0: return
    for fold in range(d['-kfold']):

        # Save likelihood plot if -v flag is set

        for sd in range(d['-lcount']):
            if os.path.isfile(d['-o'][1] + "/." + str(c + 1) + "_" + str(fold + 1) + "_" + str(sd)):
                os.system("mv" + " " + d['-o'][1] + "/." + str(c + 1) + "_" + str(fold + 1) + "_"  + str(sd) + " " + tDir + "/Fold_" + str(fold + 1) + "/.plot")
        os.system("gnuplot" + " " + "-e" + " " + "'filename=\"" + tDir + "/Fold_" + str(fold + 1) + "/" + plotFileHiddenName + "\"; var=\"" + tDir + "/Fold_" + str(fold + 1)  + "/" + plotLikelihoodImage + "\"'" + " " + d['-v'])
        os.system("rm" + " " + "-f" + " " + tDir + "/Fold_" + str(fold + 1) + "/" + plotFileHiddenName)

def learnOne():

    # Learn model with one architecture if -minarch is set to 1. This model is then learned at the beginning
    # since no parallel processing is required.

    global best
    print "**** Learning models with 1 architecture ****\n"
    for i in range(d['-kfold']):
        tout = libctest.callTrainData(trainSets[i], c_int(d['-t']), c_uint(1), c_uint(1), c_float(d['-a']), c_double(d['-lambda']), "0")
        toArr[0][i] = getTrainOut(tout)
        libctest.freeTo(tout)
        print "> Completed training fold", (i+1), "of model with 1 architecture"
    cvals[0] = sum([libctest.cvLikelihood(testSets[x], byref(createModel(toArr[0][x]['m']))) for x in range(d['-kfold'])])/d['-kfold']
    print "\n>>Cross validation likelihood computed for model with 1 architecture\n"
    print "Best model till now has 1 architecture\n"
    saveInfo(0)
    best = (1, cvals[0])

def evalF(picklefile, count, archC):

    # Call function to perform cross validation and then learn the best model

    if d['-minarch'] == 1:
        learnOne()
        archC = 1
        print "PROGRESS:", str(int(round(float(archC)/count*100))) + "% complete\n"
    pickle.dump([], open(picklefile, "wb"))
    while 1:
        count, archC = multiEval(picklefile, count, archC)
        flag = 0
        if best[1] == 0: break
        if flag == 0: break
        d['-maxarch'] = best[0]
        mid = (d['-minarch'] - 1 + d['-maxarch'] - 1)/2
        flag = 0
        if mid >= d['-minarch'] and toArr[mid-1][0] == 0: flag = flag + 1
        if toArr[mid][0] == 0: flag = flag + 1
        if (mid + 1) < d['-maxarch'] and toArr[mid + 1][0] == 0: flag = flag + 1
        if flag == 0: break
        print "Some of the learned models have tiny architectures.\nLearning simpler models.\n\nProcessing. Please wait...\n"
        break
    print "\nPROGRESS: 100% complete\n"
    printBestModel(best[0])
    progLearn(0)
    return learnModel(best[0], picklefile)

# Learn models for given lambda

def learnDiffLambda(dt, outfile, count, ds, trainSets, testSets):
    global d, toArr, cvals, best
    d = dt
    picklefile = d['-o'][1] + "/" + tempFile
    toArr = [[0 for i in range(d['-kfold'])] for i in range(d['-maxarch'])]
    cvals = np.zeros(d['-maxarch'])    # Array to save cross validation likelihood
    best = (0, 0)
    archC = 0
    tmax = d['-maxarch']
    m = evalF(picklefile, count, archC)
    os.remove(picklefile)    # Remove temporary binary file created while execution
    os.system("rm" + " " + "-f" + " " + str(d['-o'][1]) + "/.??*")    # Remove temporary hidden files created while execution
    cvals = zip(range(1, tmax + 1), cvals)[(d['-minarch']-1):]
    return m, cvals

def learn(dt, outfile, count):
    global d, ds, trainSets, testSets, lcount
    lcount = dt['-lcount']
    ds = getData(libctest.getData(dt['-f'], outfile))
    pos = libctest.posList(ds.contents.n)

    # Learn best model directly if -minarch and -maxarch are same

    if dt['-maxarch'] == dt['-minarch']:
        d = dt
        printBestModel(dt['-maxarch'])
        m = learnModel(dt['-maxarch'], dt['-o'][1] + "/" + tempFile)
        return m, []
    for i in range(dt['-kfold']):    # Get randomized train sets and test sets for every fold
        trainSets.append(getData(libctest.getTrainSubset(ds, i, dt['-kfold'], pos)))
        testSets.append(getData(libctest.getTestSubset(ds, i, dt['-kfold'], pos)))


    if dt['-lambda'] != -1:    # Learn models by single lambda value
        m, cvals = learnDiffLambda(dt, outfile, count, ds, trainSets, testSets)
    else:    # Learn models by varying lambda
        print "\n\nTrying Lambda", 0, "\n\n"
        finalOut = dt['-o'][1]
        dt['-lambda'] = 0
        dt['-outFile'] = finalOut
        dt['-o'][1] = finalOut + "/" + defaultLambdaFile + str(0)
        try:
            os.mkdir(dt['-o'][1])
        except:
            print("ERROR: Cannot create directory " + d['-o'][1])
            exit(2)
        os.system("cp " + finalOut + "/" + tempLabelsFile + " " + dt['-o'][1] + "/")
        
        m, cvals = learnDiffLambda(dt, outfile, count, ds, trainSets, testSets)
        sf.saveDetails(m, dt['-o'][1] + "/", d['-i'], cvals, dt['-tss'], 0, dt['-plotExtra'], dt['-pCol'], dt['-sortBy'])
        bestCVL = best[1]
        i = 1
        while(1):
            print "\n\nTrying Lambda", i, "\n\n"
            dt['-lambda'] = i
            dt['-o'][1] = finalOut + "/" + defaultLambdaFile + str(i)
            try:
                os.mkdir(dt['-o'][1])
            except:
                print("ERROR: Cannot create directory " + d['-o'][1])
                exit(2)
            os.system("cp " + finalOut + "/" + tempLabelsFile + " " + dt['-o'][1] + "/")
            m1, cvals1 = learnDiffLambda(dt, outfile, count, ds, trainSets, testSets)
            sf.saveDetails(m1, dt['-o'][1] + "/", d['-i'], cvals1, dt['-tss'], 0, dt['-plotExtra'], dt['-pCol'], dt['-sortBy'])
            posMin = min(m1['m']['posCount'])
            if posMin == 0 or best[1] < bestCVL or (best[1] == bestCVL and best[0] == m['m']['arch']):    # Exit when minimum number of important features is 0 or when best corss validation likelihood is lesser compared to the one by previous lambda or when cross validation likelihood and number of architecture remains same
                del m1, cvals1
                break
            bestCVL = best[1]
            del m, cvals
            m = m1
            cvals = cvals1
            if posMin < 5: break
            i = i + 2
    for i in range(dt['-kfold']):
        libctest.freeData(trainSets[i])
        libctest.freeData(testSets[i])
    libctest.freeData(ds)
    del trainSets, testSets
    return m, cvals

def printBestModel(arch):
    print "\n\nStructure of best model\n"
    print "Number of architectures:", arch
    print "\n\nLearning best model. Please wait...\n"

def progLearn(s):
    hashes = '#' * int(round(s * 40))
    spaces = ' ' * (40 - len(hashes))
    sys.stdout.write("\r[{0}] {1}%".format(hashes + spaces, int(round(s*100))))
    sys.stdout.flush()
