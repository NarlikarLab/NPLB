
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

from getpar import getValues
import gc
import os
import evaluate as ev
import savefiles as sf
import learndata as ld
import multiprocessing as mp
from config import *

def saveSettings(d):    # Save execution settings in file
    f = open(d['-o'][1] + "/" + settingsFile, "w")
    f.write("\nInput fasta file: " + os.path.abspath(d['-f']))
    f.write("\n\nOutput prefix: " + os.path.abspath(d['-o'][1]))
    f.write("\n\nMinimum number of architectures: " + str(d['-minarch']))
    f.write("\n\nMaximum number of architectures: " + str(d['-maxarch']))
    f.write("\n\nValue of k in K-Fold cross validation: " + str(d['-kfold']))
    f.write("\n\nNumber of models to be learnt per fold: " + str(d['-lcount']))
    f.write("\n\nPseudocount: " + str(d['-a']))
    f.write("\n\nNumber of iterations: ")
    if d['-t'] == 0: f.write("Default")
    else: f.write(str(d['-t']))
    f.write("\n\nVerbose mode: ")
    if d['-v'] == 0: f.write("No")
    else: f.write("Yes")
    if d['-lambda'] == -1: f.write("\n\nLambda: Varying")
    else: f.write("\n\nLambda: " + str(d['-lambda']))
    f.write("\n\nImage saved: ")
    if d['-i'] == 0: f.write("No")
    else: f.write("Yes")
    f.write("\n\nAdditional file: ")
    if d['-plotExtra'] != '': 
        f.write(d['-plotExtra'][0])
        f.write(" , column: " + str(d['-pCol']))
        if d['-sortBy'] != 0: f.write(" ,sort by column: " + str(d['-sortBy']))
    else: f.write("Not provided")
    f.write("\n\nMaximum number of cores used: " + str(d['-proc']))
    f.write("\n")
    f.close()

def printDetails(dirD, count):    # Print details of execution of promoterLearn
    print("\n\nProcessing details:\n")
    if dirD[2] != 0: print("Output directory was not provided. Output would saved in default directory " + dirD[1] + ".")
    else:
        if dirD[0] != dirD[1]: print("Directory " + dirD[0] + " was already present.")
        print("Output saved in " + dirD[1] + ".")
    print("Output file containing labels will be saved as " + dirD[1] + "/" + clusterDetailsFile + "\n" + "Model details will be saved as " + dirD[1] + "/" +  modelOutTxt + "\nModel details will also be saved in HTML format as " + dirD[1] + "/" + modelOutHTML + "\nCross validation likelihood of all the trained models will be saved as " + dirD[1] + "/" + cvLikelihoods)
    if d['-lambda'] == -1: print("An HTML file for the best model per lambda for each fold of every architecture will be saved as " + dirD[1] + "/Lambda_<lambda>/Architecture_<architecture>/Fold_<fold>/" + modelOutHTML)
    else: print("An HTML file for the best model for each fold of every architecture will be saved as " + dirD[1] + "/Lambda_<lambda>/Architecture_<architecture>/Fold_<fold>/" + modelOutHTML)
    print("Best model will be saved as " + dirD[1] + "/" + bestModelFile)
    if d['-lambda'] == -1: print("Best model will also be saved per lambda value as " + dirD[1] + "/Lambda_<lambda>/" + bestModelFile)
    print("Execution settings will be saved as " + dirD[1] + "/" + settingsFile)
    if d['-i'] != 0:
        print("Image matrices for the best model and the original unclustered one will be saved as " + dirD[1] + "/" + imageMatrix + " and " + dirD[1] + "/" + rawDataImage + " respectively.")
        if d['-lambda'] == -1: print("Image matrices will also be saved per lambda value in " + dirD[1] + "/Lambda_<lambda>/")
    if d['-v'] != 0:
        if d['-lambda'] != -1: print("Likelihood plot for the best model learned for each fold of every architecture would be saved as " + dirD[1] + "/Arch_<architecture>/Fold_<fold_number>/" + plotLikelihoodImage)
        else: print("Likelihood plot for the best model learned for each fold of every architecture and lambda would be saved as " + dirD[1] + "/Lambda_<lambda>/Arch_<architecture>/Fold_<fold_number>/" + plotLikelihoodImage)
    print("Maximum number of models to be learnt: " + str(count))
    print("Maximum number of cores that will be used: " + str(d['-proc']))
    print("\n\nProcessing. Please wait...\n")

def printLearnDetails(dirD):    # Print details of execution of promoterClassify
    print("\n\nProcessing details:")
    if dirD[2] != 0: print("Output directory was not provided. Output would be saved in default directory " + dirD[1] + ".")
    else:
        if dirD[0] != dirD[1]: print("Directory " + dirD[0] + " was already present.")
        print("Output saved in " + dirD[1] + ".")
    print("Output file containing labels will be saved as " + dirD[1] + "/" + clusterDetailsFile + " and model details will be saved as " + dirD[1] + "/" + modelOutTxt + ". The model details will also be saved in HTML format as " + dirD[1] + "/" + modelOutHTML + ". Best model will be saved as " + dirD[1] + "/" + bestModelFile)
    if d['-i'] != 0:
        print("Image matrix of the best model would be saved as " + dirD[1] + "/" + imageMatrix)
    print("\n\nProcessing. Please wait...")

def getFeatures(datafile):    # Calculate length of sequences in Fasta file
    features = 0
    odd = 0
    with open(datafile) as infile:
        for line in infile:
            if odd == 0: odd = 1
            elif line[0] != '>': 
                features = features + len(line.strip())
            else: break
    del odd
    collect = gc.collect()
    return features

def getModel(d):    # Learn best model
    dirname = d['-o'][1]
    n = mp.cpu_count()
    if d['-proc'] < n and d['-proc'] > 0: n = d['-proc']
    if 3*d['-kfold']*d['-lcount'] < n: n = 3*d['-kfold']*d['-lcount']
    d['-proc'] = n
    features = getFeatures(d['-f'])
    if d['-tss'] == 0: d['-tss'] = features/2
    if d['-tss'] > features:
        print "ERROR: -tss is more than the length of sequences"
        exit(1)
    count = d['-maxarch'] - d['-minarch'] + 1
    printDetails(d['-o'], count)
    saveSettings(d)
    m, cvals = ev.learn(d, dirname + "/" + tempLabelsFile, count)
    print "\n\nModel learnt successfully.\nSaving details..."
    sf.saveDetails(m, dirname, d['-i'], cvals, d['-tss'], 0, d['-plotExtra'], d['-pCol'], d['-sortBy'])
    print "Goodbye!"
    gc.collect()
    return dirname

def getLabels(d):    # Assign labels based on input model
    dirname = d['-o'][1]
    features = getFeatures(d['-f'])
    printLearnDetails(d['-o'])
    m = ld.learn(d['-f'], dirname + "/" + tempLabelsFile, d['-m'], features)
    print "\nModel learnt successfully.\nSaving details..."
    sf.saveDetails(m, dirname, d['-i'], [], d['-tss'], 0, d['-plotExtra'], d['-pCol'], d['-sortBy'])
    del features, m
    print "Goodbye!"
    gc.collect()
    return dirname
    
if __name__ == '__main__':
    try:
        d = None
        d, fv = getValues()
        if fv == 1: val = getModel(d)
        else: val = getLabels(d)
    except (KeyboardInterrupt, SystemExit):
        os.system('setterm -cursor on')
        if d is not None: 
#            if '-outFile' in d: d['-o'][1] = d['-outFile']
#            os.system ("rm"+ " " + "-rf" + " " + d['-o'][1])
#            print("\nDirectory " + d['-o'][1] + " deleted\n")
            print("\nExiting...\n")
        exit(2)
    exit(0)

