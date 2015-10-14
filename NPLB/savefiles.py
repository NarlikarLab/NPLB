
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

import weblogoMod.weblogolib as wl
import plotExtras
import numpy
import os
import gc 
import sys
import pickle
from config import *
import copy

def createLogo(sequences, filename, pos, features, eps):    # Create logo using Weblogo 3.3
    if sequences == []: return
    seqs = wl.read_seq_data(sequences)
    data = wl.LogoData.from_seqs(seqs)
    options = wl.LogoOptions()
    options.title = ""
    options.size = "large"
    options.color_scheme = wl.colorscheme.monochrome
    options.stacks_per_line = 100
    options.text_font = "Arial-BoldMT"
    options.annotate = pos
    formt = wl.LogoFormat(data, options)
    fout = open(filename + ".png", "w")
    wl.png_formatter(data, formt, fout)
    fout.close()
    if eps == 1: 
        fout = open(filename + ".eps", "w")
        wl.eps_formatter(data, formt, fout)
        fout.close()

def sampleOne(l, n, tu, pos, features):    # Sample values based on model
    arr = ['A', 'C', 'G']
    arr = arr + ['T'] if tu == 0 else arr + ['U']
    l1 = map(list, zip(*map(lambda x: numpy.array(x).cumsum().searchsorted(numpy.random.sample(n)).tolist(), l)))
    l1 = map(lambda x: "".join(map(lambda y: arr[y], l1[x])), range(n))
    return l1

def makeImages(d, dirname, tu, tss, prefix, eps):    # Create logo for each architecture of given model
    lst = [map(lambda x: " ", range(1, d['features'] + 1)) for i in range(d['arch'])]
    numpy.random.seed(5)
    l = numpy.zeros(shape=(d['arch'], d['features'], d['featureValues']))
    for i in range(d['arch']):
        for j in range(d['features']):
            for k in range(d['featureValues']):
                l[i][j][k] = float(d['fvNoise'][j][k] + d['alpha'])/(d['fnoise'][j] + d['featureValues']*d['alpha'])
        for j in d['pos'][i]:
            for k in range(d['featureValues']):
                l[i][j][k] = float(d['count'][i][j][k] + d['alpha'])/(d['t'][i] + d['featureValues']*d['alpha'])
            lst[i][j] = "*"
        lst[i][tss] = "+1"
        if tss in d['pos'][i]: lst[i][tss] = "+1*"

        diffN = 25

        c = -diffN
        c1 = tss - diffN
        while c1 >= 0: 
            lst[i][c1] = str(c) + lst[i][c1]
            c = c - diffN
            c1 = c1 - diffN
        c = diffN
        c1 = tss + diffN - 1

        while c1 < d['features']:
            lst[i][c1] = str(c) + lst[i][c1]
            c = c + diffN
            c1 = c1 + diffN
    l = map(lambda x: sampleOne(l[x], d['t'][x], tu, d['pos'][x], d['features']), range(d['arch']))
    for i in range(d['arch']): 
        createLogo(l[i], dirname + "/" + prefix + str(i), lst[i], d['features'], eps) 


def makehtml(dirname, d, l):    # Create HTML file containing logos for best model learned by NPLB
    f = open(dirname + modelOutHTML, "w")
    f.write("<!DOCTYPE html>\n<html>\n<body>\n<h1>MODEL</h1>\n")
    f.write("<h3>Lambda: " + str(d['lambda']) + "</h3>\n")
    f.write("<h3>Dataset structure: " + str(d['n']) + " sequences with " + str(d['features']) + " features</h3>\n")
    f.write("<h3>Number of architectures in the best model: " + str(d['arch']) + "</h3>\n")
    f.write("<h3>Likelihood of best model: " + str(l) + "</h3>\n")
    for i in range(d['arch']):
        f.write("<h4>Architecture " + str(i+1) + ": " + str(d['t'][i]) + " sequences with " + str(d['posCount'][i]) + " important features</h4>\n")
        if d['t'][i] == 0:
            f.write("<h5>No Sequences</h5>\n")
        else:
            f.write("<h5>Sequence logo for the important positions in architecture " + str(i+1) + "</h5>\n")
            f.write("<img src=\"" + htmlFiles + "/" + str(i) + ".png\" style=\"border:thin solid black\">\n")
    f.write("<p><i>NOTE: All important positions in the logos are followed by an asterisk symbol and are coloured blue</i></p>")
    f.write("</body>\n</html>\n")
    f.close()


def makehtmlOrig(dirname, d, l, dO):    # Create HTML file containing logos for best model learned by NPLB along with the logo of raw data
    f = open(dirname + modelOutHTML, "w")
    f.write("<!DOCTYPE html>\n<html>\n<body>\n<h1>MODEL</h1>\n")
    f.write("<h3>Lambda: " + str(d['lambda']) + "</h3>\n")
    f.write("<h3>Dataset structure: " + str(d['n']) + " sequences with " + str(d['features']) + " features</h3>\n")
    f.write("<h3>Number of architectures in the best model: " + str(d['arch']) + "</h3>\n")
    f.write("<h3>Likelihood of best model: " + str(l) + "</h3>\n")
    for i in range(d['arch']):
        f.write("<h4>Architecture " + str(i+1) + ": " + str(d['t'][i]) + " sequences with " + str(d['posCount'][i]) + " important features</h4>\n")
        if d['t'][i] == 0:
            f.write("<h5>No Sequences</h5>\n")
        else:
            f.write("<h5>Sequence logo for the important positions in architecture " + str(i+1) + "</h5>\n")
            f.write("<img src=\"" + htmlFiles + "/" + str(i) + ".png\" style=\"border:thin solid black\">\n")
    f.write("<h5>Logo for the raw data</h5>\n")
    f.write("<img src=\"" + htmlFiles + "/" + rawDataImgPref + "0.png\" style=\"border:thin solid black\">\n")
    f.write("<p><i>NOTE: All important positions in the logos are followed by an asterisk symbol and are coloured blue</i></p>")
    f.write("</body>\n</html>\n")
    f.close()
    
def maketxt(dirname, d):    # Create text file containing details about the best model

    f = open(dirname + modelOutTxt, "w")
    f.write("MODEL\n\n")
    f.write("Lambda: " + str(d['m']['lambda']) + "\n\n")
    f.write("Dataset structure: " + str(d['m']['n']) + " sequences with " + str(d['m']['features']) + " features\n")
    f.write("Number of architectures in the best model: " + str(d['m']['arch']) + "\n\n")
    for i in range(d['m']['arch']):
        f.write("Architecture " + str(i+1) + ": " + str(d['m']['t'][i]) + " sequences with " + str(d['m']['posCount'][i]) + " important features\n")
        for j in range(d['m']['posCount'][i]):
            f.write(str(d['m']['pos'][i][j]+1) + " (");
            f.write(str(float(d['m']['count'][i][d['m']['pos'][i][j]][0] + d['m']['alpha'])/(d['m']['t'][i] + d['m']['featureValues']*d['m']['alpha'])) + " {" + str(d['m']['count'][i][d['m']['pos'][i][j]][0]) + "/" + str(d['m']['t'][i]) + "}")
            for k in range(1, d['m']['featureValues']):
                f.write(", " + str(float(d['m']['count'][i][d['m']['pos'][i][j]][k] + d['m']['alpha'])/(d['m']['t'][i] + d['m']['featureValues']*d['m']['alpha'])) + " {" + str(d['m']['count'][i][d['m']['pos'][i][j]][k]) + "/" + str(d['m']['t'][i]) + "}")
            f.write(")\n")
        f.write("\n")
    f.close()
    f = open(dirname + tempLabelsModelFile, "w")
    for i in d['lp']: f.write(str(i) + "\n")
    f.close()
    
    os.system("paste" + " " + dirname + tempLabelsModelFile + " " + dirname + tempLabelsFile + " " + ">" + " " + dirname + clusterDetailsFile)
    os.system("rm" + " " + dirname + tempLabelsModelFile  + " " + dirname + tempLabelsFile)

def makeImage(dirname, model, rfile, tss, imgfile, imgfileeps, inpfile):    # Create image matrix of input model
    os.system("cut" + " " + "-f1" + " " + dirname + inpfile + " " + ">" + " " + dirname + hiddenLabels)
    os.system("rev" + " " + dirname + inpfile + " " + "|" + " " + "cut" + " " + "-f1" + " | rev " + ">" + " " + dirname + hiddenOldData)
    
    indices = [[] for i in range(model['arch'])]
    j = 0
    with open(dirname + hiddenLabels) as infile:
        for line in infile:
            tmp = int(line)
            indices[model['arch'] - tmp] = indices[model['arch'] - tmp] + [j]
            j = j + 1

    f = open(dirname + hiddenData, "w")
    for i in indices:
        j = 0
        k = 0
        with open(dirname + hiddenOldData) as infile:
            for line in infile:
                try:
                    if k == i[j]:
                        f.write(line)
                        j = j + 1
                    k = k + 1
                except: pass
    f.close()
    if sys.platform == "darwin":
        os.system("sed" + " " + "-i" + " '' " + "'s/A/0\t/g;s/a/0\t/g;s/C/1\t/g;s/c/1\t/g;s/G/2\t/g;s/g/2\t/g;s/T/3\t/g;s/t/3\t/g;'" + " " + dirname + hiddenData)    # Modify input Fasta file to replace A, C, G, and T with 0, 1, 2 and 3 respectively on OS X.
    else:
        os.system("sed" + " " + "-i" + " " + "'s/A/0\t/g;s/a/0\t/g;s/C/1\t/g;s/c/1\t/g;s/G/2\t/g;s/g/2\t/g;s/T/3\t/g;s/t/3\t/g;'" + " " + dirname + hiddenData)    # Modify input Fasta file to replace A, C, G, and T with 0, 1, 2 and 3 respectively on Linux.
    f = open(dirname + hiddenDrawLines, "w")    # Save lines to be drawn on image matrix

    # Save labels for both axes of image matrix

    f1 = open(dirname + hiddenDrawLabels1, "w")
    c = 0
    for i in (model['t'][::-1])[:-1]:
        c = c + i
        f.write("-0.5\t" + str(c) + "\n" + str(model['features']-0.5) + "\t" + str(c) + "\n\n")
        f1.write(str(c) + "\n")
    f.close()
    f1.close()
    f = open(dirname + hiddenDrawLabels2, "w")
    c = 0
    for i in reversed(range(model['arch'])):
        f.write("A" + str(i + 1) + "\t" + str((c + c + model['t'][i])/2) + "\n\n")
        c = c + model['t'][i]
    f.close()
    lst = []
    gap = max(int(round(imgMatrixNumGap*model['features'])), 1)
    c = -gap
    c1 = tss - gap
    while c1 >= 0:
        lst = [(str(c1), str(c))] + lst
        c = c - gap
        c1 = c1 - gap

    lst = lst + [(str(tss), "+1")]
    c = gap
    c1 = tss + gap - 1
    while c1 < model['features']:
        lst = lst + [(str(c1), "+" + str(c))]
        c = c + gap
        c1 = c1 + gap
    f = open(dirname + hiddenDrawXTics, "w")
    for (i1, i2) in lst:
        f.write(i1 + "\t" + i2 + "\n")
    f.close()
    os.system("gnuplot" + " " + "-e" + " " + "'filename=\"" + dirname + hiddenData + "\"; var=\"" + dirname + imgfile + "\"; var1=\"" + dirname + hiddenDrawLines + "\"; var2=\"" + dirname + hiddenDrawLabels1 + "\"; var3=\"" + dirname + hiddenDrawLabels2 + "\"; var4=\"" + dirname + hiddenDrawXTics + "\"'" + " " + rfile[0] + " 2> /dev/null")
    if imgfileeps != "": os.system("gnuplot" + " " + "-e" + " " + "'filename=\"" + dirname + hiddenData + "\"; var=\"" + dirname + imgfileeps + "\"; var1=\"" + dirname + hiddenDrawLines + "\"; var2=\"" + dirname + hiddenDrawLabels1 + "\"; var3=\"" + dirname + hiddenDrawLabels2 + "\"; var4=\"" + dirname + hiddenDrawXTics + "\"'" + " " + rfile[1] + " 2> /dev/null")

    os.system("rm" + " " + "-f" + " " + dirname + "/.??*")


def savecvls(dirname, cvals):    # Save cross validation likelihood of the models learned
    if cvals == []: return
    maxArch = len(cvals)
    f = open(dirname + cvLikelihoods, "w")
    f.write("Cross validation likelihood of the best models\n\n")
    for i in range(maxArch):
        f.write(str(cvals[i][0]) + " architectures: ")
        if cvals[i][1] == 0: f.write("Not calculated\n")
        else: f.write(str(cvals[i][1]) + "\n")
    f.close()

def saveDetails(d, dirname, rfile, cvals, tss, flag, pEx, pCol, sBy, eps):
    dirname = dirname + "/"
    tmp_d_m_pos = d['m']['pos'][0]
    if ((tmp_d_m_pos[0] == 0 or tmp_d_m_pos[0] == 1) and (tmp_d_m_pos[1] == 0 or tmp_d_m_pos[1] == 1) and (tmp_d_m_pos[2] == 0 or tmp_d_m_pos[2] == 1)):
        for i in range(d['m']['arch']):
            d['m']['pos'][i] = filter(lambda x: d['m']['pos'][i][x] == 1, range(d['m']['features']))

    if flag == 0: pickle.dump(d['m'], open(dirname + bestModelFile, "wb"))
    os.system("rm" + " " + "-rf" + " " + dirname + htmlFiles)
    try:
        os.mkdir(dirname + htmlFiles)
    except OSError:
        print "ERROR: Cannot create directory", dirname + htmlFiles
        exit(2)

    if pEx != '':
        if pEx[0] != '' and pCol != 0 and sBy != 0 and flag == 0:
            cv = plotExtras.checkValid(pEx[0], sBy, d['m']['n'])
            if cv == -1 or cv == 2:
                print "Could not sort by values in column", sBy
                print "Please check -plotExtra file and/or -sortBy column number"
            else:
                d = plotExtras.rearrange(d, pEx, sBy)
    savecvls(dirname, cvals)
    makeImages(d['m'], dirname + htmlFiles, 0, tss, "", eps)
    if flag != 0: makehtml(dirname, d['m'], d['l'])
    if flag == 0: 
        # Save information about the raw data
        dOrig = {}
        dOrig['features'] = d['m']['features']
        dOrig['arch'] = 1
        dOrig['featureValues'] = d['m']['featureValues']
        dOrig['fvNoise'] = map(lambda z: map(lambda y: sum(map(lambda x: d['m']['count'][x][z][y], range(d['m']['arch']))), range(d['m']['featureValues'])), range(d['m']['features']))
        dOrig['pos'] = [[]]
        dOrig['alpha'] = d['m']['alpha']
        dOrig['fnoise'] = [d['m']['n'] for i in range(d['m']['features'])]
        dOrig['t'] = [d['m']['n']]
        makeImages(dOrig, dirname + htmlFiles, 0, tss, rawDataImgPref, eps)
        makehtmlOrig(dirname, d['m'], d['l'], dOrig)
        if rfile != 0: 
            os.system("sed -e 's/^/1\t/' " + dirname + tempLabelsFile + " > " + dirname + rawClusterDetailsFile)
            if eps == 0: makeImage(dirname, dOrig, rfile, tss, rawDataImage, "", rawClusterDetailsFile)
            else: makeImage(dirname, dOrig, rfile, tss, rawDataImage, rawDataImageEPS, rawClusterDetailsFile)
        maketxt(dirname, d)
    if rfile != 0 and flag == 0: 
        if eps == 0: makeImage(dirname, d['m'], rfile, tss, imageMatrix, "", clusterDetailsFile)
        else: makeImage(dirname, d['m'], rfile, tss, imageMatrix, imageMatrixEPS, clusterDetailsFile)
    if pEx != '':
        if pEx[0] != '' and pCol != 0 and flag == 0: plotExtras.plotExt(d, pEx, pCol, dirname)
    collected = gc.collect()
