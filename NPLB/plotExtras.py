
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

import os
import numpy as np
import copy
from config import * 

def getNumbers(n):    # Calculate layout for pie charts of all architectures of best model
    for i in range(1,(n+1)):
        for j in range(i, (i+3)):
            if i*j >= n: return (i, j)

def checkValid(filename, colNum, total):    # Check format of input file. Input file must be a BED file.
    ln = 0
    n = 0
    if not os.path.isfile(filename):
        print "ERROR: Invalid file " + filename
        return -1
    try:
        with open(filename) as infile:
            for line in infile:
                if line[0:3] != 'chr' and line[0] == '\#': continue
                elif line[0:3] != 'chr': 
                    print "WARNING: Header must begin with #"
                    continue
                ln = ln + 1
                val = line.split()[colNum - 1]
                try:
                    tmp = float(val)
                    n = n + 1
                except ValueError: pass
        if total != ln:
            print "ERROR: Unequal number of sequences in file " + filename + " and input FASTA file"
            return -1
    except IndexError:
        print "ERROR: Invalid column " + str(colNum) + " in line: " + str(ln)
        return -1
    except:
        print "ERROR: Could not open file " + filename
        return -1
    if n == ln: return 1
    else: return 2

def piechart(arch, labels, filename, colNum, pname):
    outfile = pname + pieChartImage
    inpfile = pname + pieChartHiddenFile
    pfile = filename[2]
    filename = filename[0]
    (v3, v4) = getNumbers(arch + 1)
    dct = {}
    with open(filename) as infile:
        for line in infile:
            if line[0:3] != 'chr': continue
            val = line.split()[colNum - 1]
            if val not in dct:
                dct[val] = 0
    l = len(dct.keys())
    l1 = (" ").join(dct.keys())
    tmp = [(i, 0) for i in dct.keys()]
    segments = [dict(tmp) for _ in xrange(arch)]
    i = 0
    with open(filename) as infile:
        for line in infile:
            if line[0:3] != 'chr': continue
            val = line.split()[colNum - 1]
            segments[labels[i] - 1][val] = segments[labels[i] - 1][val] + 1
            i = i + 1
    total = [labels.count(i) for i in range(1, (arch + 1))]
    for i in range(1, (arch + 1)):
        f = open((inpfile + str(i)), "w")
        k = 0
        pv = 0
        for j in segments[i - 1]:
            f.write("0\t0\t1\t" + str(pv) + "\t")
            pv = min(int(round(float(segments[i - 1][j])/total[i - 1]*360)) + pv, 360)            # Save coordinates, radius and angle
            f.write(str(pv) + "\t" + str(k) + "\n")
            k = k + 1
        f.close()
    os.system("gnuplot -e 'var1=\"" + inpfile + "\"; var2=" + str(arch) + "; var3=" + str(v3) + "; var4=" + str(v4) + "; var5=\"" + outfile + "\"; var6=" + str(l) + "; var7=\"" + l1 + "\"' " + pfile)
    os.system("rm -rf " + inpfile + "*")

def boxplot(arch, labels, filename, colNum, dirname):
    outfile = dirname + boxPlotImage
    inpfile = dirname + boxPlotHiddenFile
    pfile = filename[1]
    filename = filename[0]

    lst = [[] for _ in xrange(arch)]
    i = 0
    with open(filename) as infile:
        for line in infile:
            if line[0:3] != 'chr': continue
            val = float(line.split()[colNum - 1])
            (lst[labels[i]-1]).append(val)
            i = i + 1
    i = arch + 1
    j = 0
    lst.reverse()
    f = open(inpfile, 'w')
    mn = []
    mx = []
    for l in lst:
        i = i - 1
        j = j + 1
        q1 = np.percentile(l, 25)
        q3 = np.percentile(l, 75)
        t1 = max(min(l), q1 - 1.5 * (q3 - q1))
        t2 = min(max(l), q3 + 1.5 * (q3 - q1))
        f.write(str(i) + " " + str(t1) + " " + str(q1) + " " + str(np.median(l)) + " " + str(q3) + " " + str(t2) + " PA" + str(i) + "\n")        # Save first quartile, median, third quartile and architecture count.
        mn.append(t1)
        mx.append(t2)
    f.close()
    mn1 = min(mn)
    mx1 = max(mx) 
    mx = mx1 + ((mx1-mn1)/10)
    mn = mn1 - ((mx1-mn1)/10)
    if mn == mx: 
        mx = mx + float(mx)/10
        mn = mn - float(mn)/10

    os.system("gnuplot" + " " + "-e" + " " + "'filename=\"" + inpfile + "\"; var1=\"" + outfile + "\"; var2=" + str(arch + 1) + "; var3=" + str(mn) + "; var4=" + str(mx) + "' " + pfile)
    os.system("convert -rotate 90 " + outfile + " " + outfile)

def rearrange(d, pEx, colNum):
    filename = pEx[0]
    lst = [[] for _ in xrange(d['m']['arch'])]
    i = 0
    with open(filename) as infile:
        for line in infile:
            if line[0:3] != 'chr': continue
            val = float(line.split()[colNum - 1])
            (lst[d['lp'][i]-1]).append(val)
            i = i + 1
    lst = map(lambda x: (np.median(lst[x]), x), range(d['m']['arch']))
    lst.sort()
    lst.reverse()
    (l1, l2) = zip(*lst)
    l2 = list(l2)
    l3 = map(lambda x: 0, range(d['m']['arch']))
    for i in range(d['m']['arch']): l3[l2[i]] = i
    for i in range(d['m']['n']): d['lp'][i] = l3[d['lp'][i] - 1] + 1
    tmp_pos = copy.deepcopy(d['m']['pos'])
    tmp_count = copy.deepcopy(d['m']['count'])
    tmp_t = copy.deepcopy(d['m']['t'])
    tmp_posCount = copy.deepcopy(d['m']['posCount'])
    for i in range(d['m']['arch']):
        d['m']['pos'][i] = tmp_pos[l2[i]]
        d['m']['count'][i] = tmp_count[l2[i]]
        d['m']['t'][i] = tmp_t[l2[i]]
        d['m']['posCount'][i] = tmp_posCount[l2[i]]
    return d

def plotExt(d, pEx, pCol, dirname):
    cv = checkValid(pEx[0], pCol, d['m']['n'])
    if cv == -1:
        print "Could not plot extras"
        return
    elif cv == 1: 
        boxplot(d['m']['arch'], d['lp'], pEx, pCol, dirname)
        print "\nBoxplot saved as " + dirname + boxPlotImage + "\n" 
    elif cv == 2: 
        piechart(d['m']['arch'], d['lp'], pEx, pCol, dirname)
        print "\nPie chart saved as " + dirname + pieChartImage + "\n"
