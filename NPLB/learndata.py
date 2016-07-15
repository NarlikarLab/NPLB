
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

import multiprocessing as mp
from cstructures import *
from config import *
import pickle
import os

# Assign labels to data set based on scores calculated by input model

def learn(datafile, outfile, modelfile, features, tempfile):
    ds = libctest.getData(datafile, outfile)
    try:
        with open(modelfile, "rb") as fp:
            d = pickle.load(fp)
    except:
        print "ERROR: Could not open file", modelfile
        exit(1)

    for i in range(d['arch']):
        tmp = [0 for k in range(d['features'])]
        for j in d['pos'][i]:
            tmp[j] = 1
        d['pos'][i] = tmp

    if(features != d['features']):
        print "ERROR: Sequences length unequal for data file", datafile, "and model", modelfile
        print "Sequence length in datafile:", features
        print "Sequence length in model:", d['features']
        exit(1)
    m = createModel(d)
    to = libctest.callLearnData(ds, byref(m), tempfile);
    tout = getTrainOut(to)
    libctest.freeToXModel(to)
    libctest.freeData(ds)
    del m, d
    return tout
