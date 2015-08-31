
/********************** NPLB **********************/
/**
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
**/

/**************************************************/

#ifndef _modelops_h
#define _modelops_h

typedef struct dataset{
  int **data;
  int features;
  int featureValues;
  int n;
  int tu;
}dataSet;

typedef struct modelstruct{
  int arch;
  int features;
  int featureValues;
  int n;
  int ***count;
  int *t;
  double lambda;
  int **pos;
  int *posCount;
  int **fvNoise;
  int *fnoise;
  float alpha;
  int tu;
}model;

typedef struct trainout{
  model *m;
  int *lp;
  double likelihood;
}trainOut;

model* createModel(int, dataSet*, int*, float, unsigned int*, double);
model* initializeModel(int, int, int, int);
void freeModel(model*);
void freeTo(trainOut*);
void freeToXModel(trainOut*);
void copyModel(model*, model*);
dataSet* getData(char*, char*);
void freeData(dataSet*);
int* sampleLabels(int, int, unsigned int*);
void nRandom(int*, int, int, unsigned int*);
void nRandom1(int*, int, int, unsigned int*);
void quickSort(int*, int, int);
void initPositions(model*, unsigned int*);
int sample(double*, int, double);
int maxPos(double*, int);
void getCount(model*, int*, dataSet*);
void getNoise(model*);
void saveCVResults(model*, char*, double, int, double);
void copyLabels(int*, int*, int);
#endif
