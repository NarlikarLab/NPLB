
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

#ifndef _traindata_h
#define _traindata_h

typedef struct traindatastruct{
  model *m;
  dataSet *ds;
  int *lp;
  int times;
  double likelihood;
  unsigned int seed;
  char *filename;
}trainDataStruct;

trainOut* callTrainData(dataSet*, int, unsigned int, int, float, double, char*);
void sig_handler1(int);
double calculateLikelihood(model*);
double calculateLikelihoodLabel(model*, double, int, int, int*);
void calculateLikelihoodPos(model*, double*, int, int);
void addRemoveDataPoint(model*, dataSet*, int*, int, int);
int sampleNewLabel(model*, dataSet*, int, unsigned int*, int, int);
void *trainData(void *vtds);
trainOut* callLearnData(dataSet*, model*, char*);
int sampleImpFeature(model*, int, unsigned int*, int, double*, trainDataStruct*, double*, int*, long long*, FILE*);
int solvePosExp(model*, int, int, unsigned int*, int);
double cvLikelihood(dataSet*, model*);
#endif
