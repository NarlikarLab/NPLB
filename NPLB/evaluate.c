
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

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include<signal.h>
#include "messages.h"
#include "modelops.h"
#include "traindata.h"
#include "evaluate.h"

void sig_handler(int signo){
  if(signo == SIGINT)
    printf("\n");
    exit(1);
}

/* Randomized sequence of indices */
int* posList(n){
  int *pos;
  unsigned int seed;
  seed = 1;
  pos = (int*)malloc(sizeof(int)*n);
  if(!pos) printMessages(0, NULL);
  nRandom1(pos, n, n, &seed);
  return pos;
}

/* Get subset for training of given fold */
dataSet* getTrainSubset(dataSet *ds, int subset, int kfold, int *pos){
  int startPos, endPos, i, j, tr;
  dataSet *trainSet;
  startPos = subset*ds->n/kfold;
  endPos = subset == (kfold-1) ? ds->n : (subset+1)*ds->n/kfold;
  trainSet = (dataSet*)malloc(sizeof(dataSet));
  if(!trainSet) printMessages(0, NULL);
  trainSet->n = ds->n - endPos + startPos;
  trainSet->tu = ds->tu;
  trainSet->data = (int**)malloc(sizeof(int*)*trainSet->n);
  if(!(trainSet->data)) printMessages(0, NULL);
  trainSet->features = ds->features;
  trainSet->featureValues = ds->featureValues;
  tr = 0;
  for(i = 0; i < startPos; i++){
    (trainSet->data)[tr] = (int*)malloc(sizeof(int)*trainSet->features);
    if(!(trainSet->data)[tr]) printMessages(0, NULL);
    for(j = 0; j < trainSet->features; j++)
      (trainSet->data)[tr][j] = (ds->data)[pos[i]][j];
    tr++;
  }
  for(i = endPos; i < ds->n; i++){
    (trainSet->data)[tr] = (int*)malloc(sizeof(int)*trainSet->features);
    if(!(trainSet->data)[tr]) printMessages(0, NULL);
    for(j = 0; j < trainSet->features; j++)
      (trainSet->data)[tr][j] = (ds->data)[pos[i]][j];
    tr++;
  }
  return trainSet;
}

/* Get subset for training of given fold */
dataSet* getTestSubset(dataSet *ds, int subset, int kfold, int *pos){
  int startPos, endPos, i, j, te;
  dataSet *testSet;
  startPos = subset*ds->n/kfold;
  endPos = subset == (kfold-1) ? ds->n : (subset+1)*ds->n/kfold;
  testSet = (dataSet*)malloc(sizeof(dataSet));
  if(!testSet) printMessages(0, NULL);
  testSet->n = endPos - startPos;
  testSet->tu = ds->tu;
  testSet->data = (int**)malloc(sizeof(int*)*testSet->n);
  if(!(testSet->data)) printMessages(0, NULL);
  testSet->features = ds->features;
  testSet->featureValues = ds->featureValues;
  te = 0;
  for(i = startPos; i < endPos; i++){
    (testSet->data)[te] = (int*)malloc(sizeof(int)*testSet->features);
    if(!(testSet->data)[te]) printMessages(0, NULL);
    for(j = 0; j < testSet->features; j++){
      (testSet->data)[te][j] = (ds->data)[pos[i]][j];
    }
     te++;
  }
  return testSet;
}
