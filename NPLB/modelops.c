
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
#include<string.h>
#include<math.h>
#include "messages.h"
#include "modelops.h"

/* Create model  */
model* createModel(int arch, dataSet *ds, int *lp, float alpha, unsigned int *seed, double lambda){
  model *m;
  m = initializeModel(arch, ds->features, ds->featureValues, ds->n);
  m->n = ds->n;
  m->arch = arch;
  m->lambda = lambda;
  m->features = ds->features;
  m->featureValues = ds->featureValues;
  m->alpha = alpha;
  getCount(m, lp, ds);
  initPositions(m, seed);
  getNoise(m);
  return m;
}

/* Allocate space to model structure */
model* initializeModel(int arch, int features, int featureValues, int n){
  model *m;
  int i, j, k;
  m = (model*)malloc(sizeof(model));
  if(!m) printMessages(0, NULL);
  m->pos = (int**)malloc(sizeof(int*)*arch);
  if(!m->pos) printMessages(0, NULL);
  m->count = (int***)malloc(sizeof(int**)*arch);
  if(!m->count) printMessages(0, NULL);
  m->t = (int*)malloc(sizeof(int)*arch);
  if(!m->t) printMessages(0, NULL);
  m->fvNoise = (int**)malloc(sizeof(int*)*features);
  if(!m->fvNoise) printMessages(0, NULL);
  m->fnoise = (int*)malloc(sizeof(int)*features);
  if(!m->fnoise) printMessages(0, NULL);
  m->posCount = (int*)malloc(sizeof(int)*arch);
  if(!m->posCount) printMessages(0, NULL);
  for(i = 0; i < arch; i++){
    (m->pos)[i] = (int*)malloc(sizeof(int)*features);
    if(!(m->pos)[i]) printMessages(0, NULL);
    (m->count)[i] = (int**)malloc(sizeof(int*)*features);
    if(!(m->count)[i]) printMessages(0, NULL);
    for(j = 0; j < features; j++){
      (m->pos)[i][j] = 0;
      (m->count)[i][j] = (int*)malloc(sizeof(int)*featureValues);
      if(!(m->count)[i][j]) printMessages(0, NULL);
      for(k = 0; k < featureValues; k++) (m->count)[i][j][k] = 0;
    }
    (m->t)[i] = 0;
  }
  for(i = 0; i < features; i++){
    (m->fvNoise)[i] = (int*)malloc(sizeof(int)*featureValues);
    if(!(m->fvNoise)[i]) printMessages(0, NULL);
    for(j = 0; j < featureValues; j++) (m->fvNoise)[i][j] = 0;
    (m->fnoise)[i] = 0;
  }
  return m;
}

/* Copy values fro one model structure to another */
void copyModel(model *m, model *mc){
  int i, j, k;
  mc->arch = m->arch;
  mc->features = m->features;
  mc->featureValues = m->featureValues;
  mc->n = m->n;
  mc->lambda = m->lambda;
  mc->alpha = m->alpha;
  for(i = 0; i < mc->arch; i++){
    (mc->posCount)[i] = (m->posCount)[i];
    for(j = 0; j < mc->features; j++){
      (mc->pos)[i][j] = (m->pos)[i][j];
      for(k = 0; k < mc->featureValues; k++) (mc->count)[i][j][k] = (m->count)[i][j][k];
    }
    (mc->t)[i] = (m->t)[i];
  }
  for(i = 0; i < mc->features; i++){
    for(j = 0; j < mc->featureValues; j++) (mc->fvNoise)[i][j] = (m->fvNoise)[i][j];
    (mc->fnoise)[i] = (m->fnoise)[i];
  }
}

/* Deallocate space allotted to a model structure */
void freeModel(model *m){
  int i, j;
  for(i = 0; i < m->arch; i++){
    for(j = 0; j < m->features; j++){
      free((m->count)[i][j]); (m->count)[i][j] = NULL;
    }
    free((m->count)[i]); (m->count)[i] = NULL;
    free((m->pos)[i]); (m->pos)[i] = NULL;
  }
  for(i = 0; i < m->features; i++){
    free((m->fvNoise)[i]); (m->fvNoise)[i] = NULL;
  }
  free(m->posCount); m->posCount = NULL;
  free(m->pos); m->pos = NULL;
  free(m->count); m->count = NULL;
  free(m->t); m->t = NULL;
  free(m->fvNoise); m->fvNoise = NULL;
  free(m->fnoise); m->fnoise = NULL;
  m->posCount = NULL;
  free(m); m = NULL;
}


/* Read Fasta file and save data in a dataSet structure */
dataSet* getData(char *s, char *outFile){
  int features, n, tmpfeatures;
  FILE *dataFile, *fo;
  char c;
  dataSet *ds;
  int flag, ot, i, j, tu;
  dataFile = fopen(s, "r");
  if(!dataFile) printMessages(4, s);
  n = 0;
  tmpfeatures = 0;
  features = 0;
  flag = 0;
  /* Check if it is a valid Fasta file */
  while((c = fgetc(dataFile)) != EOF)
    if(c == '>' && flag == 0){ 
      n++;
      flag = 1;
    }
    else if(c == '\n' && flag == 1)
      flag = 0;
  rewind(dataFile);
  if(n == 0) printMessages(1, s);
  ot = 0;
  c = fgetc(dataFile);
  if(c != '>') printMessages(1, s);
  rewind(dataFile);
  tu = -1;
  flag = 0;
  while((c = fgetc(dataFile)) != EOF){
    if(c == '>' && flag == 0){
      if(ot == 0 && tmpfeatures > 0){
	ot = 1;
	features = tmpfeatures;
      }
      else if(ot == 1) 
	if(features != tmpfeatures) printMessages(2, s);
      flag = 1;
      tmpfeatures = 0;
    }
    else if(c == '\n' && flag == 1)
      flag = 0;
    else if(c == '\n' && flag == 0)
      continue;
    else if(((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) && flag == 0){
      switch(c){
      case 'A': tmpfeatures++; break;
      case 'a': tmpfeatures++; break;
      case 'C': tmpfeatures++; break;
      case 'c': tmpfeatures++; break;
      case 'G': tmpfeatures++; break;
      case 'g': tmpfeatures++; break;
      case 'T': tmpfeatures++; break;
      case 't': tmpfeatures++; break;
      default: printf("ERROR: Invalid symbol %c in file %s\n", c, s); exit(1);
      }
    }
    if((c == 'T' || c == 't') && tu == -1 && flag == 0) tu = 0;
    else if((c == 'U' || c == 'u') && tu == -1 && flag == 0) tu = 1;
    if((c == 'T' || c == 't') && tu == 1 && flag == 0) printMessages(3, s);
    else if((c == 'U' || c == 'u') && tu == 0 && flag == 0) printMessages(3, s);
  }

  if(ot == 1){
    if(features != tmpfeatures) printMessages(2, s);
  }
  else
    features = tmpfeatures;
  rewind(dataFile);
  ds = (dataSet*)malloc(sizeof(dataSet));
  if(!(ds)) printMessages(0, NULL);
  ds->data = (int**)malloc(sizeof(int*)*n);
  if(!(ds->data)) printMessages(0, NULL);
  for(i = 0; i < n; i++){
    (ds->data)[i] = (int*)malloc(sizeof(int)*features);
    if(!((ds->data)[i])) printMessages(0, NULL);
  }
  i = -1;
  flag = 0;
  j = 0;
  /* Save data in structure */
  while((c = fgetc(dataFile)) != EOF){
    if(c == '>' && flag == 0){
      i++;
      j = 0;
      flag = 1;
    }
    else if(c == '\n' && flag == 1)
      flag = 0;
    else if(c == '\n' && flag == 0)
      continue;
    else if(((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) && flag == 0){
      switch(c){
      case 'A': (ds->data)[i][j] = 0; j++; break;
      case 'a': (ds->data)[i][j] = 0; j++; break;
      case 'C': (ds->data)[i][j] = 1; j++; break;
      case 'c': (ds->data)[i][j] = 1; j++; break;
      case 'G': (ds->data)[i][j] = 2; j++; break;
      case 'g': (ds->data)[i][j] = 2; j++; break;
      case 'T': (ds->data)[i][j] = 3; j++; break;
      case 't': (ds->data)[i][j] = 3; j++; break;
      default: printf("ERROR: Invalid symbol %c in file %s\n", c, s); exit(1);
      }
    }
  }
  rewind(dataFile);
  fo = fopen(outFile, "w");
  flag = 0;
  ot = 0;
  while((c = fgetc(dataFile)) != EOF){
    if(c == '>' && flag == 0){
      if(ot != 0) fprintf(fo, "\n");
      else ot = 1;
      flag = 1;
    }
    else if(flag == 1 && c == '\n'){
      fprintf(fo, "\t");
      flag = 0;
    }
    else if(flag == 1){
      if(!(c <= 31 || c >=127))
	fprintf(fo, "%c", c);
    }
    else if(flag == 0 && (c == '\n' || c == ' ' || c == '\t')) continue;
    else{
      if(!(c <= 31 || c >=127))
      fprintf(fo, "%c", c);
    }
  }
  fclose(dataFile);
  fclose(fo);
  ds->n = n;
  ds->features = features;
  ds->featureValues = 4;
  ds->tu = tu;
  return ds;
}

/* Deallocate space allotted to a dataSet structure */
void freeData(dataSet *ds){
  int i;
  for(i = 0; i < ds->n; i++){
    free((ds->data)[i]); (ds->data)[i] = NULL;
  }
  free(ds->data); ds->data = NULL;
  free(ds); ds = NULL;
}

/* Deallocate space allotted to a trainOut structure */
void freeTo(trainOut *to){
  freeModel(to->m); to->m = NULL;
  free(to->lp); to->lp = NULL;
  free(to);
}

/* Deallocate space allotted to a trainOut structure except for the model structure inside */
void freeToXModel(trainOut *to){
  free(to->lp); to->lp = NULL;
  to->m = NULL;
  free(to);
}

/* Randomly sample labels for model with given number of architectures */
int* sampleLabels(int n, int arch, unsigned int *seed){
  int i, *lp;
  lp = (int*)malloc(sizeof(int)*n);
  if(!lp) printMessages(0, NULL);
  for(i = 0; i < n; i++) lp[i] = (rand_r(seed)%arch)+1;
  return lp;
}

/* Count the number of neucleotides and the number of sequences per architecture of a model */
void getCount(model *m, int *lp, dataSet *ds){
  int i, j;
  for(i = 0; i < m->n; i++){
    for(j = 0; j < m->features; j++) (m->count)[lp[i]-1][j][(ds->data)[i][j]]++;
    (m->t)[lp[i]-1]++;
  }
}

/* Randomly decide importance of positions in each architecture of given model */
void nRandom(int *pos, int posCount, int size, unsigned int *seed){
  int i, n, num, *numArr;
  numArr = (int*)malloc(sizeof(int)*size);
  if(!numArr) printMessages(0, NULL);
  n = size;
  for(i = 0; i < size; i++)
    numArr[i] = i;
  for(i = 0; i < posCount; i++){
    num = rand_r(seed)%n;
    pos[numArr[num]] = 1;
    numArr[num] = numArr[n-1];
    n--;
  }
  free(numArr);
}

/* Randomly decide importance of positions in each architecture of given model */
void nRandom1(int *pos, int posCount, int size, unsigned int *seed){
  int i, n, num, *numArr;
  numArr = (int*)malloc(sizeof(int)*size);
  if(!numArr) printMessages(0, NULL);
  n = size;
  for(i = 0; i < size; i++)
    numArr[i] = i;
  for(i = 0; i < posCount; i++){
    num = rand_r(seed)%n;
    pos[i] = numArr[num];
    numArr[num] = numArr[n-1];
    n--;
  }
  free(numArr);
}

void quickSort(int *x, int first, int last){
  int pivot, i, j, temp;
  if(first < last){
    pivot = first;
    i = first;
    j = last;
    while(i < j){
      while(x[i] <= x[pivot] && i < last) i++;
      while(x[j] > x[pivot]) j--;
      if(i < j){
	temp = x[i];
	x[i] = x[j];
	x[j] = temp;
      }
    }
    temp = x[pivot];
    x[pivot] = x[j];
    x[j] = temp;
    quickSort(x, first, j-1);
    quickSort(x, j+1, last);
  }
}

/* Initialize important positions per architecture of the model */
void initPositions(model *m, unsigned int *seed){
  int i;
  for(i = 0; i < m->arch; i++){
    (m->posCount)[i] = (rand_r(seed)%m->features) + 1;
    nRandom((m->pos)[i], (m->posCount)[i], m->features, seed);
  }
}

/* Compute number of features and the neucleotides in them that are considered noise in the architectures */
void getNoise(model *m){
  int i, j, k;
  for(i = 0; i < m->arch; i++)
    for(j = 0; j < m->features; j++)
      for(k = 0; k < m->featureValues; k++){
	(m->fvNoise)[j][k] = (m->fvNoise)[j][k] + (!(m->pos)[i][j])*(m->count)[i][j][k];
	(m->fnoise)[j] = (m->fnoise)[j] + (!(m->pos)[i][j])*(m->count)[i][j][k];
      }

}

/* First occurance of maximum value in array */
int maxPos(double *p, int n){
  int i, index;
  double max;
  index = 0;
  max = p[0];
  for(i = 1; i < n; i++)
    if(p[i] > max){
      max = p[i];
      index = i;
    }

  return index;
}

/* Sample values from weighted array */
int sample(double *p, int n, double r){
  int low, high, mid, i;
  double max, sum;
  low = 0;
  high = n;
  mid = n/2;
  max = p[0];
  for(i = 1; i < n; i++)
    max = p[i] > max ? p[i] : max;
  for(i = 0; i < n; i++)
    p[i] = pow(exp(1), p[i] - max);
  sum = 0;
  for(i = 0; i < n; i++)
    sum = sum + p[i];
  for(i = 0; i < n; i++)
    p[i] = p[i]/sum;
  for(i = 1; i < n; i++)
    p[i] = p[i] + p[i-1];
  while(low <= high){
    mid = (low + high)/2;
    if(r == p[mid]){
      if(mid >= n) return n-1;
      return mid;
    }
    else if(r > p[mid])
      low = mid + 1;
    else if(r < p[mid])
      high = mid - 1;
  }
  if(low >= n) return n-1;
  return low;
}

/* Save cross validation results */
void saveCVResults(model *m, char *s, double likelihood, int n, double l){
  int i;
  FILE *fp = fopen(s, "w");
  if(!(fp)){
    printf("ERROR: Could not open file %s. Exiting\n", s);
    exit(1);
  }
  fprintf(fp, "Length of sequences: %d\n", m->features);
  fprintf(fp, "Number of architectures: %d\n", m->arch);
  fprintf(fp, "Number of important positions per architecure: [%d", (m->posCount)[0]);
  for(i = 1; i < m->arch; i++)
    fprintf(fp, ", %d", (m->posCount)[i]);
  fprintf(fp, "]\n");
  fprintf(fp, "Number of sequences in trained set: %d\n", m->n);
  fprintf(fp, "Likelihood of trained set: %lf\n", likelihood);
  fprintf(fp, "Number of sequences in test set: %d\n", n);
  fprintf(fp, "Cross validation likelihood: %lf\n", l);
  fclose(fp);
}

/* Copy values from one array to another */
void copyLabels(int *lp1, int *lp2, int n){
  int i;
  for(i = 0; i < n; i++)
    lp2[i] = lp1[i];
}
