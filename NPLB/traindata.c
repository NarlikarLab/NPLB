
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
#include<math.h>
#include<string.h>
#include "messages.h"
#include "modelops.h"
#include "traindata.h"

/* Calls function to train model */
trainOut* callTrainData(dataSet *ds, int times, unsigned int seed, int arch, float alpha, double lambda, char *filename){
  trainDataStruct *tds;
  trainOut *to;
  tds = (trainDataStruct*)malloc(sizeof(trainDataStruct));
  if(!tds) printMessages(0, NULL);
  to = (trainOut*)malloc(sizeof(trainOut));
  if(!to) printMessages(0, NULL);

  tds->seed = seed;
  tds->filename = filename;
  tds->ds = ds;
  tds->times = times;

  tds->lp = sampleLabels((tds->ds)->n, arch, &(seed));



  tds->m = createModel(arch, tds->ds, tds->lp, alpha, &(tds->seed), lambda);
  trainData((void*)tds);
  to->m = tds->m;
  to->lp = tds->lp;
  (to->m)->tu = (tds->ds)->tu;
  to->likelihood = tds->likelihood;
  tds->ds = NULL;
  tds->m = NULL;
  tds->lp = NULL;
  free(tds);
  return to;
}

/* Assigns labels to data based on given model */
trainOut* callLearnData(dataSet *ds, model *m){
  int *lp;
  int i, j, k, pos;
  double *p;
  trainOut *to;
  p = (double*)malloc(sizeof(double)*m->arch);
  if(!p) printMessages(0, NULL);
  lp = (int*)malloc(sizeof(int)*ds->n);
  if(!lp) printMessages(0, NULL);
  for(i = 0; i < ds->n; i++){
    for(j = 0; j < m->arch; j++){
      p[j] = 0;
      pos = 0;
      /* Computes score */
      for(k = 0; k < m->features; k++){
	if((m->pos)[j][k] == 1){
	  p[j] = p[j] + log((m->count)[j][k][(ds->data)[i][k]] + m->alpha) - log((m->t)[j] + m->featureValues*m->alpha);
	  pos++;
	}
	else
	  p[j] = p[j] + log((m->fvNoise)[k][(ds->data)[i][k]] + m->alpha) - log((m->fnoise)[k] + m->featureValues*m->alpha);
      }
    }
    /* Assigns labels with maximum score */
    lp[i] = maxPos(p, m->arch) + 1;
  }
  free(p);
  /* Update model structure values according to the labels  */
  m->n = ds->n;
  for(i = 0; i < m->arch; i++){
    for(j = 0; j < m->features; j++)
      for(k = 0; k < m->featureValues; k++)
	(m->count)[i][j][k] = 0;
    (m->t)[i] = 0;
  }
  for(i = 0; i < m->features; i++){
    for(j = 0; j < m->featureValues; j++)
      (m->fvNoise)[i][j] = 0;
    (m->fnoise)[i] = 0;
  }
  getCount(m, lp, ds);
  getNoise(m);
  to = (trainOut*)malloc(sizeof(trainOut));
  if(!to) printMessages(0, NULL);
  to->lp = lp;
  to->likelihood = calculateLikelihood(m);
  to->m = m;
  return to;
}

/* Train data */
void *trainData(void *vtds){
  double maxLikelihood, tmpLikelihood, delta, checkPoint;
  int i, j, newLabel, flag, oldLabel;
  model *m;
  int *lp;
  unsigned int seed;
  trainDataStruct *tds;
  long long times, ct;
  FILE *fp;
  delta = log(pow(10, -6));
  tds = (trainDataStruct*)vtds;
  if((tds->filename)[0] == '0') fp = NULL;
  else{
    fp = fopen(tds->filename, "w");
    if(!fp) printMessages(4, tds->filename);
  }
  /* Importance of positions is irrelevant for model with one architecture */
  if((tds->m)->arch == 1){
    for(i = 0; i < (tds->m)->features; i++) ((tds->m)->pos)[0][i] = 1;
    ((tds->m)->posCount)[0] = (tds->m)->features;
    tds->likelihood = calculateLikelihood(tds->m);
    return 0;
  }

  /* Importance of positions is irrelevant when lambda is 0. Note it has been tested on data where all positions subsequently become important when lambda is set to 0.  */

  if((tds->m)->lambda == 0){
    for(i = 0; i < (tds->m)->arch; i++){
      for(j = 0; j < (tds->m)->features; j++) ((tds->m)->pos)[i][j] = 1;
     ((tds->m)->posCount)[i] = (tds->m)->features;
    }
  }
  m = initializeModel((tds->m)->arch, (tds->ds)->features, (tds->ds)->featureValues, (tds->ds)->n);
  copyModel(tds->m, m);
  lp = (int*)malloc(sizeof(int)*m->n);
  if(!lp) printMessages(0, NULL);
  copyLabels(tds->lp, lp, m->n);
  seed = tds->seed;
  maxLikelihood = calculateLikelihood(m);
  if(fp != NULL) fprintf(fp, "%lf\n", maxLikelihood);
  tmpLikelihood = maxLikelihood;
  j = 0;
  times = 0;
  ct = m->n + m->arch*m->features;
  /* Minimum number of iterations before exiting the loop */
  checkPoint = ((ct > 1000) ? 1000 : ct);
  ct = 10*ct;

  /* Loop till maximum likelihood stops increases for ct consecutive turns */

  while(times < ct){
    j++;
    if(j < checkPoint) times = 0;
    flag = 0;
    for(i = 0; i < m->n; i++){
      
      /* Apply Collapsed Gibbs Sampling to determine labels */

      addRemoveDataPoint(m, tds->ds, lp, i, -1);
      newLabel = sampleNewLabel(m, tds->ds, i, &(seed), 0, lp[i]-1);
      flag = newLabel != lp[i];
      oldLabel = lp[i];
      lp[i] = newLabel;
      addRemoveDataPoint(m, tds->ds, lp, i, 1);
      if(newLabel != oldLabel){

	/* Update likelihood on change of label */

	tmpLikelihood = calculateLikelihoodLabel(m, tmpLikelihood, oldLabel-1, newLabel-1, ((tds->ds)->data)[i]);
	if(fp != NULL) fprintf(fp, "%lf\n", tmpLikelihood);
      	if(tmpLikelihood > maxLikelihood){
      	  times = 0;
      	  maxLikelihood = tmpLikelihood;
      	  copyModel(m, tds->m);
      	  copyLabels(lp, tds->lp, m->n);
	  
	  
      	}
      	else times++;
      }
      else times++;
    }

    /* Sample the importance of features as 0 or 1 */

    if(m->lambda != 0)
      for(i = 0; i < m->arch; i++)
	sampleImpFeature(m, i, &(seed), 0, &tmpLikelihood, tds, &maxLikelihood, lp, &times, fp);
    if(j == tds->times) break;
  }
  copyLabels(tds->lp, lp, m->n);
  copyModel(tds->m, m);
  tmpLikelihood = maxLikelihood;
  flag = 0;

  /* Apply EM like method for convergence */

  while(flag < m->n){
    flag = 0;
    for(i = 0; i < m->n; i++){
      newLabel = sampleNewLabel(m, tds->ds, i, NULL, 1, -1);
      if(newLabel != lp[i]){
  	addRemoveDataPoint(m, tds->ds, lp, i, -1);
  	oldLabel = lp[i];
  	lp[i] = newLabel;
  	addRemoveDataPoint(m, tds->ds, lp, i, 1);
	tmpLikelihood = calculateLikelihoodLabel(m, tmpLikelihood, oldLabel-1, newLabel-1, ((tds->ds)->data)[i]);
	if(fp != NULL) fprintf(fp, "%lf\n", tmpLikelihood);
  	if(tmpLikelihood > maxLikelihood){
  	  maxLikelihood = tmpLikelihood;
  	  copyModel(m, tds->m);
  	  copyLabels(lp, tds->lp, m->n);
  	}
      }
      else flag++;
    }
    if(m->lambda != 0)
      for(i = 0; i < m->arch; i++){
	flag -= sampleImpFeature(m, i, NULL, 1, &tmpLikelihood, tds, &maxLikelihood, lp, NULL, fp);
      }
  }
  tds->likelihood = maxLikelihood;
  free(lp);
  freeModel(m);
  if(fp != NULL) fclose(fp);
  m = NULL;
  return 0;
}

/* Likelihood computation for the model */

double calculateLikelihood(model *m){
  double s;
  double tmp1;
  int i, j, k;
  s = 0;
  tmp1 = 0;
  for(i = 0; i < m->arch; i++)
    s += (m->posCount)[i];
  s = s * (-(m->lambda));
  for(i = 0; i < m->arch; i++){
    for(j = 0; j < m->features; j++)
      for(k = 0; k < m->featureValues; k++)
  	if((m->pos)[i][j] == 1)
  	  s = s + ((m->count)[i][j][k]*(log((m->count)[i][j][k] + m->alpha) - log((m->t)[i] + m->featureValues*m->alpha)));
    s = s + log((m->t)[i] + m->alpha) - log(m->n + m->arch*m->alpha);
  }
  for(i = 0; i < m->features; i++)
    for(j = 0; j < m->featureValues; j++){
      s = s + ((m->fvNoise)[i][j])*(log((m->fvNoise)[i][j] + m->alpha) - log((m->fnoise)[i] + m->featureValues*m->alpha));
      tmp1 = tmp1 + ((m->fvNoise)[i][j])*(log((m->fvNoise)[i][j] + m->alpha) - log((m->fnoise)[i] + m->featureValues*m->alpha));
    }

  return s;
}

/* Update likelihood based on change on importance of positions */

double calculateLikelihoodLabel(model *m, double likelihood, int oldArch, int newArch, int *data){
  int i, j;
  double tmp1, tmp2;
  tmp1 = 0;
  tmp2 = 0;
  for(i = 0; i < m->features; i++){
    if((m->pos)[oldArch][i] == 1){
      for(j = 0; j < m->featureValues; j++){
	likelihood = likelihood - ((m->count)[oldArch][i][j] + (data[i] == j))*(log((m->count)[oldArch][i][j] + (data[i] == j) + m->alpha) - log((m->t)[oldArch] + 1 + m->featureValues*m->alpha));
	likelihood = likelihood + ((m->count)[oldArch][i][j])*(log((m->count)[oldArch][i][j] + m->alpha) - log((m->t)[oldArch] + m->featureValues*m->alpha));
      }
      if((m->pos)[newArch][i] == 0)
	for(j = 0; j < m->featureValues; j++){
	  likelihood = likelihood - ((m->fvNoise)[i][j] - (data[i] == j))*(log((m->fvNoise)[i][j] - (data[i] == j) + m->alpha) - log((m->fnoise)[i] - 1 + m->featureValues*m->alpha));
	  likelihood = likelihood + ((m->fvNoise)[i][j])*(log((m->fvNoise)[i][j] + m->alpha) - log((m->fnoise)[i] + m->featureValues*m->alpha));
	}
    }
    if((m->pos)[newArch][i] == 1){
      for(j = 0; j < m->featureValues; j++){
    	likelihood = likelihood - ((m->count)[newArch][i][j] - (data[i] == j))*(log((m->count)[newArch][i][j] - (data[i] == j) + m->alpha) - log((m->t)[newArch] - 1 + m->featureValues*m->alpha));
    	likelihood = likelihood + ((m->count)[newArch][i][j])*(log((m->count)[newArch][i][j] + m->alpha) - log((m->t)[newArch] + m->featureValues*m->alpha));
      }
      if((m->pos)[oldArch][i] == 0){
	for(j = 0; j < m->featureValues; j++){
	  likelihood = likelihood - ((m->fvNoise)[i][j] + (data[i] == j))*(log((m->fvNoise)[i][j] + (data[i] == j) + m->alpha) - log((m->fnoise)[i] + 1 + m->featureValues*m->alpha));
	  likelihood = likelihood + ((m->fvNoise)[i][j])*(log((m->fvNoise)[i][j] + m->alpha) - log((m->fnoise)[i] + m->featureValues*m->alpha));
	}
      }
    }
  }
  likelihood = likelihood - (log((m->t)[oldArch] + 1 + m->alpha) - log(m->n + m->arch*m->alpha));
  likelihood = likelihood + (log((m->t)[oldArch] + m->alpha) - log(m->n + m->arch*m->alpha));
  likelihood = likelihood - (log((m->t)[newArch] - 1 + m->alpha) - log(m->n + m->arch*m->alpha));
  likelihood = likelihood + (log((m->t)[newArch] + m->alpha) - log(m->n + m->arch*m->alpha));
  return likelihood;
}


/* Update likelihood based on neucleotide counts in important and unimportant positions of the given architecture */

void calculateLikelihoodPos(model *m, double *likelihood, int zi, int pos){
  int i;
  double tmp1, tmp2;
  tmp1 = tmp2 = 0;
  if((m->pos)[zi][pos] == 1){
    *likelihood = *likelihood - (m->lambda);
    for(i = 0; i < m->featureValues; i++){
      *likelihood = *likelihood - ((m->fvNoise)[pos][i] + (m->count)[zi][pos][i])*(log((m->fvNoise)[pos][i] + (m->count)[zi][pos][i] + m->alpha) - log((m->fnoise)[pos] + (m->t)[zi] + m->featureValues*m->alpha));
      tmp1 += ((m->fvNoise)[pos][i] + (m->count)[zi][pos][i])*(log((m->fvNoise)[pos][i] + (m->count)[zi][pos][i] + m->alpha) - log((m->fnoise)[pos] + (m->t)[zi] + m->featureValues*m->alpha));
      *likelihood = *likelihood + ((m->fvNoise)[pos][i])*(log((m->fvNoise)[pos][i] + m->alpha) - log((m->fnoise)[pos] + m->featureValues*m->alpha));
      *likelihood = *likelihood + (m->count)[zi][pos][i]*(log((m->count)[zi][pos][i] + m->alpha) - log((m->t)[zi] + m->featureValues*m->alpha));
    }
  }
  else{
    *likelihood = *likelihood + (m->lambda);
    for(i = 0; i < m->featureValues; i++){
      *likelihood = *likelihood - ((m->fvNoise)[pos][i] - (m->count)[zi][pos][i])*(log((m->fvNoise)[pos][i] - (m->count)[zi][pos][i] + m->alpha) - log((m->fnoise)[pos] - (m->t)[zi] + m->featureValues*m->alpha));
      *likelihood = *likelihood - (m->count)[zi][pos][i]*(log((m->count)[zi][pos][i] + m->alpha) - log((m->t)[zi] + m->featureValues*m->alpha));
      *likelihood = *likelihood + ((m->fvNoise)[pos][i])*(log((m->fvNoise)[pos][i] + m->alpha) - log((m->fnoise)[pos] + m->featureValues*m->alpha));
    }
  }
}

/* Add/Remove data point and update neucleotide counts accordingly */

void addRemoveDataPoint(model *m, dataSet *ds, int *lp, int index, int ar){
  int i, label;
  label = lp[index]-1;
  for(i = 0; i < m->features; i++){
    (m->count)[label][i][(ds->data)[index][i]] = (m->count)[label][i][(ds->data)[index][i]] + ar;
  }
  (m->t)[label] = (m->t)[label] + ar;
  m->n = m->n + ar;
  for(i = 0; i < m->features; i++)
    if((m->pos)[label][i] == 1) continue;
    else{
      (m->fvNoise)[i][(ds->data)[index][i]] = (m->fvNoise)[i][(ds->data)[index][i]] + ar;
      (m->fnoise)[i] = (m->fnoise)[i] + ar;
    }
}

/* Sample label based on the computed weights */

int sampleNewLabel(model *m, dataSet *ds, int index, unsigned int *seed, int flag, int label){
  int i, j;
  double *values;
  values = (double*)malloc(sizeof(double)*m->arch);
  if(!values) printMessages(0, NULL);
  for(i = 0; i < m->arch; i++){
    values[i] = log((m->t)[i] + m->alpha) - log(m->n + m->arch*m->alpha);
    for(j = 0; j < m->features; j++)
      if((m->pos)[i][j] == 1) values[i] = values[i] + log((m->count)[i][j][(ds->data)[index][j]] + m->alpha) - log((m->t)[i] + m->featureValues*m->alpha); 
      else values[i] = values[i] + log((m->fvNoise)[j][(ds->data)[index][j]] + m->alpha) - log((m->fnoise)[j] + m->featureValues*m->alpha);
  }
  if(flag == 0) j = sample(values, m->arch, ((double)rand_r(seed))/(RAND_MAX))+1;
  else j = maxPos(values, m->arch) + 1;
  free(values);
  return j;
}


/* Sample importance of features */

int sampleImpFeature(model *m, int zi, unsigned int *seed, int flag, double *tmpLikelihood, trainDataStruct *tds, double *maxLikelihood, int *lp, long long *times, FILE *fp){
  int i, j;
  int newLabel;
  int changed;
  changed = 0;
  for(i = 0; i < m->features; i++){
    newLabel = solvePosExp(m, zi, i, seed, flag);
    if(newLabel != (m->pos)[zi][i]){
      changed = 1;

      /* Update neucleotide counts when change in importance */

      for(j = 0; j < m->featureValues; j++){
	if(newLabel == 1) (m->fvNoise)[i][j] = (m->fvNoise)[i][j] - (m->count)[zi][i][j];
	else (m->fvNoise)[i][j] = (m->fvNoise)[i][j] + (m->count)[zi][i][j];
      }
      if(newLabel == 1) (m->fnoise)[i] = (m->fnoise)[i] - (m->t)[zi];
      else (m->fnoise)[i] = (m->fnoise)[i] + (m->t)[zi];
      if(newLabel == 1) (m->posCount)[zi]++;
      else (m->posCount)[zi]--;
      (m->pos)[zi][i] = newLabel;
      calculateLikelihoodPos(m, tmpLikelihood, zi, i);
      if(fp != NULL) fprintf(fp, "%lf\n", *tmpLikelihood);
      if(*tmpLikelihood > *maxLikelihood){
      	if(flag == 0) *times = 0;
      	*maxLikelihood = *tmpLikelihood;
      	copyModel(m, tds->m);
      	copyLabels(lp, tds->lp, m->n);
      }
      else if(flag == 0) *times = *times + 1;
    }
    else if(flag == 0) *times = *times + 1;
  }
  return changed;
}

/* Compute score for imporance of position in a given architecture */

int solvePosExp(model *m, int zi, int pos, unsigned int *seed, int flag){
  int i, j, k;
  double *s, ts;
  s = (double*)malloc(sizeof(double)*2);
  if(!s) printMessages(0, NULL);
  s[0] = -(m->lambda);
  s[1] = 0;
  for(i = 0; i < m->featureValues; i++){
    s[0] = s[0] + (m->count)[zi][pos][i]*(log((m->count)[zi][pos][i] + m->alpha) - log((m->t)[zi] + m->featureValues*m->alpha));
    s[0] = s[0] + ((m->fvNoise)[pos][i] - (!(m->pos)[zi][pos])*(m->count)[zi][pos][i])*(log((m->fvNoise)[pos][i] - (!(m->pos)[zi][pos])*(m->count)[zi][pos][i] + m->alpha) - log((m->fnoise)[pos] - (!(m->pos)[zi][pos])*(m->t)[zi] + m->featureValues*m->alpha));
    s[0] = s[0] - ((m->fvNoise)[pos][i] + (m->pos)[zi][pos]*(m->count)[zi][pos][i])*(log((m->fvNoise)[pos][i] + (m->pos)[zi][pos]*(m->count)[zi][pos][i] + m->alpha) - log((m->fnoise)[pos] + (m->pos)[zi][pos]*(m->t)[zi] + m->featureValues*m->alpha));
  }
  ts = s[0];
  if(ts == 0 && (m->t)[zi] > 0){
    j = 0;
    k = 0;
    for(i = 0; i < m->arch; i++)
      if((m->pos)[i][pos] == 1 && i != zi && (m->t)[i] > 0) j++;
    for(i = 0; i < m->arch; i++)
      if(i != zi && (m->t)[i] > 0) k++;
    if(j != k){
      fprintf(stderr, "ERROR: Important position score becoming zero!\n");
    }
  }
  if(flag == 0) i = sample(s, 2, ((float)rand_r(seed))/(RAND_MAX));
  else i = maxPos(s, 2);
  free(s);
  if(i == 0) return 1;
  else return 0;
  return i;
}


/* Compute cross validation likelihood */

double cvLikelihood(dataSet *testSet, model *m){
  int i, j, k, pos;
  double s, max, sum, *p;
  p = (double*)malloc(sizeof(double)*m->arch);
  if(!p) printMessages(0, NULL);
  sum = 0;
  for(i = 0; i < testSet->n; i++){
    for(j = 0; j < m->arch; j++){
      p[j] = 0;
      pos = 0;
      for(k = 0; k < testSet->features; k++){
	p[j] += (m->pos)[j][k]*(log((m->count)[j][k][(testSet->data)[i][k]] + m->alpha) - log((m->t)[j] + m->featureValues*m->alpha)) + !(m->pos)[j][k]*(log((m->fvNoise)[k][(testSet->data)[i][k]] + m->alpha) - log((m->fnoise)[k] + m->featureValues*m->alpha));
      }
      p[j] += log((m->t)[j] + m->featureValues*m->alpha) - log(m->n + m->arch*m->featureValues*m->alpha);
    }
    max = p[0];
    for(j = 1; j < m->arch; j++)
      max = p[j] > max ? p[j] : max;
    s = 0;
    for(j = 0; j < m->arch; j++)
      s = s + pow(exp(1), p[j] - max);
    s = log(s) + max;
    sum += s;
  }
  free(p);
  return sum;
}
