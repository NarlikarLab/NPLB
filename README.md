# NPLB
No Promoter Left Behind - learn de novo promoter architectures from genome-wide TSSs

## Prerequisites
The following packages need to be installed in order to run NPLB:

- Python 2.6+ (Not compatible with Python 3.x)
- python-numpy
- python-ctypes
- python-multiprocessing
- python-re
- python-pickle
- gnuplot 4.6+
- Truetype fonts

## Installation
NPLB is freely available at https://github.com/computationalBiology/NPLB/. Execute the following commands to download and install NPLB:
```sh
wget "https://github.com/computationalBiology/NPLB/archive/v1.0.0.tar.gz"
tar -xvf v1.0.0.tar.gz 
cd NPLB-1.0.0/NPLB/
make
```

To execute NPLB from anywhere export the path to NPLB to the PATH variable:
```sh
export PATH=$PATH:/path/to/NPLB/
```

**Note:** gnuplot font path should be set the directory containing truetype fonts.
It can be set as:
```sh
GDFONTPATH=$GDFONTPATH/path/to/truetype/fonts/
```
## Usage of promoterLearn
To learn new models using promoterLearn:
```sh
/path/to/NPLB/promoterLearn [options]
```
### Options
The various options for running *promoterLearn* are as follows:
- ````-f filename````
Compulsory. Data file for which hidded promoter architectures are to be identified. File must be in fasta format and must consist of sequences of equal length.
- ````-o directory````
Valid directory name. If it exists then a new one is created with given name along with an extension number. Default directory: NPLBoutput_extension in the current working directory.
- ````-a positive_real_number````
Pseudocount. Default value: 1.
- ````-t times````
Maximum number of iterations over the entire dataset. Default: It iterates at most 1000 times and then loops until the likelihood stops increasing for 10*n consecutive iterations where n is the sum of number of sequences and number of architectures times length of sequences.
- ````-i 0 or 1````
Flag whether to save image matrix in the given directory. Default value: 1(save)
- ````-v 0 or 1````
Flag whether to save likelihood plots of the models learned in each fold. Default value: 0(don't save)
- ````-kfold k````
Value of k in K-fold cross validation. Default value: 5
- ````-lambda l````
Non negative real number. It is the regularization parameter used for determining the importance of features. Default value: Learns models by varying lambda value and chooses the best one.
- ````-lcount lc````
Number of models to be learnt while training. Only the best model is considered. Default value: 5
- ````-minarch mn````
Minimum number of architectures possible for the given dataset. Default value: 1
- ````-maxarch mx````
Maximum number of architectures possible for the given dataset. Default value: 20
- ````-proc pr````
Maximum number of processors to be used for computation. Default value is the number of processors the system has.
- ````-plotExtra filename````
Plots data from given tab separated file as pie charts or boxplots (depending on the type of data) for each architecture of the best model. **Note:** -pCol must also be set.
- ````-pCol pc````
Natural number. Plots the data in the given column of the given file in the form of pie charts or boxplots with respect to the architectures of the best model. **Note:** -plotExtra must specify a valid filename
- ````-sortBy sb````
Natural number. Sort architectures in increasing order of median values calculated from values in column -sortBy of file -plotExtra. **Note:** -plotExtra must specify a valid filename and -pCol must be set.

### Examples
The following examples illustrate the usage of the options.
- To run with all the default options on dataset example.fa:
```sh 
promoterLearn -f example.fa
```
- To run for architectures ranging from 2 to 5 with 100 iterations in directory Output:
```sh
promoterLearn -f example.fa -minarch 2 -maxarch 5 -t 100 -o Output
```
- To run it for architectures ranging from 15 to 25 and save likelihood plots:
```sh
promoterLearn -f example.fa -minarch 15 -maxarch 25 -v 1
```
- To plot pie charts in fourth column of file example.info.bed and rearrange architectures according to the sixth column:
```sh
promoterLearn -f example.fa -plotExtra example.info.bed -pCol 4 -sortBy 6
```
- To find the list of options:
```sh
promoterLearn
```
## Usage of promoterClassify
To use an existing model to learn new PAs using promoterClassify:
```sh
/path/tp/NPLB/promoterClassify [options]
```
### Options
The various options for running *promoterClassify* are as follows:
- ````-f filename````
Compulsory. Data file for which hidded promoter architectures are to be identified. File must be in fasta format and must consist of sequences of equal length.
- ````-m model````
Compulsory. Model learned by executing promoterLearn on a given dataset. Model must have sequences of the same length as the input fasta file.
- ````-o directory````
Valid directory name. If it exists then a new one is created with given name along with an extension number. Default directory: NPLBoutput_extension in the current working directory.
- ````-i 0 or 1````
Flag whether to save image matrix in the given directory. Default value: 1(save)
- ````-plotExtra filename````
Plots data from given tab separated file as pie charts or boxplots (depending on the type of data) for each architecture of the best model. **Note:** -pCol must also be set.
- ````-pCol pc````
Natural number. Plots the data in the given column of the given file in the form of pie charts or boxplots with respect to the architectures of the best model. **Note:** -plotExtra must specify a valid filename
- ````-sortBy sb````
Natural number. Sort architectures in increasing order of median values calculated from values in column -sortBy of file -plotExtra. **Note:** -plotExtra must specify a valid filename and -pCol must be set.

### Examples
The following examples illustrate the usage of the options.
- To find labels for dataset example2.fa from model present in directory Output as bestmodel.p:
```sh
promoterClassify -f example2.fa -m Output/bestmodel.p
```
- To find labels for dataset example2.fa from model bestmodel.p and plot pie chart or boxplot for column 4 of example3.info.bed:
```sh
promoterClassify -f example2.fa -m bestmodel.p -plotExtra example3.info.bed -pCol 4
```
- To find the list of options:
```sh
promoterClassify
```

License
----

GNU General Public License v3.0

