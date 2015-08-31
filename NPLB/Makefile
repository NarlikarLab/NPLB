
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

CC = gcc -g -O3
CFLAGS1 = -Wall -Wextra -ansi -pedantic-errors
CFLAGS = -fPIC -Wall
LDFLAGS = -shared
LM = -lm
RM = rm -f
TARGET_LIB = libdatacal.so
QSRCS = messages.c modelops.c traindata.c evaluate.c
VAR = -L$(CURDIR) -Wl,-rpath=$(CURDIR)
OBJS = $(QSRCS:.c=.o)

.PHONY: all

all: ${TARGET_LIB}
	cat cst1 > cstructures.py; echo "libctest = cdll.LoadLibrary(\"$(CURDIR)/libdatacal.so\")" >> cstructures.py; cat cst2 >> cstructures.py
	echo "python $(CURDIR)/pyt.py 1 "\"\$$@\"" $(CURDIR)/" > promoterLearn; chmod +x promoterLearn
	echo "python $(CURDIR)/pyt.py 0 "\"\$$@\"" $(CURDIR)/" > promoterClassify; chmod +x promoterClassify
	chmod +x chartPie.pg; chmod +x makeImage.pg; chmod +x makePlot.pg; chmod +x plotBox.pg

$(TARGET_LIB): $(OBJS) 
	   $(CC) ${LDFLAGS} -o $@ $^ -lm

.PHONY: clean
clean:
	-${RM} ${TARGET_LIB} ${OBJS} $(QSRCS:.c=.d) $(QFILE) *~ *.pyc cstructures.py promoterLearn promoterClassify