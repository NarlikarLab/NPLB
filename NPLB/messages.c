
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

/* Print error messages and exit program */
void printMessages(int i, char *s){
  switch(i){
  case 0: printf("ERROR: Invalid malloc operation\n"); exit(1); break;
  case 1: printf("ERROR: Invalid fasta file %s\n", s); exit(1); break;
  case 2: printf("ERROR: Unequal number of sequences in file %s\n", s); exit(1); break;
  case 3: printf("ERROR: Invalid fasta file %s. Cannot have both T's and U's\n", s); exit(1); break;
  case 4: printf("ERROR: Could not open file %s\n", s); exit(1); break;
  }
}
