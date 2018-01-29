#!/lustre/work/apps/Anaconda3/python3

#######################################################################
#Need to change the spaces to tabs

# overlap.py
#
# version 0.1
# 2017-09-06
#
#  This script reads in RepeatMasker output and then identifies
#   overlapping insertions.  This version of the script prints
#   out all non-overlapping annotations and when two anotations
#   overlap, the longest annotation is printed.
#
#  input: <repeat_masker_file>
#  output: <repeat_masker_file_no_overlaps>
#
#  input and output files are currently designated in the script
#
#
########################################################################

#  TO DO LIST
#----------------------
#
# [2] When comparing insertions, do we 
#       (A) merge the two annotations (if they are the same repeat)
#       (B) remove one of the annotations
#       (C) consider them as two repeats and assign overlaping area to o/o
#      Will need to consider things like percent identity, TE family, etc...
#
# [3] How are multiple overlapping repeats treated.
#      EX.
#               #----1----#
#                  #----------2---------#
#                       #-------3----------#
#
#     The current program compares insertions in a serial manner. So even
#       though #1 and #3 overlap, they are not directly compared to e/o
#       in the current iteration of the script (v0.1)
#
# [4] There are some issues with the way the last line of the file is being
#       handled.  It is being conisdered as a RM hit
#
# [5] We delete the first three lines of the RM input since it contains
#       header info, rather than assuming and deleting those lines we should
#       add a check to see if they are there or not.
#

import os
import re
import sys
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#from collections import defaultdict
#from Bio.Alphabet import generic_dna
import argparse
import random

#functions
def write_hit():
        "This function prints the hit"
        OUTPUT.write(
                str(SW_score)           + "\t" +        str(perc_div)           + "\t" +
                str(perc_del)           + "\t" +        str(perc_ins)           + "\t" +
                str(query_sequence)     + "\t" +        str(q_begin)            + "\t" +
                str(q_end)              + "\t" +        str(q_left)             + "\t" +
                str(orient)             + "\t" +        str(matching_repeat)    + "\t" +
                str(class_family)       + "\t" +        str(r_begin)            + "\t" +
                str(r_end)              + "\t" +        str(r_left)                       + "\t" + 
                str(ID)                 + "\n"
        )


def write_hit_i():
        "This function prints the hit_i"
        OUTPUT.write(
                str(SW_score_i)         + "\t" +        str(perc_div_i)         + "\t" +
                str(perc_del_i)         + "\t" +        str(perc_ins_i)         + "\t" +
                str(query_sequence_i)   + "\t" +        str(q_begin_i)          + "\t" +
                str(q_end_i)            + "\t" +        str(q_left_i)           + "\t" +
                str(orient_i)           + "\t" +        str(matching_repeat_i)  + "\t" +
                str(class_family_i)     + "\t" +        str(r_begin_i)          + "\t" +
                str(r_end_i)            + "\t" +        str(r_left_i)           + "\t" +
                str(ID_i)
        )



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 0 - Setting Arguments
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_args():
        parser=argparse.ArgumentParser()
        parser.add_argument('-i', '--repeatmasker', type=str, 
                                                help='Repeat Masker output file', required=True)
        parser.add_argument('-l', '--location', type=str, help='location to start', 
                                                required=True)
        parser.add_argument('-o', '--output', type=str, help='Output file name')
        args=parser.parse_args()
        RM=args.repeatmasker
        LOC=args.location
        OUT=args.output

        return RM, LOC, OUT

RM, LOC, OUT = get_args()

#Sanity Checks!!!!
print('The repeatMasker output file is', RM)
print('Directory to begin is', LOC)
print('The output file name is', OUT)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1 - initializing
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#set to the workingdirectory and open the output file
os.chdir(LOC)

#Generate Output File
OUTPUT=open(LOC + '/' + OUT,'w')

#Read in input file
rm_input = open(LOC + '/' + RM, 'r')

#Generate category specific output files:
CAT1OUT=open(LOC + '/' + 'cat1.out', 'w')
CAT2OUT=open(LOC + '/' + 'cat2.out', 'w')
CAT3OUT=open(LOC + '/' + 'cat3.out', 'w')
CAT4OUT=open(LOC + '/' + 'cat4.out', 'w')
CAT5OUT=open(LOC + '/' + 'cat5.out', 'w')

#Create a Lengths file that prints how many basepairs are in each overlap.
LENGTHS=open(LOC + '/' + 'lengths.csv', 'w')

#setting an iterative counter to identify the number of lines in the RM input
HIT_ENTRY = 1

#declaring empty list; all repeat makser hits will be stored in this list
HIT_ARRAY = []

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2. Storing RM data in memory
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# cycle through repeat masker input and load each record into a list 
#   (stored in mem)
for HIT in rm_input:

        #RM records start on line 3 after the header info
        if HIT_ENTRY > 3:

                #format record so that its tab delimited and no leading/trailing spaces
                HIT = re.sub(' +', '\t', HIT).rstrip().lstrip()

                #appending record to list of all other RM records
                HIT_ARRAY.append(HIT)

        HIT_ENTRY = HIT_ENTRY + 1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 3 Cycling through RM to identify overlaps
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#iterator
i = 0

#examine each "HIT" (RM insertion) in the "HIT_ARRAY" (repeat_masker file)
for HIT in HIT_ARRAY:

        #cleaning HIT of all white spaces and making tab delimited (to split on)
        HIT = HIT.lstrip()
        HIT=re.sub('\s+','\t',HIT)

        #splitting the current "HIT" on tabs to get the RM data
        (SW_score,              perc_div,               perc_del,
         perc_ins,              query_sequence,         q_begin,
         q_end,                 q_left,                 orient,
         matching_repeat,       class_family,           r_begin,
         r_end,                 r_left,                 ID)     =HIT.split("\t")
        #splitting the next record in the RM data (these will be compared to e/o)
        (SW_score_i,            perc_div_i,             perc_del_i,
         perc_ins_i,            query_sequence_i,       q_begin_i,
         q_end_i,               q_left_i,               orient_i,
         matching_repeat_i,     class_family_i,         r_begin_i,
        r_end_i,                r_left_i,               ID_i)   =HIT_ARRAY[i+1].split("\t")

        #forcing the data into proper types
        q_begin, q_end, q_begin_i, q_end_i=int(q_begin), int(q_end), int(q_begin_i), int(q_end_i)
        perc_div, perc_div_i=float(perc_div), float(perc_div_i)

        # q_begin, q_end, etc.. variables store the start and stop coordinates of the repeat
        # query_sequence is chromosome

        #calulate the lengths of each repeat (and other metrics) to determine which
        # insertion will be saved or how they will be treated
        length_q=q_end-q_begin
        length_q_i=q_end_i-q_begin_i

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Step 4 Identify overlaps to Category
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Overlapping repeat annotations can occur in 5 different ways
        #  defined as "categories" here

        CAT1=(q_begin < q_begin_i) and (q_end > q_begin_i) and (q_end < q_end_i) and (query_sequence == query_sequence_i)
        CAT2=(q_begin > q_begin_i) and (q_end < q_end_i) and ( query_sequence == query_sequence_i)
        CAT3=(q_begin < q_begin_i) and (q_end > q_end_i) and (query_sequence == query_sequence_i)
        CAT4=(q_begin > q_begin_i) and (q_end > q_end_i) and (q_begin < q_end_i) and (query_sequence == query_sequence_i)
        CAT5=(q_begin == q_begin_i) and (q_end == q_end_i) and (query_sequence == query_sequence_i)

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # CATEGORY 1 OVERLAP
        #
        # chrX   |q_begin----------- q -----------q_end|
        # chrX            |q_begin_i------------- qi --------q_end_i|

        if CAT1:

                #BASIC STATS
                #Print the number of base pairs in overlap to LENGTHS file
                CAT1_LENGTH=q_end-q_begin_i
                LENGTHS.write(str(CAT1_LENGTH) + "\n")

                #Print overlap to category 1 file
                CAT1OUT.write(HIT + "\n" + HIT_ARRAY[i+1] + "\n")

                if perc_div < perc_div_i:
                        q_begin_i=q_end+1

                elif perc_div_i < perc_div:
                        q_end=q_begin_i-1

                else:
                        if length_q > length_q_i:
                                q_begin_i=q_end+1

                        elif length_q_i > length_q:
                                q_end=q_begin_i-1

                        else:
                                k=random.randint(0,1)

                                if k == 0:
                                        q_begin_i=q_end+1
                                else:
                                        q_end=q_begin_i-1
                length_q=q_end-q_begin
                length_q_i=q_end_i-q_begin_i

                if length_q>1:
                        write_hit()
                if length_q_i>1:
                        write_hit_i()

                #delete q_i from future comparisons
                del HIT_ARRAY[i+1]

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # CATEGORY 2 OVERLAP
        #
        # chrX                 |q_begin------- q ------q_end|
        # chrX            |q_begin_i----------  qi --------q_end_i|


        elif CAT2:
#               print ("Comparing ",SW_score, "to ", SW_score_i, ". ","cat_2")

                #Print the number of base pairs in overlap to LENGTHS file
                CAT2_LENGTH=q_end-q_begin
                LENGTHS.write(str(CAT2_LENGTH) + "\n")

                #Print overlap to category 2 file
                CAT2OUT.write(HIT + "\n" + HIT_ARRAY[i+1] + "\n")

                #if cat2: print the longest RM hit
                if length_q > length_q_i:
                        OUTPUT.write(HIT)
                else:
                        OUTPUT.write(HIT_ARRAY[i+1])
                #remove overlap from additional comparisons
                del HIT_ARRAY[i+1]

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # CATEGORY 3 OVERLAP
        #
        # chrX   |q_begin----------- q -------------------q_end|
        # chrX            |q_begin_i----qi ------q_end_i|

        elif CAT3:

                #Print the number of base pairs in overlap to LENGTHS file
                CAT3_LENGTH=q_end_i-q_begin_i
                LENGTHS.write(str(CAT3_LENGTH) + "\n")

                #Print overlap to category 2 file
                CAT3OUT.write(HIT + "\n" + HIT_ARRAY[i+1] + "\n")

                #if cat3: print longest
                if length_q > length_q_i:
                        OUTPUT.write(HIT)
                else:
                        OUTPUT.write(HIT_ARRAY[i+1])
                #remove overlap from additional comparisons
                del HIT_ARRAY[i+1]

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # CATEGORY 4 OVERLAP
        #
        # chrX          |q_begin----------- q -----------q_end|
        # chrX   |q_begin_i----- qi -----q_end_i|

        elif CAT4:

                #We need to see if this overlap is significat
                # ex. 2bp vs 2,000 bp overlap
                #OVERLAP=q_end_i-q_begin
                CAT4_LENGTH=q_end_i-q_begin
                LENGTHS.write(str(CAT4_LENGTH) + "\n")

                CAT4OUT.write(HIT + "\n" + HIT_ARRAY[i+1] + "\n")

                #if the overlap is greater than 10 bp (it is significant)
                # and we need identify the longest one
                if perc_div < perc_div_i:
                        q_end_i=q_begin-1

                elif perc_div_i < perc_div:
                        q_begin=q_end_i+1

                else:
                        if length_q > length_q_i:
                                q_end_i=q_begin-1

                        elif length_q_i > length_q:
                                q_begin=q_end_i+1

                        else:
                                k=random.randint(0,1)

                                if k == 0:
                                        q_end_i=q_begin-1
                                else:
                                        q_begin=q_end_i+1


                length_q=q_end-q_begin
                length_q_i=q_end_i-q_begin_i

                if length_q>1:
                        write_hit()
                if length_q_i>1:
                        write_hit_i()

                #delete q_i from future comparisons
                del HIT_ARRAY[i+1]

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # CATEGORY 5 OVERLAP
        #
        # chrX   |q_begin---------- q -----------q_end|
        # chrX   |q_begin_i---------qi --------q_end_i|

        elif CAT5:
#               print ("Comparing", SW_score, "to", SW_score_i, ".", "cat_5")

                #Print the number of base pairs in overlap to LENGTHS file
                CAT5_LENGTH=q_end-q_begin
                LENGTHS.write(str(CAT5_LENGTH) + "\n")

                #Print overlap to category 2 file
                CAT5OUT.write(HIT + "\n" + HIT_ARRAY[i+1] + "\n")

                #if cat 5: choose the element most similar to consensus
                if perc_div > perc_div_i:
                        OUTPUT.write(HIT)
                else:
                        OUTPUT.write(HIT_ARRAY[i+1])
                #remove overlap from additional comparisons
                del HIT_ARRAY[i+1]

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # DO NOT OVERLAP
        #
        else:
        # We are assuming thata any insertions not in CAT1-5 are non overlaps
        #  and printing them out
#               print ("Comparing", SW_score, "to", SW_score_i, ".", "no overlap")
                        OUTPUT.write(HIT)


        #We have printed the RM hit, but need to add new line to output file
        # or all would be on the same line
        OUTPUT.write("\n")

        #Before continuting with the loop, make sure not to continue all the way so
        # that the last element in the RM file is "q" since there will be no "qi"
        #this is minus two because we cant do a comparison with the last entry, and array starts at entry 0
        if i == len(HIT_ARRAY)-2:
                break

        #iterate to next insertion
        #????????????????????????????
        # how is this iteration step affecting the previous if statement
        i = i + 1


#we are having issues printing the last insertion - Faking it here:
OUTPUT.write("\n")
OUTPUT.write(HIT_ARRAY[i+1])

