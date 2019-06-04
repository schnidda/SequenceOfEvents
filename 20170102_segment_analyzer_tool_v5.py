#new segment analyzer tool [2016 11 10] by schnidda
#version 0.5 2017-01-04
#implemented meltingTemp from Biopython

import sys
import csv

#Biopython
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

#defining several functions for sequence handling
def complement(seq):
  complement = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A'}
  bases =list(seq)
  bases= [complement[base] for base in bases]
  return ''.join(bases)

def reverse_complement(s):
    return complement(s[::-1])

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
        return i+1

#self-made deltaG calculation with NN model (looked up from wikipedia https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics on 2017-01-04)
def calc_deltaG(seq):
    deltaG=0
    actual_position=0

    #values from SantaLucia '98 at 1M Sodium in kJ/mol
    # calculated from deltaG(total)=deltaH-T(K)*deltaS, reference 22
    nearest_neighbours = ['AA','TT','AT','TA', 'CA', 'GT', 'CT', 'GA', 'CG', 'GC', 'GG', 'CC', 'TG', 'AC', 'AG', 'TC']
    neighbour_energies = [-4.26,-4.26,-3.67,-2.50, -6.12, -6.09, -5.40, -5.51, -9.07, -9.36, -7.66,-7.66, -6.12, -6.09, -5.40, -5.51]

    for pair in nearest_neighbours:
        actual_position=nearest_neighbours.index(pair)
        pair_occurence=seq.count(pair)
        pair_energy = neighbour_energies[actual_position]
        deltaG+=(pair_occurence*pair_energy)
    if seq[0] in {'A','T'}:
        deltaG+=4.31
    elif seq[0] in {'C', 'G'}:
        deltaG+=4.05
    if seq[-1] in {'A','T'}:
        deltaG+=4.31
    elif seq[-1] in {'C','G'}:
        deltaG+=4.05

    #convert to kcal/mol
    deltaG*=0.239006
    return deltaG


def end_stacks(seq):
    #values from Kilchherr et al '15 at 20 mM Magnesium
    stacking_neighbours = ['AA','TT','AT','TA', 'CA', 'GT', 'CT', 'GA', 'CG', 'GC', 'GG', 'CC', 'TG', 'AC', 'AG', 'TC']
    stacking_energies = [-1.36, -1.36, -2.35, -1.01, -0.81, -2.03, -1.60, -1.39, -2.06, -3.42, -1.64, -1.64, -0.81, -2.03, -1.60, -1.39]
    start_stack_position=stacking_neighbours.index(seq[0:2])
    start_stack_energy=stacking_energies[start_stack_position]

    end_stack_postion=stacking_neighbours.index(seq[-2:])
    end_stack_energy=stacking_energies[end_stack_postion]

    return (start_stack_energy, end_stack_energy)

#needs one file with scaffold sequence, one with staples + segments from gillespie lookup and the third one with only staples to find self complements
filename = sys.argv[1]
filename2 = sys.argv[2]
fp= open(sys.argv[1])
fp2 = open(sys.argv[2])
fp3 = open(sys.argv[3])

#necessary to strip the brackets [] and ''!!!
scaff_sequ=str(fp.readlines())
scaff_sequ=scaff_sequ[2:-2]
scaff_sequ=scaff_sequ.upper()

length_of_segments1=7

comp_oligos_sequ=[]
levelpoints_5prime=[]
levelpoints_3prime=[]
for line in csv.reader(fp3,dialect="excel-tab"):
    oligo_self_comp=line[0]
    print oligo_self_comp
    levelpoint5=line[1]
    levelpoint3=line[2]
    comp_oligos_sequ.append(oligo_self_comp)
    levelpoints_5prime.append(levelpoint5)
    levelpoints_3prime.append(levelpoint3)

print "5prime: ", levelpoints_5prime
print "3prime: ", levelpoints_3prime
number_of_oligos=file_len(filename2)

#changed to write mode from append
result1=open("Segment_analysis_length.txt","w")
#result1=open("Segment_analysis_length"+str(length_of_segments1)+".txt","w")

result_all=open("Segment_oligos_results.txt", "w")
#result_all=open("Segment_oligos_results"+str(length_of_segments1)+".txt", "w")

result_2segments=open("Segment_oligos_results_2segments.txt", "w")


oligo_index=0
segment_index=0
num_comp_occur=0
num_comp_occur1=0
num_comp_occur2=0
num_comp_occur3=0
maximum_segment_length=35
scaffold_length=7560
fiveprime=[0,1]
BaseReduction=1
PreLastSegment=""
LastSegment=""
IndexOfLine=-1
#need some description here!
result1.write("number   length  sequence    complement  position    occurence   Tm  deltaG  binding_time\n")
result_2segments.write("pos_5_1 seq_5_1 Tm_5_1 Tm1_5_1 Tm2_5_1 Tm3_5_1 deltaG_5_1 deltaG_5_1_5_1 deltaG_5_1_3_1 deltaG_5_1_5_3_1 add_stack_5_1_5 add_stack_5_1_3 length_5_1 occurence_5_1 occurence_comp_5_1 occurence_5_1_1 occurence_comp_5_1_1 occurence_5_1_2 occurence_comp_5_1_2 occurence_5_1_3 occurence_comp_5_1_3 pos_5_2 seq_5_2 Tm_5_2 Tm1_5_2 Tm2_5_2 Tm3_5_2 deltaG_5_2 deltaG_5_2_5_1 deltaG_5_2_3_1 deltaG_5_2_5_3_1 add_stack_5_2_5 add_stack_5_2_3 length_5_2 occurence_5_2 occurence_comp_5_2 occurence_5_2_1 occurence_comp_5_2_1 occurence_5_2_2 occurence_comp_5_2_2 occurence_5_2_3 occurence_comp_5_2_3 pos_3_2 seq_3_2 Tm_3_2 Tm1_3_2 Tm2_3_2 Tm3_3_2 deltaG_3_2 deltaG_3_2_5_1 deltaG_3_2_3_1 deltaG_3_2_5_3_1 add_stack_3_2_5 add_stack_3_2_3 length_3_2 occurence_3_2 occurence_comp_3_2 occurence_3_2_1 occurence_comp_3_2_1 occurence_3_2_2 occurence_comp_3_2_2 occurence_3_2_3 occurence_comp_3_2_3 pos_3_1 seq_3_1 Tm_3_1 Tm1_3_1 Tm2_3_1 Tm3_3_1 deltaG_3_1 deltaG_3_1_5_1 deltaG_3_1_3_1 deltaG_3_1_5_3_1 add_stack_3_1_5 add_stack_3_1_3 length_3_1 occurence_3_1 occurence_comp_3_1 occurence_3_1_1 occurence_comp_3_1_1 occurence_3_1_2 occurence_comp_3_1_2 occurence_3_1_3 occurence_comp_3_1_3 ")
for line in csv.reader(fp2, dialect="excel-tab"):
    length_of_actual_segment=int(line[1])
    print "--new oligo--"
    print line
    if length_of_actual_segment<maximum_segment_length:
        oligo_sequ1=line[2]#.rstrip()
        print oligo_sequ1
        rev_comp = reverse_complement(oligo_sequ1)
        rev_comp1 = reverse_complement(oligo_sequ1[BaseReduction:])
        rev_comp2 = reverse_complement(oligo_sequ1[:-BaseReduction])
        rev_comp3 = reverse_complement(oligo_sequ1[BaseReduction:-BaseReduction])

        #Tm calculations
        if(len(rev_comp3)>1):
            mySeq=Seq(rev_comp)
            mySeq1=Seq(rev_comp1)
            mySeq2=Seq(rev_comp2)
            mySeq3=Seq(rev_comp3)

            meltingT='%0.2f' %mt.Tm_NN(mySeq, Na=5, Tris=5, Mg=20, saltcorr=7)
            meltingT1='%0.2f' %mt.Tm_NN(mySeq1, Na=5, Tris=5, Mg=20, saltcorr=7)
            meltingT2='%0.2f' %mt.Tm_NN(mySeq2, Na=5, Tris=5, Mg=20, saltcorr=7)
            meltingT3='%0.2f' %mt.Tm_NN(mySeq3, Na=5, Tris=5, Mg=20, saltcorr=7)
        else:
            meltingT=0
            meltingT1=0
            meltingT2=0
            meltingT3=0

        index=scaff_sequ.find(rev_comp)
        index1=scaff_sequ.find(rev_comp1)
        index2=scaff_sequ.find(rev_comp2)
        index3=scaff_sequ.find(rev_comp3)

        print "--new segment--"
        print index, scaff_sequ[index-1:index+1+length_of_actual_segment]
        print rev_comp, rev_comp1, rev_comp2

        result1.write("\n"+"#"+str(line[0])+"    "+str(line[1])+"    "+oligo_sequ1+" "+str(rev_comp)+"   "+str(index)+"  ")

        actual_energy=calc_deltaG(rev_comp)
        if((0<index) and (index+length_of_actual_segment<scaffold_length)):
            actual_energy2=calc_deltaG(scaff_sequ[index-1:index+length_of_actual_segment])
            actual_energy3=calc_deltaG(scaff_sequ[index:index+1+length_of_actual_segment])
            actual_energy4=calc_deltaG(scaff_sequ[index-1:index+1+length_of_actual_segment])
            intermediate_seq=scaff_sequ[index-1:index+1+length_of_actual_segment]
            #intermediate_seq=scaff_sequ[-1]+scaff_sequ[index+length_of_actual_segment+1]

            #control prints
            #print scaff_sequ
            #print len(scaff_sequ)
            #print intermediate_seq, intermediate_seq[0:2], intermediate_seq[-2:]
            additional_stacking=end_stacks(intermediate_seq)

        elif(index==0):
            intermediate_seq=scaff_sequ[-1]+scaff_sequ[index+length_of_actual_segment]
            actual_energy2=calc_deltaG(intermediate_seq)
            actual_energy3=calc_deltaG(scaff_sequ[index:index+1+length_of_actual_segment])
            intermediate_seq=scaff_sequ[-1]+scaff_sequ[index+length_of_actual_segment+1]
            actual_energy4=calc_deltaG(intermediate_seq)
            additional_stacking=end_stacks(intermediate_seq)

        elif(index+length_of_actual_segment>=scaffold_length):
            intermediate_seq=scaff_sequ[index-1:scaffold_length]+scaff_sequ[0:((index+length_of_actual_segment)%scaffold_length)]
            actual_energy2=calc_deltaG(intermediate_seq)
            intermediate_seq=scaff_sequ[index:scaffold_length]+scaff_sequ[0:((index+length_of_actual_segment+1)%scaffold_length)]
            actual_energy3=calc_deltaG(intermediate_seq)
            intermediate_seq=scaff_sequ[index-1:scaffold_length]+scaff_sequ[0:((index+length_of_actual_segment+1)%scaffold_length)]
            actual_energy4=calc_deltaG(intermediate_seq)
            additional_stacking=end_stacks(intermediate_seq)

        print line[0], line[1]
        print actual_energy, actual_energy2, actual_energy3, actual_energy4
        num_occur=scaff_sequ.count(rev_comp)
        num_occur1=scaff_sequ.count(rev_comp1)
        num_occur2=scaff_sequ.count(rev_comp2)
        num_occur3=scaff_sequ.count(rev_comp3)

        for element in comp_oligos_sequ:
            num_comp_occur+=element.count(rev_comp)
        for element in comp_oligos_sequ:
            num_comp_occur1+=element.count(rev_comp1)
        for element in comp_oligos_sequ:
            num_comp_occur2+=element.count(rev_comp2)
        for element in comp_oligos_sequ:
            num_comp_occur3+=element.count(rev_comp3)

        #result1.write(str(num_occur)+"  "+str(meltingT)+"  "+str(actual_energy)+"\n")
        result_all.write(str(line[1])+" "+str(num_occur)+"  "+str(num_comp_occur)+" "+str(num_occur1)+"  "+str(num_comp_occur1)+" "+str(num_occur2)+"  "+str(num_comp_occur2)+" "+str(num_occur3)+"  "+str(num_comp_occur3))

        if segment_index in fiveprime:
            result_2segments.write(str(index)+"  "+str(oligo_sequ1)+"   "+str(meltingT)+"  "+str(meltingT1)+"  "+str(meltingT2)+"  "+str(meltingT3)+"  "+str(actual_energy)+" "+str(actual_energy2)+" "+str(actual_energy3)+" "+str(actual_energy4)+"   "+str(additional_stacking[0])+" "+str(additional_stacking[1])+" "+str(line[1])+" "+str(num_occur)+"  "+str(num_comp_occur)+" "+str(num_occur1)+"  "+str(num_comp_occur1)+" "+str(num_occur2)+"  "+str(num_comp_occur2)+"   "+str(num_occur3)+"  "+str(num_comp_occur3)+"   ")
        else:
            PreLastSegment=LastSegment
            LastSegment=str(str(index)+"  "+str(oligo_sequ1)+"   "+str(meltingT)+"  "+str(meltingT1)+"  "+str(meltingT2)+"  "+str(meltingT3)+"  "+str(actual_energy)+" "+str(actual_energy2)+" "+str(actual_energy3)+" "+str(actual_energy4)+"   "+str(additional_stacking[0])+" "+str(additional_stacking[1])+" "+str(line[1])+" "+str(num_occur)+"  "+str(num_comp_occur)+" "+str(num_occur1)+"  "+str(num_comp_occur1)+" "+str(num_occur2)+"  "+str(num_comp_occur2)+"   "+str(num_occur3)+"  "+str(num_comp_occur3)+"   ")

        if segment_index==0:
            print "Index", IndexOfLine
            result1.write(str(num_occur)+" "+str(meltingT)+"   "+str(actual_energy)+"  "+str(levelpoints_5prime[IndexOfLine]))
        else:
            result1.write(str(num_occur)+" "+str(meltingT)+"   "+str(actual_energy)+"      ")

        num_comp_occur=0
        num_comp_occur1=0
        num_comp_occur2=0
        num_comp_occur3=0

        result_all.write("\n")

        segment_index+=1


    else:
        if IndexOfLine==-1:
            result1.write("#"+str(line[0])+"    "+str(line[1])+"    "+str(line[2]))
        else:
            result1.write(str(levelpoints_3prime[IndexOfLine])+"\n"+"#"+str(line[0])+"    "+str(line[1])+"    "+str(line[2]))
        result_2segments.write(PreLastSegment+LastSegment+"\n")
        #result_all.write("\n")
        segment_index=0
        IndexOfLine+=1

if IndexOfLine==139:
    result1.write(str(levelpoints_3prime[IndexOfLine]))
    oligo_index+=1

#adding the last segments of the last oligo manually!
result_2segments.write(PreLastSegment+LastSegment+"\n")
