#!/usr/bin/python3

import re 
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", metavar = "FILE", help = "Please enter the file (with extension) that you wish to parse", default = "orfs.fa")
	parser.add_argument("-o", "--output", metavar = "FILE", help = "This shall provide an output for observation provided the given input")
	parser.add_argument("-c", "--choice", type=str, help = "Choose which enzyme to use out of Trypsin, Endoproteinase_LysC, Endoproteinase_Arg and V8_Proteinase", choices = ['Trypsin', 'Endoproteinase_LysC', 'Endoproteinase_Arg', 'V8_Proteinase'])

args = parser.parse_args() 

def read_fasta(file):
    """Read and store the contents of files"""
# The point of this parser is to provide a function which allows for the storing of names and sequences provided
# by the appropriate 'name' and 'seq' variables. This is perpetuated by the presumption that the 'name' will be
# the title, or name of the sequence to be parsed, as indicated by the '>' symbol, indicating that it is a FASTA
# sequence. Should '>' not be the first character of a line, it will automatically be stored as sequence data.
# This function relies on the later definition of the variable 'file.' This should be an input provided by what
# the user would otherwise specify.

name, seq = None, [] #To store the sequences in a list for line in file:
line = line.rstrip() #To read the lines
if re.search(r">", line): #to read the FASTA formatted line
	if name : yield (name, ''.join(seq))
	name, seq = line, [] 
else:
	seq.append(line) #To read the sequence line if name: yield(name, "\n".join(seq))


def Trypsin(sequence): 
"""Cleave where K or R are present, except where either is adjacent to P"""

	cleavage = [];
	proteinNumber = 0;
	Lysine_K = False; #Supposing that Lysine is not the previous letter, it should continue to read the sequence
	Arginine_R = False; #Supposing Arginine was not the previous letter, it would continue to read the sequence
	concatenatedSequences = []

	#Loop through every letter of the sequence
	for protein in sequence:

		if ((protein != 'P') & (Lysine_K == True) or (Arginine_R == True)):
			cleavage.append(proteinNumber - 1);

		if protein == 'K':
			Lysine_K = True; #If the programme finds K, it will cleave before K
		else:
			Lysine_K = False; #If K has not been found, it will continue to iterate 
		if protein == 'R':
			Arginine_R = True; #If R has been found, there should be cleavage
		else:
			Arginine_R = False; #If R has not been found, it will continue to iterate 
		proteinNumber += 1

	lastCut = 0;
	for point in cleavage:

		concatenatedSequences.append(sequence[lastCut:point]); #This should concatenate the sequences prior to the cleavage

		lastCut = point;

	return concatenatedSequences

def Endoproteinase_LysC(sequence): 
"""Cleave where K is present, except when K is adjacent to P"""

	cutPoints = [];
	proteinNumber = 0;
	Lysine_K = False; #Supposing that Lysine is not the previous letter, it should continue to read the sequence
	concatenatedSequences = [];

	#Loop through every letter of the sequence
	for protein in sequence:

		if ((protein != 'P') & (Lysine_K == True)):
			cutPoints.append(proteinNumber - 1);

		if protein == 'K':
			Lysine_K = True; #If the programme finds K, it will cleave before K
		else:
			Lysine_K = False; #If K has not been found, it will continue to iterate over the sequence

		proteinNumber += 1

	lastCut = 0;
	for point in cutPoints:

		concatenatedSequences.append(sequence[lastCut:point]); #This should concatenate the sequences prior to the cleavage

		lastCut = point;

	return concatenatedSequences

def Endoproteinase_Arg(sequence): 
"""Cleave where R is present, except when it is adjacent to P"""

	cutPoints = [];
	proteinNumber = 0;
	Arginine_R = False; #Supposing that Arginine is not the previous letter, it should continue to read the sequence
	concatenatedSequences = [];

	#Loop through every letter of the sequence
	for protein in sequence:

		if ((protein != 'P') & (Arginine_R == True)):
			cutPoints.append(proteinNumber - 1);

		if protein == 'R':
			Arginine_R = True; #If the programme finds R, it will cleave before R
		else:
			Arginine_R = False; #If R has not been found, it will continue to iterate over the sequence

		proteinNumber += 1

	lastCut = 0;
	for point in cutPoints:

		concatenatedSequences.append(sequence[lastCut:point]); #This should concatenate the sequences prior to the cleavage

		lastCut = point;

	return concatenatedSequences

def V8_Proteinase(sequence):
"""Cleave where E is present, except when adjacent to P."""

	cutPoints = [];
	proteinNumber = 0;
	Glutamine_E = False; #Supposing that Glutamine is not the previous letter, it should continue to read the sequence
	concatenatedSequences = [];

	#Loop through every letter of the sequence
	for protein in sequence:

		if ((protein != 'P') & (Glutamine_E == True)):
			cutPoints.append(proteinNumber - 1);

		if protein == 'E':
			Glutamine_E = True; #If the programme finds E, it will cleave before E
		else:
			Glutamine_E = False; #If E has not been found, it will continue to iterate over the sequence

		proteinNumber += 1

	lastCut = 0;
	for point in cutPoints:

		concatenatedSequences.append(sequence[lastCut:point]); #This should concatenate the sequences prior to the cleavage

		lastCut = point;

	return concatenatedSequences


enzymes = {Trypsin, Endoproteinase_LysC, Endoproteinase_Arg, V8_Proteinase} #Create a dictionary containing each of the enzymes

with open(args.input) as file:
	for name, seq in read_fasta(file): + "\n" 

		if args.choice == "Trypsin":
			for character in range(len(Trypsin(seq))):
				t_output = (name + "\n" + (Trypsin(seq)[character])
				output = open(args.output, "a+") 
				output.write(t_output)
                    
		if args.choice == "Endoproteinase_LysC":
			for character in range(len(Endoproteinase_LysC(seq))):
				c_output = (name + "\n" + (Endoproteinase_LysC(seq)[character]) + "\n")
				ec_output = open(args.output, "a+") 
				ec_output.write(c_output)
                    
		if args.choice == "Endoproteinase_Arg":
			for character in range(len(Endoproteinase_Arg(seq))):
				arg_output = (name + "\n" + (Endoproteinase_Arg(seq)[character]) + "\n")
				ea_output = open(args.output, "a+") 
				ea_output.write(arg_output)
                    
		if args.choice == "V8_Proteinase":
			for character in range(len(V8_Proteinase(seq))):
				v_output = (name + "\n" + (V8_Proteinase(seq)[character]) + "\n")
				v8_output = open(args.output, "a+") 
				v8_output.write(v_output)
                            
                            
                            
