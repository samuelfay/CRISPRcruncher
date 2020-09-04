#!/usr/local/bin/python3
import sys
import Bio
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


codonTable = CodonTable.standard_dna_table
fwdTable = codonTable.forward_table

translate = Seq.translate

def getBasePairPossibilities(codons):
    l = []
    for c in codons:
        l.append([k for k,v in fwdTable.items() if v == c])
    return l

def compareSequences(s1, s2):

def searchForEnzymes(sequence):

def permutate(list, previousIndicies, currentIndex, outputDict):


enzymeLineList = open("enzymes.txt", "r").read().splitlines()
enzymes = {}
for line in enzymeLineList:
    splitLine = line.split("\t")
    enzymes[splitLine[1]] = splitLine[0]



input = "GATCGATGGGCCTATATAGGATCGAAAATC"
desiredReplacement = Seq(input, generic_dna)
print(codonTable)

codons = desiredReplacement.translate()
print(codons)

permutationsList = getBasePairPossibilities(codons)
print(permutationsList)

validEnzymes = {}
permutate(permutationsList, [], 0, validEnzymes)
