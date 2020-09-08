#!/usr/local/bin/python3
import sys
import Bio
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# Class that looks through all of the permutations of the input
# sequence that keeps the overall codon output the same. It then
# looks through a whitespace separated file on enzymes and their
# sequence to find possible enzyme sites, given a desired scrolling
# window length.
class EnzymeChooser:
    def __init__(self, inputSeq, enzymeFile, windowLength):
        # dictionary of codons to amino acids
        self.fwdTable = CodonTable.standard_dna_table.forward_table
        self.inputSeq = inputSeq
        self.enzymeFile = enzymeFile
        self.windowLength = windowLength
        # dictionary of the enzymes we get from the file
        self.enzymes = {}
        # list of lists with the possible codons for each codon
        self.permutationsList = []
        # output dict with enzymes as the keys, then the value is
        # a tuple: (replacement sequence, enzyme starting index)
        self.validEnzymes = {}

    # make a dictionary of enzymes and their sequences
    def __makeEnzymesDict(self):
        enzymeLineList = open(self.enzymeFile, "r").read().splitlines()
        for line in enzymeLineList:
            splitLine = line.split()
            # key: enzyme, value: sequence
            self.enzymes[splitLine[0]] = splitLine[1]

    # make the list of all codons that could be interchanged
    # at each codon place while keeping the overall codons the same
    def __makePermuationsList(self):
        codons = Seq.translate(Seq(input, generic_dna))
        for c in codons:
            self.permutationsList.append( \
                [k for k,v in self.fwdTable.items() if v == c])

    # compares two sequences using the compare dictionary
    #def compareSequences(self, s1, s2):

    # searching for enzymes using a scrolling window of windowLength
    # returns a list of enzymes and the location that they start at
    # in the sequence
    def __searchForEnzymes(self, sequence):
        windowStart = 0
        windowEnd = self.windowLength
        if (windowEnd > len(sequence)):
            windowEnd = len(sequence)
        while (windowEnd <= len(sequence)):
            print("HI")

    # recursive function that goes through possible permutations
    # and stores the previous permutations in previousIndicies
    def __permutate(self, previousIndicies, currentIndex):
        # We have completed a permuation
        if (currentIndex == len(self.permutationsList)):
            newSequence = ""
            # make the new complete permutation based on the previous indicies
            for i in range(len(previousIndicies)):
                newSequence += self.permutationsList[i][previousIndicies[i]]
            enzymesListAndLocations = searchForEnzymes(newSequence)
            for enzymeAndLocation in enzymesList:
                # check if we have already found this enzyme so far
                if enzymeAndLocation[0] in self.validEnzymes:
                    self.validEnzymes[enzymeAndLocation[0]].append( \
                        (newSequence, enzymeAndLocation[1])
                # if we haven't make a new entry
                else:
                    self.validEnzymes[enzymeAndLocation[0]] = \
                        [(newSequence, enzymeAndLocation[1])]
        # Keep going with different options for currentIndex
        else:
            for i in range(len(list[currentIndex])):
                self.permutate(list, previousIndicies + [i], currentIndex + 1)

    # only public function, runs all the other functions to
    # get the final product
    def getEnzymes(self):
        self.__makeEnzymesDict()
        self.__makePermuationsList()



input = "GATCGATGGGCCTATATAGGATCGAAAATC"
enzymesFile = "enzymes.txt"
windowLength = 15
chooser = EnzymeChooser(input, enzymesFile, windowLength)
validEnzymes = chooser.getEnzymes()
print(validEnzymes)
