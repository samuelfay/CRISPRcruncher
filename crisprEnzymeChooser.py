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
    # dictionary for the nucleotide ambiguity codes
    nucleotideEqualities = {
        "N": ["A", "C", "G", "T"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "W": ["A", "T"],
        "S": ["C", "G"],
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "M": ["A", "C"],
        "K": ["G", "T"],
    }
    def __init__(self, inputSeq, enzymeFile, windowLength):
        # dictionary of codons to amino acids
        self.fwdTable = CodonTable.standard_dna_table.forward_table
        self.inputSeq = inputSeq
        self.enzymeFile = enzymeFile
        self.windowLength = windowLength
        # list of lists with the possible codons for each codon
        self.permutationsList = []
        # output dict with enzymes as the keys, then the value is
        # a tuple: (binding sequence, shortened binding sequence, least-
        # different input sequence permutation (defaults to empty string),
        # enzyme starting index) shortened means the leading and trailing
        # N's have been stripped from the sequence.
        # if the starting index is -1, then the enzyme hasn't been found
        self.enzymes = {}

    # make a dictionary of enzymes and their sequences
    def __makeEnzymesDict(self):
        enzymeLineList = open(self.enzymeFile, "r").read().splitlines()
        for line in enzymeLineList:
            splitLine = line.split()
            # key: enzyme, value: sequence
            self.enzymes[splitLine[0]] = (splitLine[1], \
                splitLine[1].split("N"), "", -1)

    # make the list of all codons that could be interchanged
    # at each codon place while keeping the overall codons the same
    def __makePermuationsList(self):
        codons = Seq.translate(Seq(input, generic_dna))
        for c in codons:
            self.permutationsList.append( \
                [k for k,v in self.fwdTable.items() if v == c])

    # compares two sequences using the compare dictionary
    def __compareSequences(self, s1, s2):
        if len(s1) != len(s2):
            return False
        for i in range(len(s1)):
            if s1[i] == s2[i] or s2[i] in nucleotideEqualities[s1[i]]:
                continue
            else:
                return False
        return True

    # since we want to make the permuation that does the least
    # work for each enzyme, we want the sequence that is closest
    # to the original sequence
    def __calculateSequenceDifferences(self, sequence):
        differences = 0
        for i in range(len(sequence)):
            if sequence[i] != self.inputSeq:
                differences += 1
        return differences

    # searching for enzymes using a scrolling window of windowLength
    # returns a list of enzymes and the location that they start at
    # in the sequence
    def __searchForEnzymes(self, sequence):
        windowStart = 0
        windowEnd = self.windowLength
        sequenceDiff = self.__calculateSequenceDifferences(sequence)
        if (windowEnd > len(sequence)):
            windowEnd = len(sequence)
        # while the moving window doesn't hit the end
        while (windowEnd <= len(sequence)):
            # get the current base pairs in the window
            slice = sequence[windowStart:windowEnd]
            # check for every enzyme starting at every character in slice
            for i in range(len(slice)):
                lengthLeft = self.windowLength - i
                # iterate through every enzyme
                for key in self.enzymes:
                    currentEnzyme = self.enzymes[key]
                    # if we have already found this enzyme starting on this
                    # character in the sequence then we don't need to find it
                    # again
                    if currentEnzyme[3] == i + windowStart:
                        continue
                    enzymeLength = len(currentEnzyme[1])
                    # if the enzyme won't fit, we move on to the next
                    if enzymeLength > lengthLeft:
                        continue
                    if self.__compareSequences(slice[i:i+enzymeLength]):
                        # if we haven't found one yet or if the current sequence
                        # is closer to the original, then we replace it
                        if currentEnzyme[2] == "" or \
                            self.__calculateSequenceDifferences( \
                            currentEnzyme[2]) > sequenceDiff:
                            currentEnzyme[2] = sequence
                            currentEnzyme[3] = i + windowStart
            windowStart += 1
            windowEnd += 1


    # recursive function that goes through possible permutations
    # and stores the previous permutations in previousIndicies
    def __permutateAndSearch(self, previousIndicies, currentIndex):
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
