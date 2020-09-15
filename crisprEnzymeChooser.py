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
        # list for current permutations
        self.permutationsList = []
        # output dict with enzymes as the keys, then the value is
        # a list: [binding sequence, shortened binding sequence, least-
        # different input sequence window permutation (defaults to empty string),
        # enzyme starting index, overall best sequence with length of the input,
        # isReversed]
        # shortened means the leading and trailing
        # N's have been stripped from the sequence.
        # if the starting index is -1, then the enzyme hasn't been found
        self.enzymes = {}
        # variables for the start and end of the window
        self.windowStart = 0
        self.windowEnd = 0

    # make a dictionary of enzymes and their sequences
    def __makeEnzymesDict(self):
        enzymeLineList = open(self.enzymeFile, "r").read().splitlines()
        for line in enzymeLineList:
            splitLine = line.split()
            # key: enzyme, value: sequence
            self.enzymes[splitLine[0]] = [splitLine[1], \
                splitLine[1].strip("N"), "", -1, "", 0]

    # make the list of all codons that could be interchanged
    # at each codon place while keeping the overall codons the same
    def __makePermuationsList(self, slice):
        # clear list
        self.permutationsList = []
        codons = Seq.translate(Seq(slice, generic_dna))
        #print(codons)
        for c in codons:
            self.permutationsList.append( \
                [k for k,v in self.fwdTable.items() if v == c])
        #print(self.permutationsList)

    # compares two sequences using the compare dictionary
    # enzyme seq can have ambiguities, dna seq can't
    def __compareSequences(self, dnaSeq, enzymeSeq):
        if len(dnaSeq) != len(enzymeSeq):
            return False
        for i in range(len(dnaSeq)):
            if dnaSeq[i] == enzymeSeq[i] or (enzymeSeq[i] in \
                EnzymeChooser.nucleotideEqualities and dnaSeq[i] in \
                EnzymeChooser.nucleotideEqualities[enzymeSeq[i]]):
                continue
            else:
                return False
        return True

    # since we want to make the permuation that does the least
    # work for each enzyme, we want the sequence that is closest
    # to the original sequence
    def __calculateSequenceDifferences(self, sequence):
        differences = 0
        for i in range(self.windowLength):
            if sequence[i] != self.inputSeq[i + self.windowStart]:
                differences += 1
        return differences


    # searching for enzymes using a scrolling window of windowLength
    # returns a list of enzymes and the location that they start at
    # in the sequence
    def __searchForEnzymes(self, sequence):
        sequenceDiff = self.__calculateSequenceDifferences(sequence)
        # check for every enzyme starting at every character in slice
        seqObj = Seq(sequence)
        reverse = str(seqObj.reverse_complement())
        #print(sequence)
        #print(reverse)
        for i in range(len(sequence)):
            lengthLeft = self.windowLength - i
            # iterate through every enzyme
            for key in self.enzymes:
                currentEnzyme = self.enzymes[key]
                enzymeLength = len(currentEnzyme[1])
                # if the enzyme won't fit, we move on to the next
                if enzymeLength > lengthLeft:
                    continue
                if self.__compareSequences(sequence[i:i+enzymeLength], currentEnzyme[1]):
                    # if we haven't found one yet or if the current sequence
                    # is closer to the original, then we replace it
                    if currentEnzyme[2] == "" or \
                        self.__calculateSequenceDifferences( \
                        currentEnzyme[2]) > sequenceDiff:
                        # we found one!
                        currentEnzyme[2] = sequence
                        currentEnzyme[3] = i + self.windowStart
                        currentEnzyme[4] = self.inputSeq[0:self.windowStart] + \
                            sequence + self.inputSeq[self.windowEnd:len(self.inputSeq)]
                        currentEnzyme[5] = 0
                if self.__compareSequences(reverse[i:i+enzymeLength], currentEnzyme[1]):
                    # if we haven't found one yet or if the current sequence
                    # is closer to the original, then we replace it
                    if currentEnzyme[2] == "" or \
                        self.__calculateSequenceDifferences( \
                        currentEnzyme[2]) > sequenceDiff:
                        # we found one, but it's on the other side!
                        currentEnzyme[2] = sequence
                        currentEnzyme[3] = self.windowEnd - i - 1
                        currentEnzyme[4] = self.inputSeq[0:self.windowStart] + \
                            sequence + self.inputSeq[self.windowEnd:len(self.inputSeq)]
                        currentEnzyme[5] = 1

    # recursive function that goes through possible permutations
    # of each scrolling window
    # and stores the previous permutations in previousIndicies
    def __permutate(self, previousIndicies, currentIndex):
        # We have completed a permuation
        if (currentIndex == len(self.permutationsList)):
            newSequence = ""
            # make the new complete permutation based on the previous indicies
            for i in range(len(previousIndicies)):
                newSequence += self.permutationsList[i][previousIndicies[i]]
            #print("Current sequence:" + newSequence)
            #print("Current codons: " + Seq.translate(Seq(newSequence, generic_dna)))

            self.__searchForEnzymes(newSequence)
        # Keep going with different options for currentIndex
        else:
            for i in range(len(self.permutationsList[currentIndex])):
                self.__permutate(previousIndicies + [i], currentIndex + 1)

    # only public function, runs all the other functions to
    # get the final product
    def getEnzymes(self):
        self.__makeEnzymesDict()
        self.windowStart = 0
        self.windowEnd = self.windowLength
        # while the moving window doesn't hit the end
        while (self.windowEnd <= len(self.inputSeq)):
            print("SEARCHING WINDOW FROM " + str(self.windowStart) + " TO " \
                + str(self.windowEnd))
            # get the current base pairs in the window
            slice = self.inputSeq[self.windowStart:self.windowEnd]
            #print(slice)
            self.__makePermuationsList(slice)
            self.__permutate([], 0)
            self.windowStart += 3
            self.windowEnd += 3
        return self.enzymes



input = "GATCGATGGGCCTATATAGGATCGAAAATC"
enzymesFile = "enzymes.txt"
windowLength = 15
chooser = EnzymeChooser(input, enzymesFile, windowLength)
enzymes = chooser.getEnzymes()
found = 0
for key in enzymes:
    if enzymes[key][3] != -1 and len(enzymes[key][1]) >= 6:
        found += 1
        print(key)
        print(enzymes[key])
print("Found " + str(found) + " enzymes.")
