#!/usr/local/bin/python3
import sys
import Bio
from Bio.Data import CodonTable
from Bio.Seq import Seq
from flask import Flask, render_template, request, url_for, flash, redirect
from flask import Markup, Response

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
    def __init__(self, inputSeq, enzymeFile, windowLength, minEnzymeLength):
        # dictionary of codons to amino acids
        self.fwdTable = CodonTable.standard_dna_table.forward_table
        self.inputSeq = inputSeq
        self.enzymeFile = enzymeFile
        self.windowLength = windowLength
        self.numWindows = (len(self.inputSeq) - self.windowLength) / 3
        self.minEnzymeLength = minEnzymeLength
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
                splitLine[1].strip("N"), "", -1, "", False]

    # make the list of all codons that could be interchanged
    # at each codon place while keeping the overall codons the same
    def __makePermuationsList(self, slice):
        # clear list
        self.permutationsList = []
        codons = Seq.translate(Seq(slice))
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
                # get the length without the N'
                trueLength = self.__enzymeLength(currentEnzyme[1])
                # if the enzyme won't fit, or, its shorter than the min length
                # we move on to the next
                if enzymeLength > lengthLeft or trueLength < self.minEnzymeLength:
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
                        currentEnzyme[5] = False
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
                        currentEnzyme[5] = True

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

    # helper function for getting startIndex
    # used for list sorting in the public funciton
    def __getStartIndex(self, enzyme):
        if enzyme[5]:
            return enzyme[3] - len(enzyme[1])
        else:
            return enzyme[3]

    # marks up the sequence for output in the table.
    # colors the changed sequence and makes the enzyme site bolded
    # returns a new version of the input sequence with html in it of the Markup class.
    def __markup(self, sequence, startIndex, enzymeLength, isReversed):
        out = ""
        # modify the start index if it's reversed
        if isReversed:
            startIndex = startIndex - enzymeLength + 1
        for i in range(len(self.inputSeq)):
            # if this is the start of the enzyme sites
            if i == startIndex:
                out += "<b>"
            # if we are past the enzyme site
            elif i == startIndex + enzymeLength:
                out += "</b>"
            # if this nucleotide was changed
            if self.inputSeq[i] != sequence[i]:
                out += "<span style=\"color:#e93c73\">" + sequence[i] + "</span>"
            # otherwise just add the nucleotide
            else:
                out += sequence[i]
        return Markup(out)

    # get true length of enzyme (without "N"'s')
    def __enzymeLength(self, enzyme):
        return len(enzyme) - enzyme.count("N")

    def __isPalindrome(self, sequence):
        print(sequence)
        if sequence == Seq(sequence).reverse_complement():
            return True
        return False
    # only public function, runs all the other functions to
    # get the final product
    def getEnzymes(self, outList):
        self.__makeEnzymesDict()
        self.windowStart = 0
        self.windowEnd = self.windowLength
        currentWindow = 0
        # while the moving window doesn't hit the end
        while (self.windowEnd <= len(self.inputSeq)):
            flash(str(currentWindow) + "/" + str(self.numWindows))
            currentWindow += 1
            # get the current base pairs in the window
            slice = self.inputSeq[self.windowStart:self.windowEnd]
            #print(slice)
            self.__makePermuationsList(slice)
            self.__permutate([], 0)
            self.windowStart += 3
            self.windowEnd += 3
        tempDict = self.enzymes
        #print(tempDict)
        self.enzymes = {}
        # remove enzymes that weren't found, and enzymes below the min length
        for key in tempDict:
            if (tempDict[key][3] != -1 and
            self.__enzymeLength(tempDict[key][1]) >= self.minEnzymeLength):
                self.enzymes[key] = tempDict[key]
        # insert underline and bolded (html format) for the enzyme site
        for key in self.enzymes:
            e = self.enzymes[key]
            # get the site to see if the sequence is a __isPalindrome
            enzymeSite = e[4][e[3]:e[3]+len(e[1])]
            isPalindrome = self.__isPalindrome(enzymeSite)
            # markup the sequence in html
            e[4] = self.__markup(e[4], e[3], len(e[1]), e[5])
            outList.append(self.enzymes[key] + [key] + [isPalindrome])
        outList.sort(key = lambda x:self.__getStartIndex(x))

# function to write the user data to a text files
def writeUserData(name, affiliation, email, organism):
    out = open("../usersInfo.txt", "a")
    out.write("\n" + name + "\t" + affiliation + "\t" + email)
    out.write("\t" + organism)
    out.close()

# function to write the user feedback to a text file
def writeFeedback(feedback):
    out = open("../feedback.txt", "a")
    out.write(feedback + "\n\n\n")
    out.close()

    # input = "GATCGATGGGCCTATATAGGATCGAAAATC"
enzymesFile = "enzymes.txt"
windowLength = 15
# chooser = EnzymeChooser(input, enzymesFile, windowLength)
# enzymes = chooser.getEnzymes()
# found = 0
# for key in enzymes:
#     if enzymes[key][3] != -1 and len(enzymes[key][1]) >= 6:
#         found += 1
#         print(key)
#         print(enzymes[key])
# print("Found " + str(found) + " enzymes.")

app = Flask(__name__)
app.config['SECRET_KEY'] = 'POGGERS'

@app.route("/", methods=("GET", "POST"))
def index():
    # if they submitted the form
    if request.method == 'POST':
        name = request.form['name']
        affiliation = request.form['affiliation']
        email = request.form['email']
        organism = request.form['organism']
        sequence = request.form['sequence'].upper()
        minLength = request.form['minLength']
        errorMessage = ""
        # input validation
        if (len(name) == 0 or len(affiliation) == 0 or len(organism) == 0
                or len(sequence) == 0):
            errorMessage += "Name, Affiliation, Organism, and Sequence are " \
                "required.\n"
        if len(sequence) >= 100 or len(sequence) < 15:
            errorMessage += "The sequence must be from 15 to 99 base pairs.\n"
        if len(sequence) % 3 != 0:
            errorMessage += "Please make sure sequence is in frame and has a length" \
            " that is a multiple of 3.\n"
        for char in sequence:
            if char != "A" and char != "T" and char != "C" and char != "G":
                errorMessage += "The sequence must only contain the base " \
                    "pairs A, T, C, and G.\n"
                break
        try:
            int(minLength)
        except ValueError:
            errorMessage += "Minimum Length must be an integer.\n"
        if len(errorMessage) > 0:
            errorMessage = errorMessage.rstrip("\n")
            return render_template("index.html", hasResults=False,
                hasErrors=True, errorMessage=errorMessage)
        else:
            writeUserData(name, affiliation, email, organism)
            enzymes = []
            chooser = EnzymeChooser(sequence, enzymesFile, 15, int(minLength))
            #return Response(chooser.getEnzymes(), mimetype= 'text/event-stream')
            flash("this could take a few minutes")
            chooser.getEnzymes(enzymes)
            #rint(enzymes)
            return render_template("index.html", hasResults=True,
                enzymes=enzymes)
    return render_template("index.html", hasResults=False, hasErrors=False)

@app.route("/feedback", methods=("GET", "POST"))
def feedback():
    # if they submitted the form
    if request.method == 'POST':
        feedback = request.form['feedback']
        errorMessage = ""
        # input validation
        if len(feedback.split()) > 200:
            errorMessage += "Please keep the feedback to 200 words.\n"
        if len(feedback) < 1:
            errorMessage += "Feedback too short.\n"
        if len(errorMessage) > 0:
            errorMessage = errorMessage.rstrip("\n")
            return render_template("feedback.html",
                hasErrors=True, errorMessage=errorMessage)
        else:
            writeFeedback(feedback)
            return render_template("feedback.html", submitted=True)
    return render_template("feedback.html", hasErrors=False)

@app.route("/instructions")
def instructions():
    return render_template("instructions.html")

if __name__ == "__main__":
    #app.run(debug=True)
    app.run("host=0.0.0.0")
