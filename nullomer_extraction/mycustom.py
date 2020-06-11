class GenericFile(object):
    def __init__(self, fileNameS):
            try:
                fileF = open(fileNameS)
                self.fileS = fileF.read()
                fileF.close()
            except:
                self.fileS = fileNameS.read()

import pdb

class FastaFile(GenericFile):
    """
    This class represents a FastaFile.
    """

    def __init__(self, fasta, fileName=True):
        """

        :rtype : object
        :param fasta: This can be either a file name (string) or a HTML file upload element.
        :param fileName: This is a boolean value. If False, then it is the fasta string
        :return: does not return anything.
        """
        if fileName == True:
            super(FastaFile, self).__init__(fasta)
        elif fileName == False:
            if ">" not in fasta:
                raise ValueError("It is not in Fasta Format")
            self.fileS = fasta
        else:
            raise Exception

        self.sequences = None
        self.sequences = self.getSequences()


    def __getitem__(self, item):
        if self.sequences is None:
            self.getSequences()
        return self.sequences[item]

    def __len__(self):
        return len(self.getSequences())

    def __repr__(self):
        return "< Fasta file: %s sequences >" % len(self)

    def __str__(self):
        return "< Fasta file: %s sequences >" % len(self)

    def getSequences(self):
        if self.sequences is None:
            firstSplitL = self.fileS.strip().split(">")[1:]
            eachByLineL = [i.splitlines() for i in firstSplitL]

            sequences = []
            for i in eachByLineL:
                header = i[0]
                sequence = "".join(i[1:]).upper()
                sequence = sequence.replace(" ", "")
                sequence = sequence.replace("\t", "")
                if sequence == "":
                    raise ValueError("There is no sequence !")
                sequences.append(Sequence(name=header, sequence=sequence))

            self.sequences = sequences

        return self.sequences

    def getMaxLength(self):
        """
        This returns the max length of a sequence in this object
        """
        sequences = self.getSequences()
        maxLength = 0
        for sequence in sequences:
            seqLength = len(sequence)
            if seqLength > maxLength:
                maxLength = seqLength
        return maxLength

    def areDuplicatesPresent(self):
        sequences = self.getSequences()
        sequenceL = [i.sequence for i in sequences]
        if len(sequenceL) > len(set(sequenceL)):
            return True
        return False

    def getNames(self):

        names = [sequence.name for sequence in self.sequences]
        return names


class Sequence(object):
    """
    This class represents a DNA/RNA sequence.
    Each object needs to have a name.
    """
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def __repr__(self):
        return "%s : %s" %(self.name, self.sequence)

    def __getitem__(self, item):
        return self.sequence[item]

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return "%s : %s " % (self.name, self.sequence)

