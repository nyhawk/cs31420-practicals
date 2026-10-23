# -*- coding: utf-8 -*-
stopCodons = ["TAA", "TAG", "TGA"]
def reverseComplement(sequence):
    # generate reverse complement of given sequence

    # dictionary of complimenting bases
    complement = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G"
        }
    
    reverseComplementSeq = ""

    for nucleotide in reversed(sequence):
        reverseComplementSeq += complement[nucleotide]

    return reverseComplementSeq


def writeToFile(name, start, stop, length, frame):
    # format data
    data = f"ORF{name}, frame: {frame}, start position: {start}, stop position: {stop}, ORF length: {length}\n"


    with open ("orfData.txt", "a") as file:
        file.writelines(data)

def findOrfs(sequence: str):
    # generate other strand for 6-frame translation
    reverseSeq = reverseComplement(sequence)

# frames for 6-frame translation
    frames = [
        ("+1", sequence, 0),
        ("+2", sequence, 1),
        ("+3", sequence, 2),
        ("-1", reverseSeq, 0),
        ("-2", reverseSeq, 1),
        ("-3", reverseSeq, 2),
    ]
 
    startPos = 0
    stopPos = 0
    orfLength = 0
    orfNum = 0
    startFound = False

    #  for each frame
    for frameName, seq, offset in frames:
        # iterate through sequence
        for pos in range (offset, len(seq)-2, 3):
            codon = seq[pos:pos+3]
            # print (codon) DEBUG

            if (codon == "ATG"):
                # start codon found
                startPos = pos + 1
                startFound = True

            elif (codon in stopCodons):
                stopPos = pos + 3

            #    only consider as frame if a start codon has been found in the same frame
                if (startFound == True):
                     # process sequence
                    # in-frame check
                    orfLength = (stopPos) - (startPos) + 1
                    if ((orfLength % 3 == 0) and (orfLength > 75)):
                        # no remainder so divides by 3 and passes check

                        # write to file
                        orfNum += 1
                        
                        writeToFile(orfNum, startPos, stopPos, orfLength, frameName)
                        # print(f"start {startPos}, stop {stopPos}, stop codon {codon}") DEBUG print

                        # reset position variables
                        startPos = 0
                        stopPos = 0
                        startFound = False




def main():
  # first read file using practical1a code
      import practical1a
  
      sequences = practical1a.readFile("OPN1LW_datasets/ncbi_dataset/data/gene.fna")
      for sequence in sequences:
          findOrfs(sequence)

if __name__ == "__main__":
    main()