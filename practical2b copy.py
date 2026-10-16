'''
code to detect open reading frames in a DNA sequence
- orf = section of dna/rna that begins with start codon, continues as series of base ripletes and ends with stop codon

must consider 6 frame translation whilst searching for ORFs
    - a peice of dna has 6 possible reading frames
    - get 3 possible combinations by starting search for start codon from 1st, 2nd or 3rd nucleotide - +1, +2, +3
    - get other 3 going in reverse order -1, -2, -3

generate output file with:
    number of ORFs in sequence
    length (num nucleotides) 
    starting position (frame translation number and index of first nucleotide) of each ORF




'''

stop_codons = ["UAA", "UAG", "UGA"]

import practical1a, practical2a




def createORF(sequence):
    triplets = []
    for character in range(start, len(sequence), 3):
        triplets.append(sequence[character:character + 3])

        # print(f"codon list {codons}")

        # translate each codon using codon table

    openReadingFrame = ""
        
    for triplet in triplets:
        # rna fasta has T not U
        tripletReplaceT = triplet.replace("T", "U")
        
        translation = practical2a.codonTable.get(tripletReplaceT)

        if translation == "stop":
            openReadingFrame += triplet
            # ORF found
            lengthORF = len(openReadingFrame)
            print(f"Starting position of ORF: {start}\nLength of ORF: {lengthORF}")
            return openReadingFrame


            
            
        else:
            if translation is not None:
                openReadingFrame += triplet
            else:
                print(f"triplet {triplet} not found in codon table")


def findORF(sequence):
    amededSeq = "split sequence somehow for frame position thing"
    for frameNum in range (1, 6):
    
        # only start translation at start codon
        start = sequence.find("ATG")
        
        createORF(sequence, start)



   
        



def main():
    sequences = practical1a.readFile("OPN1LW_datasets/ncbi_dataset/data/gene.fna")
    for sequence in sequences:
        findORF(sequence)


if __name__ == "__main__":
    main()