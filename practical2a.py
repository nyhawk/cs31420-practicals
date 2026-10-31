# codon table to link each codon to an amino acid
codonTable = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",

    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",

    "UAU": "Y",
    "UAC": "Y",
    "UAA": "stop",
    "UAG": "stop",

    "UGU": "C",
    "UGC": "C",
    "UGA": "stop",
    "UGG": "W",

    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",

    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",

    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",

    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",

    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M", #start

    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",

    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",

    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",

    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",

    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",

    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G"
}


def translateSequence(sequence) -> str:
    # split sequence into codons
    codons = []

    # only start translation at start codon
    start = sequence.find("ATG")
    
    for character in range(start, len(sequence), 3):
        codons.append(sequence[character:character + 3])

    # print(f"codon list {codons}")

    # translate each codon using codon table
    aminoAcids = ""
    
    for codon in codons:
        # rna fasta has T not U
        codonReplaceT = codon.replace("T", "U")
        
        translation = codonTable.get(codonReplaceT)
        if translation == "stop":
            print(f"amino acids {aminoAcids}")
            return aminoAcids         

        else:
            if translation is not None:
                aminoAcids += translation
            else:
                print(f"codon {codon} not found in codon table")


def checkTranslation(translation: str, actual: list[str]):
    
    # design and implement an algorithm (scoring metric) 
    # to compare two amino acid sequences to evaluate similarity

    totalMatch = 0

    for translationLetter, actualLetter in zip(translation, actual[0]):
        if translationLetter == actualLetter:
            totalMatch += 1

    percentageLikeness = (totalMatch / len(actual[0])) * 100
    print(f"Percentace likeness of translation to actual protein sequence: {percentageLikeness}%")
    

        # check letter at same position matches



def main():
    # print(codonTable)
    # translate mRNA sequence to protein sequence using codon table

    # first read file using practical1a code
    import practical1a

    sequences = practical1a.readFile("OPN1LW_datasets/ncbi_dataset/data/rna.fna")
    for header, sequence in sequences:
        translation = translateSequence(sequence)

    actualAminoAcids = practical1a.readFile("OPN1LW_datasets/ncbi_dataset/data/protein.faa")    
    
    checkTranslation(translation, actualAminoAcids[0][1])

if __name__ == "__main__":
    main()
