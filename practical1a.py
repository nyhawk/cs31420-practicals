# -*- coding: utf-8 -*-
import re

NUCLEIC_ACID_CODES = {"A", "C", "G", "T", "N", "U", "R", "Y", "K", 
                      "M", "S", "W", "B", "D", "H", "V", "-", "*"}
AMINO_ACID_CODES = {"A", "B", "C", "D", "E", "F", "G", "H", "I", 
                    "J", "K", "L", "M", "N", "O", "P", "Q", "R",
                    "S", "T", "U", "V", "W", "Y", "Z", "X", "*", "-"}

files = [
    "OPN1LW_datasets/ncbi_dataset/data/gene.fna",
    "OPN1LW_datasets/ncbi_dataset/data/protein.faa",
    "OPN1LW_datasets/ncbi_dataset/data/rna.fna",

    "OPN1MW_datasets/ncbi_dataset/data/gene.fna",
    "OPN1MW_datasets/ncbi_dataset/data/protein.faa",
    "OPN1MW_datasets/ncbi_dataset/data/rna.fna",

    "OPN1SW_datasets/ncbi_dataset/data/gene.fna",
    "OPN1SW_datasets/ncbi_dataset/data/protein.faa",
    "OPN1SW_datasets/ncbi_dataset/data/rna.fna",

# test datasets
    "test_datasets/blank_lines.fasta",
    "test_datasets/empty_header.fasta",
    "test_datasets/empty_sequence.fasta",
    "test_datasets/invalid_sequence_symbol.fasta",
    "test_datasets/mixed_case.fasta", 
    "test_datasets/multiple_sequences.fasta",
    "test_datasets/no_header.fasta",
    "test_datasets/sequence_before_header.fasta",
    "test_datasets/special_characters.fasta",
    "test_datasets/unexpected_header_marker.fasta"
]

def getSequenceType(headerLine):
    '''
    NC_ - genes - nucleic acid
    NM_ - RNA - nucleic acid
    NP_ - amino acid

    are many more
    '''
    sequence_type = headerLine.split("_")
    match sequence_type:
        case "NC":
            return "gene"
        case "NM":
            return "rna"
        case "NP":
            return "amino_acid"

    return "unknown"


def processSequence(data, sequenceType):
    # print("Length of sequence: ", len(data))
    
    # detect unique characters    
    # create dictionary
    sequenceDict = {}
    lengthOfSequence = 0

    # lower and upper case letters should be treated as the same base
    sequence = data.upper()
    
    # process each letter in string
    for character in sequence:        
        if character in sequenceDict:
            # increment counts
            sequenceDict[character] += 1
            lengthOfSequence += 1

        # character not in dictionary, validate and add
        else:
            # validate character using appropriate set of expected codes
            if sequenceType == "gene" or sequenceType == "rna":
                if character not in NUCLEIC_ACID_CODES:
                    print(f"Invalid {sequenceType} character detected: {character}")
                    continue

            elif sequenceType == "amino_acid":
                if character not in AMINO_ACID_CODES:
                    print(f"Invalid amino acid character detected: {character}")
                    continue

            elif sequenceType == "unknown":
                # if unknown what sequence type in file, check for alphabetical chars only
                match = re.fullmatch(r"[A-Za-z\*\-]", character)
                # if match is not None: # if match found
                # if character.isascii() and character.isalpha(): # without regex 
                    
                if match is None: # no match
                    print(f"Invalid character detected: \"{character}\"")
                    continue

            # passes validation, add to dictionary
            sequenceDict[character] = 1
            lengthOfSequence += 1

    print(f"Length of sequence: {lengthOfSequence}")  
    print(f"{sequenceDict}\n")
    return sequenceDict

    
def readFile(filename) -> list[str]:
    with open(filename) as sequence_file:   
        current_header = None
        sequence = ""
        sequenceType = "unknown"
        sequences = []
        print(filename, "\n------------------------\n")    

        for line in sequence_file:          
            line = line.strip()

            if line.startswith(">"): 
                # get type of sequence data
                sequenceType = getSequenceType(line)                

                # typical file data             
                if current_header is not None and sequence:
                    print("Processing: ", current_header)
                    sequences.append(sequence)
                    processSequence(sequence, sequenceType)                
                
                # if header has no following sequence
                elif current_header is not None and sequence == "":
                    print(f"{line} has no sequence")

                # if sequence before header
                elif current_header is None and sequence:
                    print(f"Sequence found before header line.")  
                    userInput = input("Would you like to analyse this sequence? y/n ")  
                    if userInput == "y":
                        sequences.append(sequence)
                        processSequence(sequence, sequenceType)

                # don't add header line to data list
                current_header = line
                sequence = ""   


            else :
                # sequence data only, add to list
                sequence += line

    # at end of file
    if current_header is not None and sequence:
        print("Processing: ", current_header)
        sequences.append(sequence)
        processSequence(sequence, sequenceType)    
        sequence = ""    

    elif current_header is None and sequence:
        print(f"No header found in {filename}")
        processSequenceNoHeader = input("Process sequence anyway? y/n ")
        
        if processSequenceNoHeader == "y":
            print(f"\nProcessing sequence in {filename} \nSequence data treated as one sequence")
            sequences.append(sequence)
            processSequence(sequence, sequenceType)    
            sequence = ""   

    return sequences


def fastaFileReader():
    print("\n========================\nFasta file reader\n========================\n")

    for file in files: 
        readFile(file)     
