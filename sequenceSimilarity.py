'''
Write a program to calculate the similarity between two amino acid sequences.
For simplicity, you may assume that two sequences have the same number of
amino acids (but be aware that this is rarely the case). To test your program,
download one of the orthologs of Homo sapiens OPN1LW from the NCBI
database (choose one ortholog that has the same number of amino acids).

In an essence, Task 2 asks you to evaluate the similarity of two amino acid
sequences. In the most basic sense, similarity can be measured in a binary
grouping: different (0) versus identical (1). However, in bioinformatics we would be
better of quantifying it using a range of real values which indicate how much
different or similar two amino acid sequences are (e.g., 90%)
'''
import practical1a
import levenshteinDistance

# seq1 = "MATCPSGWLDK"
# seq2 = "MPWKFLDGSLAAFHJSL"


def validate(seq1, seq2):
    # check they are the same length else calculation incorrect
    if len(seq1) != len(seq2):
        print("amino acid lengths must be the same")
        return 0
    elif seq1 == "" or seq2 == "":
        "sequence cannot be empty"
        return 1
    

def similarity(seq1, seq2):
    validation = validate(seq1, seq2)
    if validation is None:
        count = 0
        for num in range (0, len(seq1)):
            if seq1[num] == seq2[num]:
                count += 1

        similartityScore = round(count/len(seq1), 5) 
        percentageSimilarity = similartityScore * 100
        print(f"Percentage similarity (using Hamming distance), rounded to 3dp: {percentageSimilarity}%")
        return percentageSimilarity

def main():
    seq1List = practical1a.readFile("OPN1LW_datasets/ncbi_dataset/data/protein.faa")
    seq2List = practical1a.readFile("OPN1LW_chimpanzee_datasets/ncbi_dataset/data/protein.faa")

    header1, seq1 = seq1List[0]
    header2, seq2 = seq2List[0]

    similarity(seq1, seq2)

   
    # timeLev = 3 ** len(seq2List[0])
    timeLev = len(seq1) * len(seq2)
    confirm = input(f"running levenshtein will take {timeLev}\nRun? y/n\n")
    if confirm == "y":
        lev = levenshteinDistance.levenshteinFullMatrix(seq1, seq2)
        print(f"Levenshtein distance: {lev}")

if __name__ == "__main__":
    main()