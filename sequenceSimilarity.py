# def distance(seq1, seq2, pos1, pos2):
#     if seq1[pos1-1] == seq2[pos2-1]:
#         return distance(seq1, seq2, pos1+1, pos2+1)
#     elif seq1[pos1] == seq1[pos2-1]:
#         distance(seq1, seq2, pos)

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

# aa1 = "MATCPSGWLDK"
# aa2 = "MPWKFLDGSLAAFHJSL"


def validate(aa1, aa2):
    # check they are the same length else calculation incorrect
    if len(aa1) != len(aa2):
        print("amino acid lengths must be the same")
        return 0

def similarity(aa1, aa2):
    validation = validate(aa1, aa2)
    if validation is None:
        count = 0
        for num in range (0, len(aa1)):
            if aa1[num] == aa2[num]:
                count += 1

        similartityScore = round(count/len(aa1), 5) 
        percentageSimilarity = similartityScore * 100
        print(f"Percentage similarity (Hamming distance), rounded to 3dp: {percentageSimilarity}%")
        return percentageSimilarity

def main():
    aa1List = practical1a.readFile("OPN1LW_datasets/ncbi_dataset/data/protein.faa")
    aa2List = practical1a.readFile("OPN1LW_chimpanzee_datasets/ncbi_dataset/data/protein.faa")

    similarity(aa1List[0], aa2List[0])

    # DO NOT DO THIS IT TAKES TOO LONG
    # timeLev = 3 ** len(aa2List[0])
    timeLev = len(aa2List[0]) * len(aa1List[0])
    confirm = input(f"running levenshtein will take {timeLev}\nRun? y/n\n")
    if confirm == "y":
        lev = levenshteinDistance.levenshteinFullMatrix(aa1List[0], aa2List[0]) #, len(aa1List[0]), len(aa2List[0]))
        print(f"Levenshtein distance: {lev}")

if __name__ == "__main__":
    main()