'''
Write a program to
a. calculate 2-mer frequency of DNA sequences of OPN1SW, OPN1MW and
OPN1LW genes.
b. calculate sequence similarity using 2-mer frequency and Manhattan
distance.
2. Make your k-mer counting program more generalizable so it would take k value
as input.
3. Make your similarity calculator more generalizable so it would take distance
metric as input.
4. Measure runtime of your program for every k value from 1 to 5.
5. Plot your results to visualize the relationship between k value and runtime (x-
axis: k value and y-axis: runtime). I recommend Matplotlib plotting library for
Python programmers.
'''
import practical1a
testSeqA = "ATGGCCTTGA"
testSeqB = "ATTGGGTAA"


files = [
    "OPN1LW_datasets/ncbi_dataset/data/gene.fna",
    "OPN1MW_datasets/ncbi_dataset/data/gene.fna",
    "OPN1SW_datasets/ncbi_dataset/data/gene.fna",
]

mers = {}

def kFrequencyAnalysis(sequence, k, seqNum):
    

    # find mers in each sequence
    for pos in range(len(sequence)-(k-1)):
        currentMerKey = sequence[pos:pos+k]
        currentMerDict = {seqNum: 1}

        # add to dict
        if currentMerKey in mers:
            if seqNum in mers[currentMerKey]:
                mers[currentMerKey][seqNum] += 1
            else:
                updateValue = mers.get(currentMerKey, {}) 
                updateValue.update(currentMerDict)
                mers[currentMerDict] = updateValue
            
            
        else:
            mers[currentMerKey] = currentMerDict

        # test
    print(mers)




# def manhattanDistance():
#     # Dm = sum of abs(a-b)

#     for pos in allMers:



def kMerCounting():
    # user input k value
    k = 2
    # read files with sequence
    seqStrings = []
    count = 0
    # for file in files:
    #     seqences = practical1a.readFile(file)
    #     count += 1

    kFrequencyAnalysis(testSeqA, k, 1)
    kFrequencyAnalysis(testSeqB, k, 2)

    

if __name__ == "__main__":
    kMerCounting()
