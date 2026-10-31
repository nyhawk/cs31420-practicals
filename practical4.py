'''
Write a program to
a. calculate 2-mer frequency of DNA sequences of OPN1SW, OPN1MW and
OPN1LW genes.
b. calculate sequence similarity using 2-mer frequency and Manhattan
distance.
DONE 2. Make your k-mer counting program more generalizable so it would take k value
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
testSeqB = "ATTGGGTATA"


files = [
    "OPN1LW_datasets/ncbi_dataset/data/gene.fna",
    "OPN1MW_datasets/ncbi_dataset/data/gene.fna",
    "OPN1SW_datasets/ncbi_dataset/data/gene.fna",
]

mers = {}

def kFrequencyAnalysis(sequence, k, seqID):
    

    # find mers in each sequence
    for pos in range(len(sequence)-(k-1)):
        currentMerKey = sequence[pos:pos+k]
        currentMerDict = {seqID: 1}

        # add to dict
        if currentMerKey in mers:
            if seqID in mers[currentMerKey]:
                mers[currentMerKey][seqID] += 1
            else:
                updateValue = mers.get(currentMerKey, {}) 
                updateValue.update(currentMerDict)
                mers[currentMerKey] = updateValue
            
            
        else:
            mers[currentMerKey] = currentMerDict

        # test
    print(mers)




def manhattanDistance():
    # Dm = sum of abs(a-b)
    

    for mer in mers:
        # get a mer eg "AT"
        # mer give just the key to the main dictionary
        # then get the values in each sequence and minus them
        print("\n", mers[mer][])
        # abs(me)



def kMerCounting():
    # user input k value
    k = 2
    # read files with sequence
    records = []
    count = 0
    for file in files:
        records = practical1a.readFile(file)

        for header, sequence in records:
            kFrequencyAnalysis(sequence, k, header)

    print("done") 
    manhattanDistance()

    # kFrequencyAnalysis(testSeqA, k, 1)
    # kFrequencyAnalysis(testSeqB, k, 2)

    

if __name__ == "__main__":
    kMerCounting()
