'''
Write a program to
DONE a. calculate 2-mer frequency of DNA sequences of OPN1SW, OPN1MW and
OPN1LW genes.
DONE b. calculate sequence similarity using 2-mer frequency and Manhattan
distance.
DONE 2. Make your k-mer counting program more generalizable so it would take k value
as input.
DONE 3. Make your similarity calculator more generalizable so it would take distance
metric as input.
4. Measure runtime of your program for every k value from 1 to 5.
5. Plot your results to visualize the relationship between k value and runtime (x-
axis: k value and y-axis: runtime). I recommend Matplotlib plotting library for
Python programmers.
'''
import practical1a, math, timeit, time
import matplotlib.pyplot as plt

testSeqA = "ATGGCCTTGA"
testSeqB = "ATTGGGTATA"

startTime = 0
endTime = 0


files = [
    "OPN1LW_datasets/ncbi_dataset/data/gene.fna",
    "OPN1MW_datasets/ncbi_dataset/data/gene.fna",
    "OPN1SW_datasets/ncbi_dataset/data/gene.fna",
]

mers = {}

def normaliseMers():
    totals = {}

    for countsBySequence in mers.values():
        for seqID, count in countsBySequence.items():
            totals[seqID] = totals.get(seqID, 0) + count

    frequencies = {}
    for mer, countsBySequence in mers.items():
        frequencies[mer] = {
            seqID: count / totals[seqID]
            for seqID, count in countsBySequence.items()
        }
    return frequencies                

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

def euclideanDistance(count1, count2):
    distance = math.sqrt((count1 - count2)**2)
    return distance


def manhattanDistance(count1, count2):
    # Dm = sum of abs(a-b)
    distance = abs(count1 - count2)
    return distance


def calculateSimilarity(seqID1, seqID2, frequencies, distanceMetric):       
    # Calculate the distance using the passed-in metric
    distance = 0
    for countsPerSeq in frequencies.values():
        count1 = countsPerSeq.get(seqID1, 0)
        count2 = countsPerSeq.get(seqID2, 0)
        distance += distanceMetric(count1, count2)
    
   
    similarity = 1 / (1 + distance)
    return similarity

def formatHeader(header):
    header = header.split(" ")
    return header[0]

def kMerCounting(k):
   
    # read files with sequence
    sequenceIDs = []
    records = []
    for file in files:
        records = practical1a.readFile(file)

        for header, sequence in records:
            shortHeader = formatHeader(header)
            kFrequencyAnalysis(sequence, k, shortHeader)
            if shortHeader not in sequenceIDs:
                sequenceIDs.append(shortHeader)

    for mer, counts in sorted(mers.items()):
        print(mer, counts) 

    return sequenceIDs

    
    
    endTime = timeit.default_timer()
    runtime = endTime - startTime
    print(f"Runtime for k={k} is {runtime}")
    
def distance():
    # user input k value
    k = int(input("enter K value\n").strip())
    
    sequenceIDs = kMerCounting(k) # do k frequency first
    
    # comparing every pair of sequences of OPN1 
    frequencies = normaliseMers()
    
    distanceMetricInp = int(input("Choose distance metric:\n"
                           "1. Manhattan distance\n"
                           "2. Euclidean distance\n"))
    distanceMetric = manhattanDistance
    match distanceMetricInp:
        case "1":
            distanceMetric = manhattanDistance
            
        case "2":
            distanceMetric = euclideanDistance
        
    for i in range(len(sequenceIDs)):
        for j in range(i + 1, len(sequenceIDs)):
            distance = calculateSimilarity(sequenceIDs[i], sequenceIDs[j], frequencies, distanceMetric)
            print(f"Distance between sequences {sequenceIDs[i]} and {sequenceIDs[j]}: {distance}\n")
 
    
def findingRuntime():
    runtimes = []
    startTime = time.time()
    for k in range(1, 6):
        kMerCounting(k)
        endTime = time.time()
        runtime = endTime - startTime
        runtimes.append(runtime)
        
    print(runtimes)
    
    # plotting values
    x = range(1, 6)
    y = runtimes
    
    plt.plot(x,y, "o:k")
    
    plt.title("Runtimes for k-mer analysis on OPN1 DNA sequences for LW, MW and SW variants for humans")
    plt.ylabel("Runtime (seconds)")
    plt.xlabel("k values")
    #plt.xticks(x, ['1', '2', '3', '4', '5'])
    plt.show()
    
    
def main():
    inputVal = int(input("Choose calculation:\n"
                "1. k-mer frequency analysis\n"
                "2. sequence similarity\n"
                "3. Runtime for k=1 to k=5\n"))
    match inputVal:
        case 1:
            # user input k value
            k = int(input("enter K value\n").strip())
            kMerCounting(k)
        case 2:
            distance()
        case 3:
            findingRuntime()
    
            
            
if __name__ == "__main__":
    main()
