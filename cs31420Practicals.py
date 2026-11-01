import practical1a, practical1b, practical2a, gcContent, orfs, sequenceSimilarity, practical4

def main():
    print("================================")
    print("CS31420 practicals")
    print("================================")
    selection = input("Select a program to run\n"
    "1 - FASTA file reader - outputs length of sequences in a FASTA file\n"
    "2 - Inheritance probability calculator\n"
    "3 - mRNA translation to amino acids\n"
    "4 - ORF finder\n"
    "5 - GC content calculator\n"
    "6 - Similarity calculator\n"
    "7 - K-mer analysis and similarity calculator\n")

    match selection.strip():
        case "1":
            practical1a.fastaFileReader()
        case "2":
            practical1b.inheritance()
        case "3":
            practical2a.main()
        case "4":
            orfs.main()
        case "5":
            gcContent.main()
        case "6":
            sequenceSimilarity.main()
        case "7":
            practical4.main()



if __name__ == "__main__":
    main()