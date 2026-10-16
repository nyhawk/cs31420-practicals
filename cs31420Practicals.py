import practical1a, practical1b, practical2a, practical2b

def main():
    print("================================")
    print("CS31420 practicals")
    print("================================")
    selection = input("Select a program to run\n"
    "1 - FASTA file reader - outputs length of sequences in a FASTA file\n"
    "2 - Inheritance probability calculator\n"
    "3 - mRNA translation to amino acids\n"
    "4 - ORF finder")

    match selection.strip():
        case "1":
            practical1a.fastaFileReader()
        case "2":
            practical1b.inheritance()
        case "3":
            practical2a.main()
        case "4":
            practical2b.main()


if __name__ == "__main__":
    main()