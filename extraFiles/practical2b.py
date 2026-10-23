'''
code to detect open reading frames in a DNA sequence
- orf = section of dna/rna that begins with start codon, continues as series of base ripletes and ends with stop codon

must consider 6 frame translation whilst searching for ORFs
    - a peice of dna has 6 possible reading frames
    - get 3 possible combinations by starting search for start codon from 1st, 2nd or 3rd nucleotide - +1, +2, +3
    - get other 3 on opposite strand -1, -2, -3

generate output file with:
    number of ORFs in sequence
    length (num nucleotides) 
    starting position (frame translation number and index of first nucleotide) of each ORF




'''

stop_codons = ["TAA", "TAG", "TGA"]

import practical1a, practical2a

def reverse_complement(sequence):
    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    rc = ""

    for nucleotide in reversed(sequence):
        rc += complement[nucleotide]

    return rc



def findORFs(sequence):
    orfs = []
    reverse_seq = reverse_complement(sequence)
    original_length = len(sequence)

    sequences = [
        ("+1", sequence, 0),
        ("+2", sequence, 1),
        ("+3", sequence, 2),
        ("-1", reverse_seq, 0),
        ("-2", reverse_seq, 1),
        ("-3", reverse_seq, 2),
    ]

    # for each possible frame
    for frame_name, seq, offset in sequences:

        # find start codon
        # step 3 as searching for codons in base triplets
        for startCodonPos in range(offset, len(seq) - 2, 3):    
            codon = seq[startCodonPos:startCodonPos+3]    
            if codon == "ATG":

                #  search for stop codon
                for stopCodonPos in range(startCodonPos + 3, len(seq) - 2, 3):    
                    stop = seq[stopCodonPos:stopCodonPos+3]
    
                    if stop in stop_codons:    
                        orf_length = stopCodonPos + 3 - startCodonPos
    
                        # Calculate coordinates in original sequence
                        if frame_name.startswith("+"):
                            start_pos = startCodonPos + 1 # as python index starts at 0
                            end_pos = stopCodonPos + 3 # stopCodonPos = start of stop codon so +2 to get to end, +1 for index start at 0
    
                        else:
                            start_pos = original_length - (startCodonPos)
                            end_pos = original_length - stopCodonPos - 2
    
                        orfs.append({
                            "frame": frame_name,
                            "start": start_pos,
                            "end": end_pos,
                            "length": orf_length
                        })
    
                        break
    
    return orfs


   
        



def main():
    # sequences = practical1a.readFile("OPN1LW_datasets/ncbi_dataset/data/gene.fna")
    # sequence = "CCCATGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAGGG"
    # sequence = "CATGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAACCGGGGGGTTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCATG"
    sequence = "GCCAGGGCCGGCTGCCGTCGGGGACAGGGCTTTCCATAGCCATGGCCCAGCAGTGGAGCCTCCAAAGGCTCGCAGGCCGCCATCCGCAGGACAGCTATGAGGACAGCACCCAGTCCAGCATCTTCACCTACACCAACAGCAACTCCACCAGAGGCCCCTTCGAAGGCCCGAATTACCACATCGCTCCCAGATGGGTGTACCACCTCACCAGTGTCTGGATGATCTTTGTGGTCACTGCATCCGTCTTCACAAATGGGCTTGTGCTGGCGGCCACCATGAAGTTCAAGAAGCTGCGCCACCCGCTGAACTGGATCCTGGTGAACCTGGCGGTCGCTGACCTAGCAGAGACCGTCATCGCCAGCACTATCAGCATTGTGAACCAGGTCTCTGGCTACTTCGTGCTGGGCCACCCTATGTGTGTCCTGGAGGGCTACACCGTCTCCCTGTGTGGGATCACAGGTCTCTGGTCTCTGGCCATCATTTCCTGGGAGAGATGGATGGTGGTCTGCAAGCCCTTTGGCAATGTGAGATTTGATGCCAAGCTGGCCATCGTGGGCATTGCCTTCTCCTGGATCTGGGCTGCTGTGTGGACAGCCCCGCCCATCTTTGGTTGGAGCAGGTACTGGCCCCACGGCCTGAAGACTTCATGCGGCCCAGACGTGTTCAGCGGCAGCTCGTACCCCGGGGTGCAGTCTTACATGATTGTCCTCATGGTCACCTGCTGCATCATCCCACTCGCTATCATCATGCTCTGCTACCTCCAAGTGTGGCTGGCCATCCGAGCGGTGGCAAAGCAGCAGAAAGAGTCTGAATCCACCCAGAAGGCAGAGAAGGAAGTGACGCGCATGGTGGTGGTGATGATCTTTGCGTACTGCGTCTGCTGGGGACCCTACACCTTCTTCGCATGCTTTGCTGCTGCCAACCCTGGTTACGCCTTCCACCCTTTGATGGCTGCCCTGCCGGCCTACTTTGCCAAAAGTGCCACTATCTACAACCCCGTTATCTATGTCTTTATGAACCGGCAGTTTCGAAACTGCATCTTGCAGCTTTTCGGGAAGAAGGTTGACGATGGCTCTGAACTCTCCAGCGCCTCCAAAACGGAGGTCTCATCTGTGTCCTCGGTATCGCCTGCATGAGGTCTGCCTCCTACCCATCCCGCCCACCGGGGCTTTGGCCACCTCTCCTTTCCCCCTCCTTCTCCATCCCTGTAAAATAAATGTAATTTATCTTTGCCAAAACCAA"
    # for sequence in sequences:
    orfs = findORFs(sequence)

    print(f"Total ORFs found: {len(orfs)}\n")

    for number, orf in enumerate(orfs, start=1):
        print(f"ORF {number}")
        print(f"  Frame: {orf['frame']}")
        print(f"  Start nucleotide: {orf['start']}")
        print(f"  End nucleotide: {orf['end']}")
        print(f"  Length: {orf['length']} bp")
        print()




if __name__ == "__main__":
    main()