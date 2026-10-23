def gcContent(sequence):
    total = len(sequence)

    gcCount = 0
    for base in sequence:
        if sequence[base] == "G" or sequence[base] == "C":
            gcCount+=1

    result = (gcCount / total) * 100
    return result

def main():
    sequence = input("Enter sequence")
    # normalise
    sequence = sequence.strip
    sequence = sequence.capitalize
    result = gcContent(sequence)
    print(result)
