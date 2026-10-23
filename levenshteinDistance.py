def levenshteinRecursive(str1, str2, m, n):
    print(str1[m-1])
    print(str2[n-1])
    print("----")
      # str1 is empty
    if m == 0:
        return n
    # str2 is empty
    if n == 0:
        return m

    if str1[m - 1] == str2[n - 1]:
        return levenshteinRecursive(str1, str2, m - 1, n - 1)
    return 1 + min(
          # Insert     
        levenshteinRecursive(str1, str2, m, n - 1),
        min(
              # Remove
            levenshteinRecursive(str1, str2, m - 1, n),
          # Replace
            levenshteinRecursive(str1, str2, m - 1, n - 1))
    )

# Drivers code
str1 = "AAATG"
str2 = "AAAATG"
distance = levenshteinRecursive(str1, str2, len(str1), len(str2))
print("Levenshtein Distance:", distance)

# def distance(seq1, seq2, pos1, pos2):
#     if seq1[pos1-1] == seq2[pos2-1]:
#         return distance(seq1, seq2, pos1+1, pos2+1)
#     elif seq1[pos1] == seq1[pos2-1]:
#         distance(seq1, seq2, pos)