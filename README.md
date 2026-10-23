# CS31420 Practicals

A collection of Python practical exercises for the CS31420 course.

## Requirements

- Python 3

## Getting started

Clone the repository and navigate to its directory:

```bash
git clone <repository-url>
cd cs31420-practicals
```

To keep any installed packages separate from your system Python, you can create a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

On Windows, activate the environment with:

```powershell
.venv\Scripts\activate
```

## Running a practical

Find the Python file for the practical you want to run, then execute it with Python 3:

```bash
python3 path/to/practical.py
```

You can also run all practicals as one program using the main Python file.

```bash
python3 path/to/cs31420Practicals.py
```

Install any additional dependencies specified by that practical before running it. Check files in the practical's directory for exercise-specific instructions.

FASTA files read in the scripts are hardcoded currently. This may be improved in the future.

# About each practical

## Practical 1a
Reading FASTA files and calculating the length of nucleotide sequences

## Practical 1b
Calculates the probability of a male or female child inheriting a defective gene using X-linked or autosomal inheritance

## Practical 2a
Translates a sequence of RNA into amino acids

## Practical 2b - filename orfs.py
Find Open Reading Frames in a sequence of DNA

## Practical 3
1. Calculates GC content - gcContent.py
2. Calculates similarity between sequences of amino acids by comparing each letter at each position in a sequence

An algorithm for calculating Levenshtein distance is included but has not been written by me