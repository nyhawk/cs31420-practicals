"""Estimate inheritance and phenotype probabilities for a single gene.

Autosomal genotypes may be written as ``A/a`` (or ``Aa``). 
For sex-linked augenes, use chromosome-labelled alleles, for example ``X^D/X^d`` and
``X^D/Y``. 
Xx and xY notation also accepted.
The defective allele is ``d`` by default.
"""


from itertools import product
import re


def _normalise_defective_allele(defective_allele):
    value = (defective_allele or "").strip()
    if not value:
        return "d"
    if "^" in value:
        value = value.split("^", 1)[1]
    return value


def _autosomal_alleles(genotype):
    # Remove whitespace and split into alleles
    text = genotype.strip()
    if "/" in text:
        alleles = [part.strip() for part in text.split("/")]
    elif "," in text:
        alleles = [part.strip() for part in text.split(",")]
    elif len(text) == 2:
        alleles = list(text)
    else:
        raise ValueError("Use two alleles, such as A/a or Aa.")
    if len(alleles) != 2 or any(not allele for allele in alleles):
        raise ValueError("An autosomal genotype must contain exactly two alleles.")
    return alleles


def _normalise_sex_linked_shorthand(genotype):
    """Accept common shorthand like Xx, xX, xY, and XY and convert them to the canonical X^... format."""
    text = genotype.strip().replace(" ", "")
   
    # female alleles xx
    if re.fullmatch(r"[x][x]", text):
        return "X^d/X^d"

    # female alleles XX
    elif re.fullmatch(r"[X][X]", text):
        return "X^D/X^D"

    # female alleles Xx or xX
    elif re.fullmatch(r"[Xx][Xx]", text):
        if text[0].islower():
            return "X^d/X^D"
        return "X^D/X^d"
    
    # male alleles XY
    elif re.fullmatch(r"[X][Y]", text):
        return "X^D/Y^D"

    # male alleles xy
    elif re.fullmatch(r"[x][y]", text):
        return "X^d/Y^d"

    # male alleles Xy or xY
    elif re.fullmatch(r"[Xx][Yy]", text):
        if text[0].islower():
            return "X^d/Y^D"
        return "X^D/Y^d"

    return genotype


def _sex_chromosomes(genotype):
    # Allow shorthand input such as Xx or xY before parsing canonical chromosome forms.
    genotype = _normalise_sex_linked_shorthand(genotype)
    compact = genotype.strip().replace(" ", "")

    compact_match = re.fullmatch(
        r"([XxYy])(?:\^([A-Za-z0-9*+-]+))?([XxYy])(?:\^([A-Za-z0-9*+-]+))?",
        compact,
    )
    if compact_match:
        first_chrom, first_allele, second_chrom, second_allele = compact_match.groups()
        return [
            (first_chrom.upper(), first_allele),
            (second_chrom.upper(), second_allele),
        ]

    # Remove whitespace and split into alleles
    chromosomes = []
    for token in re.split(r"[/,\s]+", genotype.strip()):
        if not token:
            continue
        # Match X or Y chromosome with optional allele label, e.g. X^D or Y^d
        match = re.fullmatch(r"([XxYy])(?:\^([A-Za-z0-9*+-]+))?", token)
        if not match:
            # Accept shorthand like Xx or xY after normalisation, but reject unsupported forms.
            raise ValueError("Use chromosome-labelled alleles, e.g. X^d/X^D, X^D/Y, Xx, or xY.")
        chromosome, allele = match.groups()
        chromosomes.append((chromosome.upper(), allele))

    return chromosomes


def calculate_probability(mother_genotype, father_genotype, offspring_sex,
                          inheritance, chromosome, defective_allele="d"):
    """Return (probability of inheriting the allele, probability of phenotype)."""
    # Normalise all user input so comparisons are consistent.
    sex = offspring_sex.strip().lower()
    inheritance = inheritance.strip().lower()
    chromosome = chromosome.strip().lower()
    defective_allele = _normalise_defective_allele(defective_allele)

    # Validate the basic inputs first.
    if sex not in {"male", "female", "m", "f"}:
        raise ValueError("Offspring sex must be male or female.")
    is_male = sex in {"male", "m"}

    if inheritance not in {"dominant", "recessive"}:
        raise ValueError("Inheritance must be dominant or recessive.")

    # ------------------------------------------------------------------
    # X-linked inheritance
    # ------------------------------------------------------------------
    # A male gets only one X chromosome, so he inherits it only from his mother.
    # A female gets one X from each parent, so she may inherit 0, 1 or 2 defective X alleles.
    if chromosome in {"x", "x-linked", "x linked"}:
        mother = _sex_chromosomes(mother_genotype)
        father = _sex_chromosomes(father_genotype)

        # check that the mother has two X chromosomes and the father has one X and one Y chromosome
        maternal_x = [allele for chrom, allele in mother if chrom == "X"]
        paternal_x = [allele for chrom, allele in father if chrom == "X"]
        paternal_y = [allele for chrom, allele in father if chrom == "Y"]

        if len(maternal_x) != 2 or len(paternal_x) != 1 or len(paternal_y) != 1:
            raise ValueError("For X-linked inheritance, use a mother with two Xs and a father with one X and one Y (e.g. X^D/Y).")

        # Build every genotype combination for the child.
        if is_male:
            # Male child gets one maternal X and the father's Y. The paternal X is irrelevant.
            # children = [(maternal, None) for maternal in maternal_x]
            children = list(product(maternal_x, paternal_y))
        else:
            # Female child gets one maternal X and one paternal X.

            children = list(product(maternal_x, paternal_x))
            # children = [(maternal, paternal) for maternal, paternal in product(maternal_x, paternal_x)]

        # Count how many defective alleles appear in each possible outcome.
        # defective_counts keeps track of how many children have 0, 1, or 2 defective alleles.
        defective_counts = {0: 0, 1: 0, 2: 0}  # Initialize counts for 0, 1, and 2 defective alleles.
        for child in children:
            # For each child outcome, count how many X alleles match the defective allele.
            childs_defective_alleles = 0
            for allele in child:
                if allele is not None and allele == defective_allele:
                    childs_defective_alleles += 1
            # Increment the count for the number of defective alleles found in this child.
            defective_counts[childs_defective_alleles] += 1

        total = len(children)

        # calculate the probabilities

        inherited = 0.0
        inherited_count = 0
        # Probability of inheriting at least one defective allele.
   
        # for num_defective, count in sorted(defective_counts.items(), reverse=True):
        for num_defective, count in defective_counts.items():
            if num_defective > 0 and count > 0:
                inherited_count += count
      
                        
        inherited = inherited_count / total

        # Probability of expressing the phenotype.
        affected = 0.0
        for num_defective, count in defective_counts.items():
            if is_male:
                # One defective X chromosome is enough for a male to be affected.
                if num_defective > 0:
                    affected += count
            else:
                # Female phenotype depends on whether the trait is dominant or recessive.
                if inheritance == "dominant" and num_defective > 0:
                    affected += count
                elif inheritance == "recessive" and num_defective == 2:
                    affected += count

        affected /= total
        return inherited, affected

    # ------------------------------------------------------------------
    # Y-linked inheritance
    # ------------------------------------------------------------------
    # Y-linked traits are passed only from father to son.
    if chromosome in {"y", "y-linked", "y linked"}:
        # probability of inheriting defective allele is 0 for females
        if not is_male:
            return 0.0, 0.0

        father = _sex_chromosomes(father_genotype)
        paternal_y = [allele for chrom, allele in father if chrom == "Y"]

        if len(paternal_y) != 1:
            raise ValueError("The father's genotype must include one Y chromosome, e.g. X/Y^d.")

        probability = float(paternal_y[0] == defective_allele)
        return probability, probability

    # ------------------------------------------------------------------
    # Autosomal inheritance
    # ------------------------------------------------------------------
    # Each parent contributes one allele, so every child outcome has 2 alleles total.
    mother = _autosomal_alleles(mother_genotype)
    father = _autosomal_alleles(father_genotype)

    # create punnet square
    children = list(product(mother, father))

    defective_counts = {0: 0, 1: 0, 2: 0}
    # {
    # 0 defective alleles in pair: value,
    # 1 defective allele in pair: value,
    # 2 defective alleles in pair: value
    # }
    for pair in children:
        childs_defective_alleles = 0
        for allele in pair:
            if allele == defective_allele:
                childs_defective_alleles += 1
        # Increment the count for the number of defective alleles found in this child.
        defective_counts[childs_defective_alleles] += 1

    total = len(children)

    # Probability of inheriting at least one defective allele.
    inherited = 0.0
    for num_defective, count in defective_counts.items():
        if num_defective > 0:
            inherited += count
    inherited /= total

    # Probability of expressing the phenotype.
    affected = 0.0
    for num_defective, count in defective_counts.items():
        if inheritance == "dominant":
            if num_defective > 0:
                affected += count
        else:  # recessive
            if num_defective == 2:
                affected += count

    affected /= total
    return inherited, affected


def main():
    import sys

    
    mother = input("Mother's genotype: ")
    father = input("Father's genotype: ")
    sex = input("Offspring sex (m/f): ")
    inheritance = input("Inheritance (dominant/recessive): ")
    chromosome = input("Chromosome (autosomal, X, or Y): ")
    defective = input("Defective allele, 'd' or 'D': ").strip() or "d"
    try:
        allele_probability, phenotype_probability = calculate_probability(
            mother, father, sex, inheritance, chromosome, defective
        )
    except ValueError as error:
        print(f"Error: {error}")
        return
    
    print(f"Probability of {sex} inheriting at least one defective allele: {allele_probability:.1%}")
    print(f"Probability of {sex} expressing the phenotype: {phenotype_probability:.1%}")


if __name__ == "__main__":
    main()
