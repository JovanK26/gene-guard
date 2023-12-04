import random
import csv
from typing import List

num_strings = 500
string_length = 1000
common_patterns = ["ATCG", "AGCT", "TCGA", "GATC",'CGC', "ATA", "GCG", "AGA", "TGT", "CGA", "ATC", "GCG", "ACG", "GTG", "TAT", "AAT", "GAG", "AAA", "CGC"]
mutations = ["TAGCTGATCG","ACGACGCTC","ACGTCGATGT","AGCTGATCGATCGATCG"]
mutation_names=["Cystic Fibrosis","Hemophilia A","Hemophilia B","Spina Bifida"]
mutation_type=[1,0,0,0]
mutperc=[0,0]

def generate_dna_string(length, common_patterns):
    return ''.join(random.choice(common_patterns) for _ in range(length))

def generate_dataset(num_strings, string_length, common_patterns):
    dataset = [generate_dna_string(string_length, common_patterns) for _ in range(num_strings)]
    return dataset

def knuth_morris_pratt(pattern, text):
    m, n = len(pattern), len(text)
    lps = compute_lps(pattern)
    count = 0

    i, j = 0, 0
    while i < n:
        if pattern[j] == text[i]:
            i += 1
            j += 1

            if j == m:
                count += 1
                j = lps[j - 1]
        else:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    if count>1:
        count=1
    return count

def compute_lps(pattern):
    m = len(pattern)
    lps = [0] * m
    length = 0
    i = 1

    while i < m:
        if pattern[i] == pattern[length]:
            length += 1
            lps[i] = length
            i += 1
        else:
            if length != 0:
                length = lps[length - 1]
            else:
                lps[i] = 0
                i += 1

    return lps

def main():
    
    percentage_mutations=[]
    dataset = generate_dataset(num_strings, string_length, common_patterns)
    for mutation in mutations:
        cnt=0
        for dna in dataset:
            if knuth_morris_pratt(mutation,dna):
                cnt+=1
        percentage_mutations.append(100*(cnt/num_strings))

    
    for i, mutation in enumerate(mutations):
        print("Percentage of people who have "+ mutation_names[i]+": "+ str(percentage_mutations[i])+"%")
    for i,per in enumerate(percentage_mutations):
        mutperc[mutation_type[i]]+=per
    print("Percentage of healthy people: "+ str(100-mutperc[0]-mutperc[1])+"%")
    print("Percentage of category 2 people(medium risk): "+ str(mutperc[0])+"%")
    print("Percentage of category 1 people(highest risk): "+ str(mutperc[1])+"%")
    csv_filename = "generated_dataset.csv"
    with open(csv_filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["DNA String"])
        writer.writerows([[dna] for dna in dataset])

if __name__ == "__main__":
    main()
