#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import vcf
import pandas as pd
import json
import csv


list_of_chromosomes2 = {'EGFR': ['NM_001346897.2', 'NM_001346898.2', 'NM_001346899.2', 'NM_001346900.2', 'NM_001346941.2', 'NM_005228.5', 'NM_201282.2', 'NM_201283.2', 'NM_201284.2', 'XM_047419952.1', 'XM_047419953.1'], 'NRAS': ['NM_002524.5'], 'KRAS': ['NM_001369786.1', 'NM_001369787.1', 'NM_004985.5', 'NM_033360.4', 'XM_047428826.1'], 'PIK3CA': ['XM_006713658.5'], 'BRAF': ['NM_001354609.2', 'NM_001374244.1', 'NM_001374258.1', 'NM_001378467.1', 'NM_001378468.1', 'NM_001378469.1', 'NM_001378470.1', 'NM_001378471.1', 'NM_001378472.1', 'NM_001378473.1', 'NM_001378474.1', 'NM_001378475.1', 'NM_004333.6', 'XM_017012559.2', 'XM_047420766.1', 'XM_047420767.1', 'XM_047420768.1', 'XM_047420769.1', 'XM_047420770.1'], 'CTNNB1': ['NM_001098209.2', 'NM_001098210.2', 'NM_001330729.2', 'NM_001904.4', 'XM_006712985.2', 'XM_017005738.2', 'XM_024453356.2', 'XM_047447477.1', 'XM_047447478.1', 'XM_047447479.1', 'XM_047447480.1', 'XM_047447481.1', 'XM_047447482.1', 'XM_047447483.1'], 'MET': ['NM_000245.4', 'NM_001127500.3', 'NM_001324401.3', 'NM_001324402.2', 'XM_011516223.2', 'XM_047420400.1']}


list_of_patients=["LC_C1", "LC_C7", "LC_S14", "LC_S20"]

gene_names = ["EGFR", "NRAS", "KRAS", "PIK3CA", "BRAF", "CTNNB1", "MET"]


patient_mutational_profile = []
mutation_profile={}

# for gene_name in gene_names:
#     variant_line = "Homo sapiens " + gene_name + " proto-oncogene"
#     variant_line_2 = "(" + gene_name + "), transcript variant"
#     db_file = open("/home/thompson/data/transcriptome/reference/GCF_000001405.40_GRCh38.p14_rna.fna", 'r')
#     print(gene_name)

#     for line in db_file:

#         if variant_line in line:
#             print(line.split()[0][1:])
#             list_of_chromosomes2[gene_name].append(line.split()[0][1:]) 

#         elif variant_line_2 in line:
#             print(line.split()[0][1:])
#             list_of_chromosomes2[gene_name].append(line.split()[0][1:]) 

for patient in list_of_patients:

    vcf_reader = vcf.Reader(open('/home/thompson/project/healthy_tissue/' + patient + '/cancer_only_calls_80.vcf', 'r'))

    list_of_mutations=[]
    list_of_mutations2=[]

    print(patient)

    for record in vcf_reader:
        for gene in list_of_chromosomes2:
            if record.CHROM in list_of_chromosomes2[gene]:
                chrom=gene
                print(chrom, record.CHROM, record.POS, record.REF, record.ALT)
                list_of_mutations2.append(chrom) #+ "_" + str(record.POS))
                new_mutation={"Patient": patient, "Gene":chrom, "refseq_id":record.CHROM, "position":record.POS, "reference":record.REF, "alteration":record.ALT}
                list_of_mutations.append(new_mutation)
                break



    new_entry={"patient": patient, "list of mutations": (list_of_mutations)}
    # new_entry2={patient: list_of_mutations2}
    mutation_profile[patient]=list_of_mutations2
    patient_mutational_profile.append(new_entry) 



# print(patient_mutational_profile)
# for patient in mutation_profile:
#     if patient == 'LC_C1':
#         mutational_proflie_list = (mutation_profile[patient])
#         for item in mutational_proflie_list:

#             print(item)

with open('/home/thompson/project/under_60_data.csv', 'w', newline='') as file:
    writer = csv.writer(file)

    # Write header row (if needed)
    writer.writerow(['Patient', 'Mutation'])
    for patient in mutation_profile:
        if patient == 'LC_C1' or patient =='LC_S20':
            # writer.writerow([patient, ''])
            mutational_proflie_list = (mutation_profile[patient])
            for item in mutational_proflie_list:

                print(item)
                # Write data to specific columns
                writer.writerow([patient, item])  # Leave Column3 empty
                # writer.writerow(['', 'Value3', 'Value4'])  # Leave Column1 empty

# print(mutation_profile)

with open('/home/thompson/project/over_60_data.csv', 'w', newline='') as file:
    writer = csv.writer(file)

    # Write header row (if needed)
    writer.writerow(['Patient', 'Mutation'])
    for patient in mutation_profile:
        if patient == 'LC_C7' or patient =='LC_S14':
            # writer.writerow([patient, ''])
            mutational_proflie_list = (mutation_profile[patient])
            for item in mutational_proflie_list:

                print(item)
                # Write data to specific columns
                writer.writerow([patient, item])  # Leave Column3 empty
                # writer.writerow(['', 'Value3', 'Value4'])  # Leave Column1 empty

