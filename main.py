from collections import Counter, defaultdict
from itertools import combinations
from pyteomics import parser, mass
from tkinter import filedialog
from typing import List, Tuple
from datetime import datetime
from tkinter import *
from Bio import SeqIO
import tkinter as tk
import pandas as pd
import numpy as np
import shutil
import copy
import csv

"""dictionary/list of PTMs"""
translation_table = str.maketrans({'X': 'G', 'Z': 'Q', 'B': 'N'})
db = mass.Unimod()
aa_comp = dict(mass.std_aa_comp)

"""
Incomplete Dictionary of PTMs
aa_comp = dict(mass.std_aa_comp)
aa_comp['p'] = mass.Composition('HPO3')
aa_comp['g'] = mass.Composition('C6H12O6')  # does not account for replaced amide group on the amino acid
aa_comp['ox'] = mass.Composition('O')
aa_comp['ac'] = mass.Composition('C2H2O')
"""

"""caches for coverage maps"""
one_protease_cache_peptide = {}  # in tandem, peptide coverage
multi_protease_cache_peptide = {}  # scraped in series code, peptide coverage
protease_digestion_cache_protein = {}  # in tandem, protein identification rate

"""list of proteases supported by pyteomics, if more is added, the respective cleavage rule also needs to be added"""
proteaseList = ['arg-c', 'asp-n', 'bnps-skatole', 'caspase 1', 'caspase 2', 'caspase 3', 'caspase 4', 'caspase 5',
                'caspase 6', 'caspase 7', 'caspase 8', 'caspase 9', 'caspase 10', 'chymotrypsin high specificity',
                'chymotrypsin low specificity', 'clostripain', 'cnbr', 'enterokinase', 'factor xa', 'formic acid',
                'glutamyl endopeptidase', 'granzyme b', 'hydroxylamine', 'iodosobenzoic acid', 'lysc', 'ntcb',
                'pepsin ph1.3', 'pepsin ph2.0', 'proline endopeptidase', 'proteinase k', 'staphylococcal peptidase i',
                'thermolysin', 'thrombin', 'trypsin']

custom_expasy_rules = {
    'arg-c': r'R',
    'asp-n': r'\w(?=D)',
    'bnps-skatole': r'W',
    'caspase 1': r'(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',
    'caspase 2': r'(?<=DVA)D(?=[^PEDQKR])',
    'caspase 3': r'(?<=DMQ)D(?=[^PEDQKR])',
    'caspase 4': r'(?<=LEV)D(?=[^PEDQKR])',
    'caspase 5': r'(?<=[LW]EH)D',
    'caspase 6': r'(?<=VE[HI])D(?=[^PEDQKR])',
    'caspase 7': r'(?<=DEV)D(?=[^PEDQKR])',
    'caspase 8': r'(?<=[IL]ET)D(?=[^PEDQKR])',
    'caspase 9': r'(?<=LEH)D',
    'caspase 10': r'(?<=IEA)D',
    'chymotrypsin high specificity': r'([FY](?=[^P]))|(W(?=[^MP]))',
    'chymotrypsin low specificity': r'([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
    'clostripain': r'R',
    'cnbr': r'M',
    'enterokinase': r'(?<=[DE]{3})K',
    'factor xa': r'(?<=[AFGILTVM][DE]G)R',
    'formic acid': r'D',
    'glutamyl endopeptidase': r'[ED]',
    'granzyme b': r'(?<=IEP)D',
    'hydroxylamine': r'N(?=G)',
    'iodosobenzoic acid': r'W',
    'lysc': r'K',
    'ntcb': r'\w(?=C)',
    'pepsin ph1.3': r'((?<=[^HKR][^P])[^R](?=[FL][^P]))|((?<=[^HKR][^P])[FL](?=\w[^P]))',
    'pepsin ph2.0': r'((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
    'proline endopeptidase': r'(?<=[HKR])P(?=[^P])',
    'proteinase k': r'[AEFILTVWY]',
    'staphylococcal peptidase i': r'(?<=[^E])E',
    'thermolysin': r'[^DE](?=[AFILMV][^P])',
    'thrombin': r'((?<=G)R(?=G))|((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
    'trypsin': r'([KR](?=[^P]))',
    'trypsin_exception': r'((?<=[CD])K(?=D))|((?<=C)K(?=[HY]))|((?<=C)R(?=K))|((?<=R)R(?=[HR]))',
}

"""Best Combination Method and Method Helpers"""


def best_combination(lower_peptide_bound_entry, upper_peptide_bound_entry, allowed_missed_cleavages,
                     standard_protease_check_vars, additional_protease_check_vars, fasta_file_path_entry,
                     degree_entry, peptide_or_protein_lvl_var):
    """
    Display digestion result based on user inputs.

    This function calculates the protease combinations and their scores based on the selected finder level.
    It then displays the results to the user using a message box.

    Parameters:
    - lower_peptide_bound_entry (Entry): Entry widget containing lower peptide bound value.
    - upper_peptide_bound_entry (Entry): Entry widget containing upper peptide bound value.
    - allowed_missed_cleavages (Entry): Entry widget containing the number of allowed missed cleavages.
    - standard_protease_check_vars (list of IntVar): List of IntVar variables linked to standard protease checkboxes.
    - additional_protease_check_vars (list of IntVar): List of IntVar variables linked to additional protease checkboxes
    - fasta_file_path_entry (Entry): Entry widget containing the path to the selected FASTA/CSV file.
    - selected_finder_level (str): The selected protease finder level (radio button value).
    - degree_entry (Entry): Entry widget containing the degree of protease combinations (1, 2, 3, etc.).

    Returns:
    - None

    Usage:
    - best_combination(lowerPeptideBoundEntry, upperPeptideBoundEntry, allowedMissedCleavages,
                     standardProteaseCheckVars, additionalProteaseCheckVars, fastaFilePathEntry,
                     selectedFinderLevel, degreeEntry)

    Note: - This function uses the selected finder level to calculate protease combinations and displays the results
    using message boxes.
    - If no protease combinations are found, an error message is displayed.
    - This method connects the Find Best Combination button to the protease_finder methods.
    """
    # Determine the selected proteases based on the checked checkboxes.
    selected_standard_proteases = [proteaseList[idx] for idx, var in enumerate(standard_protease_check_vars) if
                                   var.get() == 1]
    selected_additional_proteases = [proteaseList[idx] for idx, var in enumerate(additional_protease_check_vars) if
                                     var.get() == 1]

    # Get the degree of protease combinations.
    degree = int(degree_entry.get())

    peptide_or_protein_lvl = peptide_or_protein_lvl_var.get()

    # Call the generalized protease finder function based on the selected degree.
    result = best_combination_helper(lower_peptide_bound_entry,
                                     upper_peptide_bound_entry,
                                     allowed_missed_cleavages,
                                     selected_standard_proteases,
                                     selected_additional_proteases,
                                     fasta_file_path_entry,
                                     degree,
                                     peptide_or_protein_lvl)

    result_file = result[0]
    top_results = result[1]

    # Display the results using message boxes.
    download_window = tk.Toplevel()
    download_window.title("Download Window")
    download_window.geometry("400x300")  # Adjusted height to accommodate the labels.

    label = Label(download_window, text="Top combinations:", font=("Roboto", 12))
    label.pack(pady=4)

    # Add labels displaying the first 5 items in top_results.
    for i in range(min(5, len(top_results))):
        result_text = f"{top_results[i]}"
        sample_label = Label(download_window, text=result_text, font=("Roboto", 12))
        sample_label.pack(pady=2)

    # Label to inform the user about downloading results.
    label = Label(download_window, text="Download Protease Combinations:", font=("Roboto", 12))
    label.pack(pady=10)

    # Button to initiate the download of the result file.
    download_button = Button(download_window, text="Download", command=lambda: download_csv_file(result_file),
                             font=('Roboto', 16), pady=4, padx=8)
    download_button.pack(pady=10)




def best_combination_helper(lower_peptide_bound_daltons_entry, higher_peptide_bound_daltons_entry,
                            allowed_missed_cleavages, selected_standard_proteases, selected_additional_proteases,
                            file_path_entry, degree: int, peptide_or_protein_lvl):
    """
    Finds the best protease combinations based on peptide coverage.

    Args:
        lower_peptide_bound_daltons_entry (tkinter.Entry): Entry widget for the lower bound of peptide mass in daltons.
        higher_peptide_bound_daltons_entry (tkinter.Entry): Entry widget for the upper bound of peptide mass in daltons.
        allowed_missed_cleavages (tkinter.Entry): Entry widget for the number of allowed missed cleavages.
        selected_standard_proteases (list of str): List of selected standard protease names.
        selected_additional_proteases (list of str): List of selected additional protease names.
        file_path_entry (tkinter.Entry): Entry widget for the path to the file containing protein sequences.
        degree (int): The degree of protease combinations to consider.
        peptide_or_protein_lvl (str): Determines whether the program calculates for peptide or protein coverage

    Returns:
        str: Path to the generated CSV file containing the best protease combinations and their peptide coverage.

    Notes:
        - The degree parameter specifies how many additional proteases to combine with the standard proteases.
        - A lower degree is recommended due to the exponential growth of combinations as the degree increases.
    """
    # Extract input values from GUI elements
    lower_bound = float(lower_peptide_bound_daltons_entry.get())
    upper_bound = float(higher_peptide_bound_daltons_entry.get())
    missed_cleavages_num = int(allowed_missed_cleavages.get())

    if peptide_or_protein_lvl == "peptide":
        proteins = file_to_list(file_path_entry.get())
        # Initialize coverage map with standard proteases
        standard_coverage_map = {}
        for protease in selected_standard_proteases:
            standard_coverage_map = protease_or(
                one_protease_digestion(protease, proteins, lower_bound, upper_bound, missed_cleavages_num),
                standard_coverage_map)

        # Generate protease combinations with the specified degree
        protease_combinations = generate_protease_combinations_peptide(selected_additional_proteases, degree,
                                                                       proteins, lower_bound, upper_bound,
                                                                       missed_cleavages_num, standard_coverage_map)
        protease_combinations_copy = copy.deepcopy(protease_combinations)
        protease_combinations.sort(key=lambda x: x[-1], reverse=True)

    elif peptide_or_protein_lvl == "protein":
        protein_list = file_to_list_2(file_path_entry.get())
        standard_identified_proteins = set()

        for protease in selected_standard_proteases:
            standard_identified_proteins = find_identified_proteins(protease, protein_list, lower_bound, upper_bound,
                                                                    missed_cleavages_num)[0]

        protease_combinations = generate_protease_combinations_protein(selected_additional_proteases, degree,
                                                                       protein_list, lower_bound,
                                                                       upper_bound, missed_cleavages_num,
                                                                       standard_identified_proteins)
        protease_combinations_copy = copy.deepcopy(protease_combinations)
        protease_combinations.sort(key=lambda x: x[-1], reverse=True)

    else:
        raise ValueError(
            f"Invalid radio button value, must be 'peptide' or 'protein'.")

    # Generate the CSV filename based on the degree of combinations
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    standard_protease_names = "_".join(selected_standard_proteases)
    output_csv_file = f'Combinations{degree}_{peptide_or_protein_lvl}_{standard_protease_names}_{timestamp}.csv'

    # Generate headers dynamically based on the degree
    headers = [f'Protease {i + 1}' for i in range(degree)] + ['Coverage']

    # Write the sorted list to a CSV file
    with open(output_csv_file, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(headers)
        csv_writer.writerows(protease_combinations_copy)

    return output_csv_file, protease_combinations


def generate_protease_combinations_peptide(selected_proteases, degree, proteins, lower_bound, upper_bound,
                                           missed_cleavages_num,
                                           standard_coverage_map=None):
    """
    Generate protease combinations and calculate their peptide coverage.

    Args:
        selected_proteases (list of str): List of selected protease names.
        degree (int): Degree of combinations (1, 2, 3, etc.).
        proteins (list of str): List of protein sequences.
        lower_bound (float): Lower bound of peptide mass in daltons.
        upper_bound (float): Upper bound of peptide mass in daltons.
        missed_cleavages_num (int): Number of allowed missed cleavages.
        standard_coverage_map (dict, optional): Existing coverage map to combine with. Default is None.

    Returns:
        list of tuples: List of protease combinations with their peptide coverage.
    """
    protease_combinations = []
    processed_combinations = set()

    if standard_coverage_map is None:
        standard_coverage_map = {}

    if degree == 1:
        # Special handling for single protease
        for protease in selected_proteases:
            combo = (protease,)
            if combo not in processed_combinations:
                processed_combinations.add(combo)
                temp_coverage_map = get_combined_coverage_map(combo, proteins, lower_bound, upper_bound,
                                                              missed_cleavages_num, standard_coverage_map)
                peptide_coverage = calculate_peptide_coverage(temp_coverage_map)
                protease_combinations.append((*combo, peptide_coverage))
                print(f"{(*combo, peptide_coverage)}")
    else:
        # General case for degrees greater than 1
        for combo in combinations(selected_proteases, degree):
            combo = tuple(sorted(combo))
            if combo not in processed_combinations:
                processed_combinations.add(combo)
                temp_coverage_map = get_combined_coverage_map(combo, proteins, lower_bound, upper_bound,
                                                              missed_cleavages_num, standard_coverage_map)
                peptide_coverage = calculate_peptide_coverage(temp_coverage_map)
                protease_combinations.append((*combo, peptide_coverage))
                print(f"{(*combo, peptide_coverage)}")

    return protease_combinations


def generate_protease_combinations_protein(selected_proteases, degree, protein_list, lower_bound, upper_bound,
                                           missed_cleavages_num, standard_protein_set=None):
    """
    Generate protease combinations and calculate their protein coverage.

    Args:
        selected_proteases (list of str): List of selected protease names.
        degree (int): Degree of combinations (1, 2, 3, etc.).
        protein_list (list of Tuple): List of protein sequences.
        lower_bound (float): Lower bound of peptide mass in daltons.
        upper_bound (float): Upper bound of peptide mass in daltons.
        missed_cleavages_num (int): Number of allowed missed cleavages.
        standard_protein_set (set, optional): Existing coverage map to combine with. Default is None.

    Returns:
        list of tuples: List of protease combinations with their protein coverage percentage.
    """

    protease_combinations = []
    processed_combinations = set()

    if standard_protein_set is None:
        standard_protein_set = set()

    if degree == 1:
        for protease in selected_proteases:
            combo = (protease,)
            if combo not in processed_combinations:
                processed_combinations.add(combo)
                protein_set = get_combined_protein_set(combo, protein_list, lower_bound, upper_bound,
                                                       missed_cleavages_num,
                                                       standard_protein_set)
                protein_coverage = round(len(protein_set) / len(protein_list), 6)
                protease_combinations.append((*combo, protein_coverage))
                print(f"{(*combo, protein_coverage)}")
    else:
        for combo in combinations(selected_proteases, degree):
            combo = tuple(sorted(combo))
            if combo not in processed_combinations:
                processed_combinations.add(combo)
                protein_set = get_combined_protein_set(combo, protein_list, lower_bound, upper_bound,
                                                       missed_cleavages_num,
                                                       standard_protein_set)
                protein_coverage = round(len(protein_set) / len(protein_list) * 100, 6)
                protease_combinations.append((*combo, protein_coverage))
                print(f"{(*combo, protein_coverage)}")

    return protease_combinations


"""Run Digestion Method and Method Helpers"""


def run_digestion(lower_peptide_bound_daltons_entry, higher_peptide_bound_daltons_entry, allowed_missed_cleavages,
                  protease_check_vars, fasta_file_path_entry):
    """
    Display a download window for digestion results.

    This function opens a new window that allows the user to download the results of the digestion process. It
    initiates the digestion process with provided parameters and selected proteases, then generates a downloadable
    result file.

    Parameters:
    - lower_peptide_bound_daltons_entry (Entry): Entry widget containing lower peptide bound value in Daltons.
    - higher_peptide_bound_daltons_entry (Entry): Entry widget containing higher peptide bound value in Daltons.
    - allowed_missed_cleavages (Entry): Entry widget containing the number of allowed missed cleavages.
    - protease_check_vars (list of IntVar): List of IntVar variables linked to protease checkboxes.
    - fasta_file_path_entry (Entry): Entry widget containing the path to the selected FASTA/CSV file.

    Returns:
    - None

    Usage:
    - run_digestion(lowerPeptideBoundDaltonsEntry, upperPeptideBoundDaltonsEntry, allowedMissedCleavages,
                       proteaseCheckVars, fastaFilePathEntry)

    Note:
    - This function opens a new window where the user can initiate the download of the digestion results.
    - This method connects the Run Digestion button from the GUI to the protease_runner method
    """
    # Determine the selected proteases based on the checked checkboxes.
    selected_proteases = [proteaseList[idx] for idx, var in enumerate(protease_check_vars) if var.get() == 1]

    # Run the run_digestion_helper function to generate the result file.
    output = run_digestion_helper(lower_peptide_bound_daltons_entry, higher_peptide_bound_daltons_entry,
                                  allowed_missed_cleavages, selected_proteases, fasta_file_path_entry)

    # Create a new download window.
    download_window = tk.Toplevel()
    download_window.title("Download Window")
    download_window.geometry("400x100")

    # Label to inform the user about downloading results.
    label = Label(download_window, text=f"Download Results of Digestion ({output[1]}):", font=("Roboto", 12))
    label.pack(pady=10)

    # Button to initiate the download of the result file.
    download_button = Button(download_window, text="Download", command=lambda: download_csv_file(output[0]),
                             font=('Roboto', 16), pady=4, padx=8)
    download_button.pack()


def run_digestion_helper(lower_peptide_bound_entry, upper_peptide_bound_entry, allowed_missed_cleavages,
                         selected_proteases, file_path_entry):
    """
    Makes the file that run_digestion() shows to the user to download

    Args:
        lower_peptide_bound_entry (tkinter.Entry): Entry widget for the lower bound of peptide mass in daltons.
        upper_peptide_bound_entry (tkinter.Entry): Entry widget for the upper bound of peptide mass in daltons.
        allowed_missed_cleavages (tkinter.Entry): Entry widget for the number of allowed missed cleavages.
        selected_proteases (list of str): List of selected protease names.
        file_path_entry (tkinter.Entry): Entry widget for the path to the file containing protein sequences.

    Returns:
        str: Filepath of the generated CSV output.

    Note: This function processes protease digestion and generates a CSV output containing information about peptide
    coverage and cleavage locations. It uses the provided input parameters to set the digestion parameters and read
    protein sequences. The function starts by performing protease digestion for the first selected protease and
    calculates peptide coverage for it. Then, it processes the cleavage locations and sequences obtained from the
    digestion and adds them to the CSV data. The process is repeated for the remaining proteases, with the coverage
    maps being combined using the protease OR operation. Finally, the function calculates overall peptide coverage,
    generates a timestamp and protease names to create the output CSV filename, and writes the CSV data to the file.
    The function returns the filepath of the generated CSV output.
    """
    # Extract input values from GUI elements
    lower_bound = float(lower_peptide_bound_entry.get())
    upper_bound = float(upper_peptide_bound_entry.get())
    missed_cleavages_num = int(allowed_missed_cleavages.get())
    protein_list = file_to_list_2(file_path_entry.get())

    str1_list = [protein_name for protein_name, _ in protein_list]

    # Initialize data structures for CSV data and protease processing
    csv_data = []
    peptide_lengths = []
    recognized_proteins = set()
    coverage_map = {}
    unique_peptides = 0

    # Process proteases and generate coverage map
    for protease in selected_proteases:
        proteins_coverage_map = one_protease_digestion(protease, str1_list, lower_bound, upper_bound,
                                                       missed_cleavages_num)
        coverage_map = protease_or(proteins_coverage_map, coverage_map)

        # Calculate coverage map and peptide coverage for the current protease
        detected_AA = 0
        total_AA = 0
        for protein_coverage_map in proteins_coverage_map:
            detected_AA += sum(1 for i in protein_coverage_map if i == 1)
            total_AA += len(protein_coverage_map)
        AA_coverage = round(detected_AA / total_AA * 100, 6)

        csv_data.append([protease, AA_coverage])
        print(f"{protease}, {AA_coverage}")
        print(proteins_coverage_map)

        results = find_identified_proteins(protease, protein_list, lower_bound, upper_bound, missed_cleavages_num,
                                           csv_data)
        recognized_proteins = recognized_proteins.union(results[0])
        unique_peptides = results[1]
        peptide_lengths = results[2]

    # Calculate coverage map and peptide coverage for all proteases
    detected_AA = 0
    total_AA = 0
    for protein_coverage_map in coverage_map:
        detected_AA += sum(1 for i in protein_coverage_map if i == 1)
        total_AA += len(protein_coverage_map)
    AA_coverage = round(detected_AA / total_AA * 100, 6)

    csv_data.append(['all proteases', AA_coverage])
    print(f"all proteases,{AA_coverage}")
    print(coverage_map)

    # Generate the CSV filename based on the current date and proteases
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    protease_names = '_'.join(selected_proteases)
    output_csv_file = f'{protease_names}_{timestamp}.csv'

    # Write the csv_data list to a CSV file
    with open(output_csv_file, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(
            ['Location', 'Peptide Sequence', round(len(recognized_proteins) / len(protein_list), 6), unique_peptides])
        csv_writer.writerows(csv_data)

    # Peptide length distribution
    length_ranges = [(1, 10), (11, 20), (21, 30), (31, 40), (41, 50),
                     (51, 60), (61, 70), (71, 80), (81, 90), (91, 100), (101, 110)]
    length_counts = {f"{start}-{end}": 0 for start, end in length_ranges}

    for length in peptide_lengths:
        for start, end in length_ranges:
            if start <= length <= end:
                length_counts[f"{start}-{end}"] += 1
                break

    print("Peptide Length Distribution:")
    for length_range, count in length_counts.items():
        print(f"{length_range}: {count}")

    # Return the filepath of the generated CSV output
    output = [output_csv_file, AA_coverage]
    return output


"""Find Proteins Method and Method Helpers"""


def find_proteins(lower_peptide_bound_daltons_entry, higher_peptide_bound_daltons_entry, allowed_missed_cleavages,
                  protease_check_vars, fasta_file_path_entry):
    """
    Display a download window showing which proteins out of the inputted protein list/file can be 100% covered by the
    selected proteases

    This function opens a new window that allows the user to download the results of the find proteins process. It
    initiates the process with provided parameters and selected proteases, then generates a downloadable
    result file.

    Parameters:
    - lower_peptide_bound_daltons_entry (Entry): Entry widget containing lower peptide bound value in Daltons.
    - higher_peptide_bound_daltons_entry (Entry): Entry widget containing higher peptide bound value in Daltons.
    - allowed_missed_cleavages (Entry): Entry widget containing the number of allowed missed cleavages.
    - protease_check_vars (list of IntVar): List of IntVar variables linked to protease checkboxes.
    - fasta_file_path_entry (Entry): Entry widget containing the path to the selected FASTA/CSV file.

    Returns:
    - None

    Usage:
    - find_proteins(lowerPeptideBoundDaltonsEntry, upperPeptideBoundDaltonsEntry, allowedMissedCleavages,
                    proteaseCheckVars, fastaFilePathEntry)

    Note:
    - This function opens a new window where the user can initiate the download of the digestion results.
    - This method connects the Find Proteins button from the GUI to the protease_runner method
    """
    # Determine the selected proteases based on the checked checkboxes.
    selected_proteases = [proteaseList[idx] for idx, var in enumerate(protease_check_vars) if var.get() == 1]

    # Run the find_proteins_helper function to generate the result file.
    output = find_proteins_helper(lower_peptide_bound_daltons_entry, higher_peptide_bound_daltons_entry,
                                  allowed_missed_cleavages, selected_proteases, fasta_file_path_entry)

    # Create a new download window.
    download_window = tk.Toplevel()
    download_window.title("Download Window")
    download_window.geometry("400x100")

    # Label to inform the user about downloading results.
    label = Label(download_window, text=f"Download({output[1]}):", font=("Roboto", 12))
    label.pack(pady=10)

    # Button to initiate the download of the result file.
    download_button = Button(download_window, text="Download", command=lambda: download_csv_file(output[0]),
                             font=('Roboto', 16), pady=4, padx=8)
    download_button.pack()


def find_proteins_helper(lower_peptide_bound_entry, upper_peptide_bound_entry, allowed_missed_cleavages,
                         selected_proteases, file_path_entry):
    """

    :param lower_peptide_bound_entry:
    :param upper_peptide_bound_entry:
    :param allowed_missed_cleavages:
    :param selected_proteases:
    :param file_path_entry:
    :return:
    """
    # Extract input values from GUI elements
    lower_bound = float(lower_peptide_bound_entry.get())
    upper_bound = float(upper_peptide_bound_entry.get())
    missed_cleavages_num = int(allowed_missed_cleavages.get())
    protein_list = file_to_list_2(file_path_entry.get())

    # Initialize data structures for CSV data and protease processing
    csv_data = []
    temp_proteases = selected_proteases[1:]
    num_proteins = 0

    for sequence, name in protein_list:
        # Generate singular coverage map for this protein using just the first protease
        protein = [sequence]
        coverage_map = one_protease_digestion(selected_proteases[0], protein, lower_bound, upper_bound,
                                              missed_cleavages_num)

        if temp_proteases:
            for protease in temp_proteases:
                coverage_map = protease_or(coverage_map,
                                           one_protease_digestion(protease, protein, lower_bound, upper_bound,
                                                                  missed_cleavages_num))
        detected_AA = 0
        total_AA = 0
        for protein_coverage_map in coverage_map:
            detected_AA += sum(1 for i in protein_coverage_map if i == 1)
            total_AA += len(protein_coverage_map)
        if detected_AA == total_AA:
            csv_data.append((name, sequence))
            num_proteins = num_proteins + 1

    # Generate the CSV filename based on the current date and proteases
    protease_names = '_'.join(selected_proteases)
    output_csv_file = f'output_{num_proteins}_{protease_names}.csv'

    # Write the csv_data list to a CSV file
    with open(output_csv_file, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['Name', 'Sequence', str(num_proteins)])
        csv_writer.writerows(csv_data)

    # Return the filepath of the generated CSV output
    output = [output_csv_file, num_proteins]
    return output


"""Compare Peptides Method and Method Helpers"""


def compare_peptides_2(lower_peptide_bound_daltons_entry, higher_peptide_bound_daltons_entry, allowed_missed_cleavages,
                       protease_check_vars, fasta_file_path_entry):
    """

    :param lower_peptide_bound_daltons_entry:
    :param higher_peptide_bound_daltons_entry:
    :param allowed_missed_cleavages:
    :param protease_check_vars:
    :param fasta_file_path_entry:
    :return:
    """
    selected_proteases = [proteaseList[idx] for idx, var in enumerate(protease_check_vars) if var.get() == 1]

    # Run the run_digestion_helper function to generate the result file.
    compare_peptides_helper_2(lower_peptide_bound_daltons_entry, higher_peptide_bound_daltons_entry,
                              allowed_missed_cleavages, selected_proteases, fasta_file_path_entry)


def compare_peptides_helper_2(lower_peptide_bound_entry, upper_peptide_bound_entry, allowed_missed_cleavages,
                              selected_proteases, file_path_entry):
    """

    :param lower_peptide_bound_entry: 
    :param upper_peptide_bound_entry: 
    :param allowed_missed_cleavages: 
    :param selected_proteases: 
    :param file_path_entry: 
    :return: 
    """""
    protease1 = [selected_proteases[0]]
    protease2 = [selected_proteases[1]]

    csv_file1 = run_digestion_helper(lower_peptide_bound_entry, upper_peptide_bound_entry, allowed_missed_cleavages,
                                     protease1, file_path_entry)[0]
    csv_file2 = run_digestion_helper(lower_peptide_bound_entry, upper_peptide_bound_entry, allowed_missed_cleavages,
                                     protease2, file_path_entry)[0]
    compare_files_2(csv_file1, csv_file2, "temp")


def compare_files_2(csv_file1, csv_file2, output_file):
    """

    :param csv_file1:
    :param csv_file2:
    :param output_file:
    :return:
    """
    # Load the CSV files into DataFrames
    df1 = pd.read_csv(csv_file1)
    df2 = pd.read_csv(csv_file2)

    # Assuming the peptide sequences are in a column named 'Peptide Sequence'
    peptides1 = set(df1['Peptide Sequence'])
    peptides2 = set(df2['Peptide Sequence'])

    # Find common peptides
    common_peptides = peptides1.intersection(peptides2)

    # Create a DataFrame for the common peptides
    common_df = pd.DataFrame(list(common_peptides), columns=['Peptide'])

    # Save the common peptides to a new CSV file
    common_df.to_csv(output_file, index=False)

    # Print the number of overlapping peptides
    print(f"Number of overlapping peptides: {len(common_peptides)}")
    print(f"{csv_file1} {len(df1) - 1}")
    print(f"{csv_file2} {len(df2) - 1}")


def compare_peptides_3(lower_peptide_bound_daltons_entry, higher_peptide_bound_daltons_entry, allowed_missed_cleavages,
                       protease_check_vars, fasta_file_path_entry):
    """

    :param lower_peptide_bound_daltons_entry:
    :param higher_peptide_bound_daltons_entry:
    :param allowed_missed_cleavages:
    :param protease_check_vars:
    :param fasta_file_path_entry:
    :return:
    """
    selected_proteases = [proteaseList[idx] for idx, var in enumerate(protease_check_vars) if var.get() == 1]

    # Run the run_digestion_helper function to generate the result file.
    output_file = compare_peptides_helper_3(lower_peptide_bound_daltons_entry, higher_peptide_bound_daltons_entry,
                                            allowed_missed_cleavages, selected_proteases, fasta_file_path_entry)

    # Create a new download window.
    download_window = tk.Toplevel()
    download_window.title("Download Window")
    download_window.geometry("400x100")

    # Label to inform the user about downloading results.
    label = Label(download_window, text=f"Download Overlap:", font=("Roboto", 12))
    label.pack(pady=10)

    # Button to initiate the download of the result file.
    download_button = Button(download_window, text="Download", command=lambda: download_csv_file(output_file),
                             font=('Roboto', 16), pady=4, padx=8)
    download_button.pack()


def compare_peptides_helper_3(lower_peptide_bound_entry, upper_peptide_bound_entry, allowed_missed_cleavages,
                              selected_proteases, file_path_entry):
    """

    :param lower_peptide_bound_entry:
    :param upper_peptide_bound_entry:
    :param allowed_missed_cleavages:
    :param selected_proteases:
    :param file_path_entry:
    :return:
    """
    protease1 = [selected_proteases[0]]
    protease2 = [selected_proteases[1]]
    protease3 = [selected_proteases[2]]

    csv_file1 = run_digestion_helper(lower_peptide_bound_entry, upper_peptide_bound_entry, allowed_missed_cleavages,
                                     protease1, file_path_entry)[0]
    csv_file2 = run_digestion_helper(lower_peptide_bound_entry, upper_peptide_bound_entry, allowed_missed_cleavages,
                                     protease2, file_path_entry)[0]
    csv_file3 = run_digestion_helper(lower_peptide_bound_entry, upper_peptide_bound_entry, allowed_missed_cleavages,
                                     protease3, file_path_entry)[0]

    # Generate the CSV filename based on the current date and proteases
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    protease_names = '_'.join(selected_proteases)
    output_file = f'overlap_{protease_names}_{timestamp}.csv'

    compare_files_3(csv_file1, csv_file2, csv_file3, output_file)
    return output_file


def compare_files_3(csv_file1, csv_file2, csv_file3, output_file):
    """
    Compare peptides from three CSV files and print overlaps.

    :param csv_file1: Path to the first CSV file
    :param csv_file2: Path to the second CSV file
    :param csv_file3: Path to the third CSV file
    :param output_file: Path to the output CSV file for common peptides
    """
    # Load the CSV files into DataFrames
    df1 = pd.read_csv(csv_file1, low_memory=False)
    df2 = pd.read_csv(csv_file2, low_memory=False)
    df3 = pd.read_csv(csv_file3, low_memory=False)

    # Assuming the peptide sequences are in a column named 'Peptide Sequence'
    peptides1 = set(df1['Peptide Sequence'])
    peptides2 = set(df2['Peptide Sequence'])
    peptides3 = set(df3['Peptide Sequence'])

    just1 = len(peptides1) - 1
    just2 = len(peptides2) - 1
    just3 = len(peptides3) - 1

    intersection12 = peptides1.intersection(peptides2)
    intersection13 = peptides1.intersection(peptides3)
    intersection23 = peptides2.intersection(peptides3)
    intersection123 = peptides1.intersection(peptides2).intersection(peptides3)

    unique1 = peptides1.difference(peptides2).difference(peptides3)
    unique2 = peptides2.difference(peptides1).difference(peptides3)
    unique3 = peptides3.difference(peptides2).difference(peptides1)

    intersection12without3 = intersection12.difference(peptides3)
    intersection13without2 = intersection13.difference(peptides2)
    intersection23without1 = intersection23.difference(peptides1)

    name1 = csv_file1.split('_')[0]
    name2 = csv_file2.split('_')[0]
    name3 = csv_file3.split('_')[0]

    with open(output_file, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerows([
            [name1, just1],  # Just peptides in peptides1, excluding overlaps
            [name2, just2],  # Just peptides in peptides2, excluding overlaps
            [name3, just3],  # Just peptides in peptides3, excluding overlaps
            [f"{name1} \\ {name2} \\ {name3}", len(unique1)],  # Unique peptides in peptides1
            [f"{name2} \\ {name1} \\ {name3}", len(unique2)],  # Unique peptides in peptides2
            [f"{name3} \\ {name1} \\ {name2}", len(unique3)],  # Unique peptides in peptides3
            [f"{name1} n {name2} n {name3}", len(intersection123)],  # Peptides in all three sets
            [f"{name1} n {name2} \\ {name3}", len(intersection12without3)],
            # Intersection of peptides1 and peptides2 excluding peptides3
            [f"{name1} n {name3} \\ {name2}", len(intersection13without2)],
            # Intersection of peptides1 and peptides3 excluding peptides2
            [f"{name2} n {name3} \\ {name1}", len(intersection23without1)],
            # Intersection of peptides2 and peptides3 excluding peptides1
            [f"{name1} n {name2}", len(intersection12)],  # Peptides in both peptides1 and peptides2
            [f"{name1} n {name3}", len(intersection13)],  # Peptides in both peptides1 and peptides3
            [f"{name2} n {name3}", len(intersection23)]  # Peptides in both peptides2 and peptides3
        ])


"""Boundary Finder Method(not used currently)"""


def boundary_finder(csv_file_path_entry, lower_peptide_bound_entry, upper_peptide_bound_entry):
    """
    Find peptide mass boundaries based on a sorted list of identified peptides.

    Args:
        csv_file_path_entry (tkinter.Entry): Entry widget for the path to the CSV file containing identified peptides.
        lower_peptide_bound_entry (tkinter.Entry): Entry widget for displaying the lower peptide mass boundary.
        upper_peptide_bound_entry (tkinter.Entry): Entry widget for displaying the upper peptide mass boundary.

    Note: This function reads a list of identified peptides from a CSV file and calculates peptide mass boundaries
    based on the sorted list of peptide masses. It first calculates the mass of each peptide and sorts the list
    based on mass. A certain percentage of peptides (2.5% from both ends) are trimmed to remove outliers. The
    function then calculates the lower and upper peptide mass boundaries using the masses of the trimmed peptides.
    The calculated boundaries are inserted into the respective entry widgets in the GUI for display.
    """
    # Read identified peptides from the CSV file
    identifiedPeptideList = file_to_list(csv_file_path_entry.get())

    # Calculate the mass of each peptide and sort them based on mass
    sortedPeptideList = sorted(identifiedPeptideList, key=lambda peptide: mass.calculate_mass(peptide, aa_comp=aa_comp))

    # Remove a certain percentage of peptides from both ends as outliers
    num_peptides_to_remove = int(round(len(sortedPeptideList) * 0.025))
    trimmedPeptideList = sortedPeptideList[num_peptides_to_remove:-num_peptides_to_remove]
    # Calculate lower and upper peptide mass boundaries based on trimmed peptides
    lower_boundary = round(mass.calculate_mass(trimmedPeptideList[0].translate(translation_table), aa_comp=aa_comp), 4)
    upper_boundary = round(mass.calculate_mass(trimmedPeptideList[-1].translate(translation_table), aa_comp=aa_comp), 4)

    # Insert calculated boundaries into GUI entry widgets for display
    lower_peptide_bound_entry.insert(0, str(lower_boundary))
    upper_peptide_bound_entry.insert(0, str(upper_boundary))


"""Basic Proteomic operations"""


def one_protease_digestion(protease, protein_list, lower_bound, upper_bound, missed_cleavages_num):
    """
    Perform digestion of proteins using a specific protease.

    Args:
        protease (str): The name of the protease used for digestion.
        protein_list (list of str): List of protein sequences to be digested.
        lower_bound (float): Lower bound of the desired peptide mass range.
        upper_bound (float): Upper bound of the desired peptide mass range.
        missed_cleavages_num (int): Maximum number of allowed missed cleavages.

    Returns:
        list of numpy arrays: A list of numpy arrays representing coverage maps for each protein sequence.
                             Each array contains binary values (0 or 1) indicating whether an amino acid
                             is covered by a valid peptide within the mass range.

    Note:
        The function utilizes cached results to optimize repeated computations for the same arguments.
        It generates coverage maps indicating which amino acids are covered by peptides after digestion.
    """
    # Check if the result is already cached
    cache_key = (protease, protein_list[0], lower_bound, upper_bound, missed_cleavages_num)
    if cache_key in one_protease_cache_peptide:
        return one_protease_cache_peptide[cache_key]

    # Initialize a list to store coverage maps for each protein
    coverage_map = []

    # Loop through each protein sequence for digestion
    for protein in protein_list:
        # Create an array to track coverage of each amino acid
        AA_coverage_map = np.zeros(len(protein), dtype=int)

        # Digest the protein using the specified protease
        for start, peptide in parser.icleave(protein, custom_expasy_rules[protease], missed_cleavages_num,
                                             regex=True):
            peptide = peptide.translate(translation_table)
            current_peptide_mass = mass.calculate_mass(sequence=peptide, aa_comp=aa_comp)
            # Check if the peptide mass is within the desired range
            if lower_bound < current_peptide_mass < upper_bound:
                # Mark the amino acids covered by the peptide
                AA_coverage_map[start:start + len(peptide)] = 1

        # Append the coverage map of the protein to the list
        coverage_map.append(AA_coverage_map)

    # Cache the result for future use
    one_protease_cache_peptide[cache_key] = coverage_map

    # Return the list of coverage maps
    return coverage_map


def multiple_proteases_digestion(proteases, protein_list, lower_bound, upper_bound, missed_cleavages_num):
    """
    Perform digestion of proteins using multiple proteases.

    Args:
        proteases (list of str): List of protease names used for digestion.
        protein_list (list of str): List of protein sequences to be digested.
        lower_bound (float): Lower bound of the desired peptide mass range.
        upper_bound (float): Upper bound of the desired peptide mass range.
        missed_cleavages_num (int): Maximum number of allowed missed cleavages.

    Returns:
        list of numpy arrays: A list of numpy arrays representing coverage maps for each protein sequence.
                             Each array contains binary values (0 or 1) indicating whether an amino acid
                             is covered by a valid peptide within the mass range.

    Note:
        The function utilizes cached results to optimize repeated computations for the same arguments.
        It generates coverage maps indicating which amino acids are covered by peptides after digestion.
    """
    # Check if the result is already cached
    cache_key = (tuple(proteases), protein_list[0], lower_bound, upper_bound, missed_cleavages_num)
    if cache_key in multi_protease_cache_peptide:
        return multi_protease_cache_peptide[cache_key]

    # Initialize a list to store coverage maps for each protein
    coverage_map = []

    # Loop through each protein sequence for digestion
    for protein in protein_list:
        # Create an array to track coverage of each amino acid
        AA_coverage_map = np.zeros(len(protein), dtype=int)

        # Digest the protein using each protease
        for protease in proteases:
            for start, peptide in parser.icleave(protein, custom_expasy_rules[protease], missed_cleavages_num,
                                                 regex=True):
                peptide = peptide.translate(translation_table)
                current_peptide_mass = mass.calculate_mass(sequence=peptide, aa_comp=aa_comp)
                # Check if the peptide mass is within the desired range
                if lower_bound < current_peptide_mass < upper_bound:
                    # Mark the amino acids covered by the peptide
                    AA_coverage_map[start:start + len(peptide)] = 1

        # Append the coverage map of the protein to the list
        coverage_map.append(AA_coverage_map)

    # Cache the result for future use
    multi_protease_cache_peptide[cache_key] = coverage_map

    # Return the list of coverage maps
    return coverage_map


def protease_or(*maps):
    """
    Combine multiple coverage maps from protease digestion simulations using the "OR" protease operation.

    Args:
        *maps (list of list of numpy arrays): Coverage maps from multiple protease digestion simulations.

    Returns:
        list of numpy arrays: Combined coverage map using an "OR" operation, where an amino acid is marked
                             as covered (1) if it's covered by peptides from any of the proteases.

    Note:
        The function takes multiple coverage maps generated from different protease digestion simulations
        and combines them using an "OR" operation. The combined coverage map indicates whether an amino acid
        is covered by peptides from any of the proteases. The resulting map represents the union of coverage
        from all provided proteases, ensuring that an amino acid is marked as covered (1) if it's covered
        by peptides from any of the proteases.
    """
    # Initialize a copy of the first coverage map for combining
    if not maps:
        raise ValueError("At least one coverage map is required.")

    combined_map = [row.copy() for row in maps[0]]

    # Loop through each additional map and combine using OR operation
    for coverage_map in maps[1:]:
        for x in range(len(coverage_map)):
            for y in range(len(coverage_map[x])):
                # If the current map covers an amino acid, mark it as covered in the combined map
                if coverage_map[x][y] == 1:
                    combined_map[x][y] = 1

    # Return the combined coverage map
    return combined_map


def file_to_list(file_path):
    """
    Read peptide sequences from a file and convert them into a list.

    Args:
        file_path (str): Path to the input file containing peptide sequences.

    Returns:
        list of str: List of peptide sequences extracted from the input file.

    Raises:
        ValueError: If the provided file format is unsupported.

    Note:
        This function reads peptide sequences from an input file and creates a list of these sequences.
        The function supports two file formats: CSV and FASTA. If the file has a '.csv' extension, it's
        assumed that each row contains a protein sequence in the first column and the . If the file has a '.fasta'
        extension, it's assumed to be in FASTA format, and peptide sequences are extracted from the records.
        The resulting list contains all the extracted peptide sequences. If the provided file format is
        unsupported, a ValueError is raised.
    """
    peptide_list = []

    # if csv file
    if file_path.endswith('.csv'):
        # Read peptide sequences from a CSV file
        with open(file_path, 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            for row in csvreader:
                peptide_list.append(row[0])  # Assuming the peptide sequence is in the first column

    # if fasta file
    elif file_path.endswith('.fasta'):
        # Read peptide sequences from a FASTA file
        for record in SeqIO.parse(file_path, "fasta"):
            peptide_list.append(str(record.seq))

    else:
        raise ValueError("Unsupported file format")

    return peptide_list


def file_to_list_2(file_path):
    """
    Read protein sequences from a file and convert them into a list including both the protein sequence and name.

    Args:
        file_path (str): Path to the input file containing peptide sequences.

    Returns:
        list of str: List of peptide sequences extracted from the input file.

    Raises:
        ValueError: If the provided file format is unsupported.

    Note:
        This function reads peptide sequences from an input file and creates a list of these sequences.
        The function supports two file formats: CSV and FASTA. If the file has a '.csv' extension, it's
        assumed that each row contains a protein sequence in the first column and the protein name in the second.
        If the file has a '.fasta' extension, it's assumed to be in FASTA format, and peptide sequences are extracted
        from the records. The resulting list contains all the extracted peptide sequences. If the provided file format
        is unsupported, a ValueError is raised.
    """
    peptide_list: List[Tuple[str, str]] = []

    # if csv file
    if file_path.endswith('.csv'):
        # Read peptide sequences from a CSV file
        with open(file_path, 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            for row in csvreader:
                peptide_list.append((row[0], row[1]))  # Assuming the peptide sequence is in the first column

    # if fasta file
    elif file_path.endswith('.fasta'):
        # Read peptide sequences from a FASTA file
        for record in SeqIO.parse(file_path, "fasta"):
            peptide_list.append((str(record.seq), str(record.id)))

    else:
        raise ValueError("Unsupported file format")

    return peptide_list


def find_identified_proteins(protease, protein_list, lower_bound, upper_bound, missed_cleavages_num, file_data=None):
    """
        Identifies proteins based on peptide mass and coverage from a given list of proteins.

        Args:
            protease (str): The name of the protease to use for digestion.
            protein_list (list of tuples): List of tuples where each tuple contains a protein sequence and its name.
            lower_bound (float): Lower bound of peptide mass in daltons.
            upper_bound (float): Upper bound of peptide mass in daltons.
            missed_cleavages_num (int): Number of allowed missed cleavages.
            file_data (list, optional): List to append location and sequence of peptides if provided.

        Returns:
            tuple: A tuple containing:
                - A set of recognized protein names.
                - The number of peptides counted.
                - A list of peptide lengths.
    """

    cache_key = (protease, protein_list[0], lower_bound, upper_bound, missed_cleavages_num)
    if cache_key in protease_digestion_cache_protein:
        return protease_digestion_cache_protein[cache_key]

    recognized_proteins = set()
    peptide_counter = Counter()
    protein_scores = Counter()
    peptide_lengths = []

    # Cache the results
    results_cache = defaultdict(list)

    # Run parser.icleave once and cache the results
    for protein, protein_name in protein_list:
        results_cache[protein] = list(
            parser.icleave(protein, custom_expasy_rules[protease], missed_cleavages_num, regex=True))

    # First loop using cached results
    for protein, protein_name in protein_list:
        for location, seq in results_cache[protein]:
            seq = seq.translate(translation_table)
            current_peptide_mass = mass.calculate_mass(seq, aa_comp=aa_comp)
            # Check if the peptide mass is within the desired range
            if lower_bound < current_peptide_mass < upper_bound:
                peptide_counter.update([seq])
                if file_data is not None:
                    file_data.append([location, seq])
                peptide_lengths.append(len(seq))

    # Second loop using cached results
    for protein, protein_name in protein_list:
        for location, seq in results_cache[protein]:
            if len(seq) >= 8 and peptide_counter[seq] == 1:
                protein_scores.update([protein_name])

    for protein_name, score in protein_scores.items():
        recognized_proteins.add(protein_name)

    protease_digestion_cache_protein[cache_key] = recognized_proteins, len(peptide_counter), peptide_lengths

    return recognized_proteins, len(peptide_counter), peptide_lengths


def calculate_peptide_coverage(coverage_map):
    detected_AA = 0
    total_AA = 0
    for AA_coverage_map in coverage_map:
        detected_AA += sum(1 for i in AA_coverage_map if i == 1)
        total_AA += len(AA_coverage_map)
    return round(detected_AA / total_AA * 100, 6)


def get_combined_coverage_map(proteases, proteins, lower_bound, upper_bound, missed_cleavages_num,
                              standard_coverage_map=None):
    protease_maps = [one_protease_digestion(protease, proteins, lower_bound, upper_bound, missed_cleavages_num) for
                     protease in proteases]
    if standard_coverage_map is None:
        coverage_map = protease_or(*protease_maps)
    else:
        coverage_map = protease_or(*protease_maps, standard_coverage_map)
    return coverage_map


def get_combined_protein_set(proteases, protein_list, lower_bound, upper_bound, missed_cleavages_num,
                             standard_protein_set=None):
    if standard_protein_set is None:
        protein_set = set()
    else:
        protein_set = standard_protein_set
    for protease in proteases:
        digestion_results = find_identified_proteins(protease, protein_list, lower_bound, upper_bound,
                                                     missed_cleavages_num)
        protein_set = protein_set.union(digestion_results[0])
    return protein_set  # , digestion_results[1]


def clone(list1):
    """
    Create a deep copy of a list.

    Args:
        list1 (list): The original list to be cloned.

    Returns:
        list: A deep copy of the input list.

    Note: This function creates a new list that is a deep copy of the input list. Unlike a shallow copy, a deep copy
    creates new copies of the objects within the original list. This ensures that changes made to the objects in the
    cloned list do not affect the original list and vice versa.
    """
    list_copy = list(list1)
    return list_copy


"""UI Methods"""


def download_csv_file(output_csv_file):
    """
    Downloads a CSV file using a file dialog.

    This function allows the user to choose a location to save a CSV file by opening a file dialog window.
    The function takes the path of the CSV file that needs to be downloaded and copied to the selected location.

    Parameters:
    - output_csv_file (str): The path of the CSV file to be downloaded.

    Returns:
    - None

    Usage:
    - download_csv_file("path/to/output.csv")

    """
    # Open a file dialog window to prompt the user to select a location and name for the downloaded CSV file.
    file_path = filedialog.asksaveasfilename(
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv")],
        title="Save CSV File",
        initialfile=output_csv_file.split("/")[-1],  # Provide a default file name based on the input path.
    )

    # Check if a valid file path was selected in the dialog.
    if file_path:
        # Copy the original CSV file to the selected location.
        shutil.copy(output_csv_file, file_path)


def browse_file(entry_widget):
    """
    Opens a file dialog to browse and select a file.

    This function displays a file dialog window that allows the user to browse and select a file.
    The selected file's path is then inserted into the provided entry widget.

    Parameters:
    - entry_widget (Entry): The Tkinter Entry widget where the selected file path will be inserted.

    Returns:
    - None

    Usage:
    - browse_file(fastaFilePathEntry)

    Note: - The file dialog allows the user to choose from various file types (FASTA and CSV) as specified in the
    filetypes parameter.

    """
    # Open a file dialog window to browse and select a file.
    file_path = filedialog.askopenfilename(
        filetypes=[("FASTA files", "*.fasta")])

    # Check if a valid file path was selected in the dialog.
    if file_path:
        # Clear the existing entry widget content and insert the selected file path.
        entry_widget.delete(0, tk.END)
        entry_widget.insert(0, file_path)


def _on_mousewheel(event, canvas):
    """
    Adjusts the scrolling of a canvas when the mousewheel is used.

    This function is intended to be used as a callback for mousewheel events on a canvas widget. It calculates the
    scrolling amount based on the mousewheel delta and adjusts the vertical view of the canvas accordingly.

    Parameters:
    - event (Event): The mousewheel event triggered by the user.
    - canvas (Canvas): The canvas widget to be scrolled.

    Returns:
    - None

    Usage:
    - canvas.bind("<MouseWheel>", lambda event: _on_mousewheel(event, canvas))

    Note:
    - This function is meant to be used internally and might not be directly called by the user.

    """
    # Calculate the scrolling amount based on the mousewheel delta and adjust the canvas view.
    # The scaling factor 1.5 is used to control the scroll speed.
    canvas.yview_scroll(-1 * int((event.delta / 120) * 1.5), "units")


def select_all_proteases(protease_check_vars):
    """
    Selects all protease checkboxes.

    This function is designed to be used when a "Select All" checkbox is clicked to mark all protease checkboxes as
    checked. It iterates through a list of Tkinter IntVar variables associated with protease checkboxes and sets
    their values to 1 (checked).

    Parameters:
    - protease_check_vars (list of IntVar): A list of IntVar variables linked to protease checkboxes.

    Returns:
    - None

    Usage: - protease_vars = [IntVar() for _ in range(number_of_proteases)] - # Initialize checkboxes and associate
    their IntVar variables. - select_all_checkbox = Checkbutton(root, text="Select All", command=lambda:
    select_all_proteases(protease_vars)) - protease_checkboxes = [Checkbutton(root, text=protease_name, variable=var)
    for protease_name, var in zip(protease_names, protease_vars)]

    Note: - This function assumes that the provided IntVar variables follow the Tkinter convention (0 for unchecked,
    1 for checked).

    """
    # Iterate through the list of IntVar variables and set each one to 1 (checked).
    for var in protease_check_vars:
        var.set(1)


def select_all_high_specificity_proteases(protease_check_vars):
    """
    Selects all high specificity proteases based on the protease names provided.

    Args:
        protease_check_vars (list of IntVar): Dictionary where keys are protease names and values are Tkinter variables
                                    representing whether the protease is selected (1) or not (0).
    """
    checklist_list = [1, 1, 1, 1, 1, 1, 1, 1,
                      1, 1, 1, 1, 1, 0,
                      0, 1, 1, 1, 1, 1,
                      1, 1, 1, 1, 1, 1,
                      1, 0, 1, 0, 1,
                      0, 1, 1]

    for var in protease_check_vars:
        var.set(checklist_list[0])
        checklist_list = checklist_list[1:]


"""Main UI"""


def create_tab(tab_parent, radio_button_string_var):
    """
    Create the contents of Tab in a GUI application.
    """
    # Labels and Entry fields for setting lower and upper mass spectrometer bounds
    lowerBoundText = Label(tab_parent, text="Lower Bound(Da) of Mass Spectrometer:", pady=5, padx=10,
                           font=("Roboto", 10), bg="#ededed")
    lowerPeptideBoundDaltonsEntry = Entry(tab_parent, font=("Roboto", 12))
    lowerPeptideBoundDaltonsEntry.insert(0, '400')

    higherBoundText = Label(tab_parent, text="Upper Bound(Da) of Mass Spectrometer:", pady=5, padx=10,
                            font=("Roboto", 10), bg="#ededed")
    upperPeptideBoundDaltonsEntry = Entry(tab_parent, font=("Roboto", 12))
    upperPeptideBoundDaltonsEntry.insert(0, '6000')

    missedCleavagesText = Label(tab_parent, text="Number of Allowed Missed Cleavages:", pady=5, padx=10,
                                font=("Roboto", 10), bg="#ededed")
    allowedMissedCleavages = Entry(tab_parent, font=("Roboto", 12))
    allowedMissedCleavages.insert(0, '2')

    # Standard Proteases checklist
    standard_protease_label = Label(tab_parent, text="Standard Proteases:", pady=4, padx=10,
                                    font=("Roboto", 10), bg="#ededed")

    standard_protease_frame = Frame(tab_parent, bg="#ededed", width=250)
    standard_canvas = Canvas(standard_protease_frame, bg="#ededed", height=150, width=250)
    standard_scrollbar = Scrollbar(standard_protease_frame, orient="vertical", command=standard_canvas.yview)
    standard_canvas.config(yscrollcommand=standard_scrollbar.set)
    standard_scrollbar.pack(side="right", fill="y")
    standard_canvas.pack(side="left", fill="both", expand=True)

    standard_protease_inner_frame = Frame(standard_canvas, bg="#ededed")
    standard_canvas.create_window((0, 0), window=standard_protease_inner_frame, anchor="nw")

    standard_protease_check_vars = [IntVar(value=0) for _ in range(len(proteaseList))]

    for idx, protease in enumerate(proteaseList):
        protease_check = Checkbutton(standard_protease_inner_frame, text=protease,
                                     variable=standard_protease_check_vars[idx],
                                     font=("Roboto", 10), bg="#ededed")
        protease_check.pack(anchor=W)

    standard_protease_inner_frame.update_idletasks()
    standard_canvas.config(scrollregion=standard_canvas.bbox("all"))
    standard_canvas.bind("<MouseWheel>", lambda event: _on_mousewheel(event, standard_canvas))

    # Additional Proteases checklist
    additional_protease_label = Label(tab_parent, text="Additional Proteases:", pady=4, padx=10,
                                      font=("Roboto", 10), bg="#ededed")

    additional_protease_frame = Frame(tab_parent, bg="#ededed", width=250)
    additional_canvas = Canvas(additional_protease_frame, bg="#ededed", height=150, width=250)
    additional_scrollbar = Scrollbar(additional_protease_frame, orient="vertical", command=additional_canvas.yview)
    additional_canvas.config(yscrollcommand=additional_scrollbar.set)
    additional_scrollbar.pack(side="right", fill="y")
    additional_canvas.pack(side="left", fill="both", expand=True)

    additional_protease_inner_frame = Frame(additional_canvas, bg="#ededed")
    additional_canvas.create_window((0, 0), window=additional_protease_inner_frame, anchor="nw")

    additional_protease_check_vars = [IntVar(value=0) for _ in range(len(proteaseList))]

    for idx, protease in enumerate(proteaseList):
        protease_check = Checkbutton(additional_protease_inner_frame, text=protease,
                                     variable=additional_protease_check_vars[idx],
                                     font=("Roboto", 10), bg="#ededed")
        protease_check.pack(anchor=W)

    additional_protease_inner_frame.update_idletasks()
    additional_canvas.config(scrollregion=additional_canvas.bbox("all"))
    additional_canvas.bind("<MouseWheel>", lambda event: _on_mousewheel(event, additional_canvas))

    select_all = Checkbutton(tab_parent, text="Select All Proteases",
                             command=lambda: select_all_proteases(additional_protease_check_vars))

    select_all_HS = Checkbutton(tab_parent, text="Select High Specificity Proteases",
                                command=lambda: select_all_high_specificity_proteases(additional_protease_check_vars))

    # Labels, entries, and buttons for setting bounds, file input, and starting processes
    fastaFileLabel = Label(tab_parent, text="Select or Paste a .fasta File:",
                           pady=8, font=("Roboto", 10), bg="#ededed")
    fastaFilePathEntry = Entry(tab_parent, font=("Roboto", 12))
    fastaFileBrowseButton = Button(tab_parent, text="Browse", command=lambda: browse_file(fastaFilePathEntry),
                                   font=("Roboto", 10))

    best_combo_button = Button(tab_parent, text="Find Best Combination",
                               command=lambda: best_combination(lowerPeptideBoundDaltonsEntry,
                                                                upperPeptideBoundDaltonsEntry,
                                                                allowedMissedCleavages,
                                                                standard_protease_check_vars,
                                                                additional_protease_check_vars,
                                                                fastaFilePathEntry,
                                                                degree_entry,
                                                                radio_button_string_var),
                               font=("Roboto", 11), pady=3, padx=4)

    # Degree Entry and Label
    degree_frame = Frame(tab_parent, bg="#ededed")
    degree_label = Label(degree_frame, text="Degree of Combinations:", pady=10, padx=0,
                         font=("Roboto", 10), bg="#ededed")
    degree_entry = Entry(degree_frame, font=("Roboto", 12), width=1)
    degree_entry.insert(0, '2')
    degree_label.grid(row=0, column=0, padx=3, pady=0, sticky="e")
    degree_entry.grid(row=0, column=1, padx=3, pady=0, sticky="w")

    # Radio buttons for protease finder level selection
    radio_frame = Frame(tab_parent, bg="#ededed")

    # Create radio buttons
    protein_level_coverage_radio = Radiobutton(radio_frame, text="Protein Level Coverage",
                                               variable=radio_button_string_var, value="protein")
    protein_level_coverage_radio.pack(side="left", padx=5)  # Pack to the left

    peptide_level_coverage_radio = Radiobutton(radio_frame, text="Peptide Level Coverage",
                                               variable=radio_button_string_var, value="peptide")
    peptide_level_coverage_radio.pack(side="left", padx=5)  # Pack to the left

    # Set the default selection
    radio_button_string_var.set("peptide")  # Default selection

    # Frame for buttons to keep them compact
    button_frame = Frame(tab_parent, bg="#ededed")

    startProgram = Button(button_frame, text="Find Proteins",
                          command=lambda: find_proteins(lowerPeptideBoundDaltonsEntry,
                                                        upperPeptideBoundDaltonsEntry,
                                                        allowedMissedCleavages,
                                                        standard_protease_check_vars,
                                                        fastaFilePathEntry),
                          font=("Roboto", 11), pady=3, padx=6)
    startProgram.grid(row=0, column=0, padx=5, pady=5)

    startProgram2 = Button(button_frame, text="Run Digestion",
                           command=lambda: run_digestion(lowerPeptideBoundDaltonsEntry,
                                                         upperPeptideBoundDaltonsEntry,
                                                         allowedMissedCleavages,
                                                         standard_protease_check_vars,
                                                         fastaFilePathEntry),
                           font=("Roboto", 11), pady=3, padx=4)
    startProgram2.grid(row=0, column=2, padx=5, pady=5)

    startProgram3 = Button(button_frame, text="Find Peptide Overlap",
                           command=lambda: compare_peptides_3(lowerPeptideBoundDaltonsEntry,
                                                              upperPeptideBoundDaltonsEntry,
                                                              allowedMissedCleavages,
                                                              standard_protease_check_vars,
                                                              fastaFilePathEntry),
                           font=("Roboto", 11), pady=3, padx=4)
    startProgram3.grid(row=0, column=1, padx=5, pady=5)

    """
    boundaryFinderButton = Button(button_frame, text="Find Boundary",
                                  command=lambda: boundary_finder(fastaFilePathEntry,
                                                                  lowerPeptideBoundDaltonsEntry,
                                                                  upperPeptideBoundDaltonsEntry),
                                  font=("Roboto", 11), pady=3, padx=6)
    # boundaryFinderButton.grid(row=0, column=3, padx=5, pady=5)
    """

    # Grid layout configuration
    lowerBoundText.grid(row=0, column=0, columnspan=2, padx=10, pady=5)
    lowerPeptideBoundDaltonsEntry.grid(row=1, column=0, columnspan=2, padx=10, pady=0)
    higherBoundText.grid(row=2, column=0, columnspan=2, padx=10, pady=5)
    upperPeptideBoundDaltonsEntry.grid(row=3, column=0, columnspan=2, padx=10, pady=0)
    missedCleavagesText.grid(row=4, column=0, columnspan=2, padx=10, pady=5)
    allowedMissedCleavages.grid(row=5, column=0, columnspan=2, padx=10, pady=0)

    standard_protease_label.grid(row=6, column=0, padx=10, pady=10)
    additional_protease_label.grid(row=6, column=1, padx=10, pady=10)
    standard_protease_frame.grid(row=7, column=0, padx=10, pady=0)
    additional_protease_frame.grid(row=7, column=1, padx=10, pady=0)

    select_all.grid(row=8, column=0)
    select_all_HS.grid(row=8, column=1)

    fastaFileLabel.grid(row=9, column=0, columnspan=2, padx=10, pady=5)
    fastaFilePathEntry.grid(row=10, column=0, columnspan=2, padx=10, pady=5, sticky="ew")
    fastaFileBrowseButton.grid(row=11, column=0, columnspan=2, padx=10, pady=0)

    degree_frame.grid(row=12, column=0, columnspan=2, padx=10, pady=0)

    radio_frame.grid(row=13, column=0, columnspan=2, padx=10, pady=3)

    best_combo_button.grid(row=14, column=0, columnspan=2, pady=0)

    button_frame.grid(row=15, column=0, columnspan=2, pady=5)


if __name__ == '__main__':
    master = Tk()
    radio_button_str = StringVar()  # Create a persistent StringVar
    create_tab(master, radio_button_str)
    master.mainloop()
