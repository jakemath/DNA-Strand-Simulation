//
//  strand.cpp
//  dna
//
//  Created by Jacob Mathai on 7/15/19.
//  Copyright © 2019 Jacob Mathai. All rights reserved.
//

#include "strand.h"

// Function to read txt file data to map codons to corresponding amino acids
map<string, string> create_lookup()
{
    map<string, string> lookup;
    std::ifstream in("codon.txt");
    string line;
    while (getline(in, line))
    {
        std::istringstream iss(line);
        vector<string> results((std::istream_iterator<string>(iss)), std::istream_iterator<string>());
        lookup[results[0]] = results[1];
    }
    return lookup;
}

// Construct random string of bases and compute basic statistics on contents
strand::strand(unsigned long long size, string type_, map<string, string> codon_lookup)
{
    // Initialize counts and values
    type = type_;
    char bases[4] = {'A',' ','G','C'};  // Array of bases
    if (type == "DNA")
        bases[1] = 'T'; // DNA uses Thymine
    else
        bases[1] = 'U'; // RNA uses Uracil
    for (short i = 0; i < 4; ++i) // Initialize count of each base to 0
        base_counts[bases[i]] = 0;
    map<string, string>::iterator i;    // Initialize count of codons and amino acids to 0
    for (i = codon_lookup.begin(); i != codon_lookup.end(); ++i)
    {
        codon_counts[i -> first] = 0;
        acid_counts[i -> second] = 0;
    }
    // Initialize basic statistics to garbage values
    total_codons_acids = 0;
    max_base_proportion = pair<char, double>('Z', 0.0); // use std::pair to associate
    min_base_proportion = pair<char, double>('Z', 101.0); // base, codon or acid to the
    max_codon_proportion = pair<string, double>("ZZZ", 0.0); // statistic
    min_codon_proportion = pair<string, double>("ZZZ", 101.0);
    max_acid_proportion = max_codon_proportion;
    min_acid_proportion = min_codon_proportion;
    
    // Construct strand
    sequence = "";
    std::random_device rd;
    if (rd() % 100 < 50)    // Set 50% chance of an A/T or A/U
        sequence.push_back(bases[rd() % 2]);
    else
        sequence.push_back(bases[3 - rd() % 2]); // 50% chance of G/C
    ++base_counts[sequence.back()]; // Increment count of added base
    if (rd() % 100 < 50)    // Repeat
        sequence.push_back(bases[rd() % 2]);
    else
        sequence.push_back(bases[3 - rd() % 2]);
    ++base_counts[sequence.back()];
    string temp;    // Hold temporary substring of previous 3 bases to determine codon
    for (unsigned long long i = 3; i <= size; ++i)  // Enter loop after first two insertions 
    {                                               // since no codons possible with 2 bases
        if (rd() % 100 < 50)
            sequence.push_back(bases[rd() % 2]);
        else
            sequence.push_back(bases[3 - rd() % 2]);
        ++base_counts[sequence.back()];
        temp = sequence.substr(sequence.size() - 3, 3);
        ++codon_counts[temp];   // Increment count of codon
        ++acid_counts[codon_lookup[temp]]; // Increment count of acid
        ++total_codons_acids;   // Increment count of total codons/acids (will always be the same)
    }
    
    // Compute proportions from computed counts and find maxima/minima
    double temp_;
    map<char, unsigned long long>::iterator j = base_counts.begin();
    for (; j != base_counts.end(); ++j) // Optima of base proportions
    {
        temp_ = ((1.0 * j -> second) / sequence.size()) * 100.0;
        base_proportions[j -> first] = temp_; // Assign proportion value for base
        if (temp_ > max_base_proportion.second) // Update max base proportion
        {
            max_base_proportion.first = j -> first;
            max_base_proportion.second = temp_;
        }
        else if (temp_ < min_base_proportion.second) // Update min base proportion
        {
            min_base_proportion.first = j -> first;
            min_base_proportion.second = temp_;
        }
    }
    
    // Optima of codon proportions
    map<string, unsigned long long>::iterator k = codon_counts.begin();
    for (; k != codon_counts.end(); ++k)
    {
        temp_ = ((1.0 * k -> second) / total_codons_acids) * 100.0;
        codon_proportions[k -> first] = temp_; // Assign proportion value for codon
        if (temp_ > max_codon_proportion.second)    // Update max codon proportion
        {
            max_codon_proportion.first = k -> first;
            max_codon_proportion.second = temp_;
        }
        else if (temp_ < min_codon_proportion.second)   // Update min codon proportion
        {
            min_codon_proportion.first = k -> first;
            min_codon_proportion.second = temp_;
        }
    }
    
    // Optima of acid proportions. Multiple codons correspond to the same amino acid
    // hence this is a different statistic from the previous
    for (k = acid_counts.begin(); k != acid_counts.end(); ++k)
    {
        temp_ = ((1.0 * k -> second) / total_codons_acids) * 100.0;
        acid_proportions[k -> first] = temp_;   // Assign proportion value for codon
        if (temp_ > max_acid_proportion.second) // Update max acid proportion
        {
            max_acid_proportion.first = k -> first;
            max_acid_proportion.second = temp_;
        }
        else if (temp_ < min_acid_proportion.second) // Update min acid proportion
        {
            min_acid_proportion.first = k -> first;
            min_acid_proportion.second = temp_;
        }
    }
}

// Print base counts and proportions
void strand::print_base_data()
{
    cout << "BASE COUNTS & PROPORTIONS " << endl << endl;
    map<char, unsigned long long>::iterator i = base_counts.begin();
    for (; i != base_counts.end(); ++i)
    {
        cout << "\t" << i -> first << " - " << i -> second << " count, ";
        cout << base_proportions[i -> first] << "%" << endl;
    }
}

// Print codon counts, proportions, and optima
void strand::print_codon_data(map<string, string> codon_lookup)
{
    cout << "CODON COUNTS & PROPORTIONS: " << endl << endl;
    map<string, unsigned long long>::iterator i = codon_counts.begin();
    for (; i != codon_counts.end(); ++i)
    {
        cout << "\t" << i -> first << " " << codon_lookup[i -> first] << " - ";
        cout << i -> second << " count, " << codon_proportions[i -> first] << "%" << endl;
    }
    cout << endl << "MOST COMMON CODON: " << max_codon_proportion.first << " ";
    cout << codon_lookup[max_codon_proportion.first] << " - ";
    cout << codon_counts[max_codon_proportion.first] << " count, ";
    cout << max_codon_proportion.second << "%";
    
    cout << endl << "LEAST COMMON CODON: " << min_codon_proportion.first;
    cout << codon_lookup[min_codon_proportion.first] << " - ";
    cout << codon_counts[min_codon_proportion.first] << " count, ";
    cout << min_codon_proportion.second << "%" << endl;
}

// Print acid counts, proportions, and optima
void strand::print_acid_data(map<string, string> codon_lookup)
{
    cout << "AMINO ACID COUNTS & PROPORTIONS: " << endl << endl;
    map<string, unsigned long long>::iterator i = acid_counts.begin();
    for (; i != acid_counts.end(); ++i)
    {
        cout << "\t" << i -> first << " - " << i -> second << " count, ";
        cout << acid_proportions[i -> first] << "%" << endl;
    }
    cout << endl << "MOST COMMON AMINO: " << max_acid_proportion.first << " - ";
    cout << acid_counts[max_acid_proportion.first] << " count, ";
    cout << max_acid_proportion.second << "%";
    
    cout << endl << "LEAST COMMON AMINO ACID: " << min_acid_proportion.first << " - ";
    cout << acid_counts[min_acid_proportion.first] << " count, ";
    cout << min_acid_proportion.second << "%" << endl;
}

// Print all the data at once
void strand::print_all_data(map<string, string> codon_lookup)
{
    cout << "SEQUENCE LENGTH = " << sequence.size() << endl << endl;
    print_base_data();
    cout << endl;
    print_codon_data(codon_lookup);
    cout << endl;
    print_acid_data(codon_lookup);
}
