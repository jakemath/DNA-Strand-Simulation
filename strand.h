//
//  strand.h
//  dna
//
//  Created by Jacob Mathai on 7/15/19.
//  Copyright © 2019 Jacob Mathai. All rights reserved.
//

#ifndef strand_h
#define strand_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <vector>
#include <map>
#include <string>

using std::cout;
using std::endl;
using std::vector;
using std::map;
using std::string;
using std::pair;

class strand {
public:
    strand();
    strand(unsigned long long size, string type_, map<string,string> codon_lookup);
    
    // Data printing functions
    void print_base_data();
    void print_codon_data(map<string, string> codon_lookup);
    void print_acid_data(map<string, string> codon_lookup);
    void print_all_data(map<string, string> codon_lookup);
    friend std::ostream& operator << (std::ostream& o, strand s) { o << s.sequence; return o; }
    
    // Simple accessors for data
    map<char, unsigned long long> get_base_counts() { return base_counts; }
    map<char, double> get_base_proportions() { return base_proportions; }
    pair<char, double> get_min_base_proportion() { return min_base_proportion; }
    pair<char, double> get_max_base_proportion() { return max_base_proportion; }
    
    map<string, unsigned long long> get_codon_counts() { return codon_counts; }
    map<string, double> get_codon_proportions() { return codon_proportions; }
    pair<string, double> get_min_codon_proportion() { return min_codon_proportion; }
    pair<string, double> get_max_codon_proportion() { return max_codon_proportion; }
    
    map<string, unsigned long long> get_acid_counts() { return acid_counts; }
    map<string, double> get_acid_proportions() { return acid_proportions; }
    pair<string, double> get_min_acid_proportion() { return min_acid_proportion; }
    pair<string, double> get_max_acid_proportion() { return max_acid_proportion; }
    
    string get_sequence() { return sequence; }
    string get_type() { return type; }
    unsigned long long get_total_codons_acids() { return total_codons_acids; }
    
private:
    string type;
    string sequence;
    unsigned long long total_codons_acids;
    
    map<char, unsigned long long> base_counts;
    map<char, double> base_proportions;
    pair<char, double> min_base_proportion;
    pair<char, double> max_base_proportion;
    
    map<string, unsigned long long> codon_counts;
    map<string, double> codon_proportions;
    pair<string, double> min_codon_proportion;
    pair<string, double> max_codon_proportion;
    
    map<string, unsigned long long> acid_counts;
    map<string, double> acid_proportions;
    pair<string, double> min_acid_proportion;
    pair<string, double> max_acid_proportion;
};

map<string, string> create_lookup(); // Function to read txt file data to map codons to corresponding amino acids

#endif /* strand_h */
