//  kir-mapper
//
//  Created by Erick C. Castelli
//  2024 GeMBio.Unesp.
//  erick.castelli@unesp.br

#ifndef functions_hpp
#define functions_hpp

#include <stdio.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <map>
#include <boost/algorithm/string.hpp>
#include <thread>
#include <pwd.h>
#include <mutex>
#include <sys/stat.h>
#include <math.h>

using namespace std;

std::string base_name(std::string const & path);
string split_region (int th, string chr, int start, int end, int stutter);
double median(vector<float> vec);
void screen_message (int size, int left, string message, int enter, int quiet);
bool ends_with(const std::string &filename, const std::string &ext);
string GetStdoutFromCommand(string cmd);
string findfilepath (string v_file);
int phred (char v_char);
string mtrim (string seq, string qual);
double filesize(const char *filename);
string splitcigar (string cigar);
void removefile (string v_file, int v_debug);
bool fileExists(const std::string& filename);
int is_num (string str);
string decompose_mpileup (string cigar);
size_t uiLevenshteinDistance(const std::string &s1, const std::string &s2);
string reverse_and_complement (string seq);
string reverse_string (string seq);
string typing_dna (string gene, string fastq1, string fastq2, int maxselect, int maxerror, string maptype, float mincov);
unsigned long long getTotalSystemMemory();
string splitNcigar (string cigar);
string treat_genotype (string genotype, int copyn);
std::map<std::string, std::string> TypeCDS( string gene, string r1, string r2);
bool isNewer(std::string const &file, std::vector<std::string> const &otherFiles);
string adjust(string gene, string r1, string r2, string out, string nthreads);
string adjust_freeb (string gene, string fq1, string fq2, string out, string threads, string type, string db);
void debug_message (string message);
void fasta_exons (string bedfile, string inputvcf, string referencefile, string gene, string outfile, string master_output);
void check_fasta(string vstart, string vend, string inputvcf, string referencefile, string v_chr, string outfile);
double StandardDeviation(std::vector<double>);
double Variance(std::vector<double>);
double Average(std::vector<double>);
#endif /* functions_hpp */
