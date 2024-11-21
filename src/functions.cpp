//  kir-mapper
//
//  Created by Erick C. Castelli
//  2024 GeMBio.Unesp.
//  erick.castelli@unesp.br
// Contributions from code on tht web


#include "functions.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/range/algorithm/count.hpp>
#include <boost/algorithm/string.hpp>
#include <filesystem>

#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <map>
#include <thread>
#include <pwd.h>
#include <mutex>
#include <sys/stat.h>
#include <math.h>
#include <string>
#include <thread>
#include <unordered_map>
#include <future>
#include <unistd.h>



#include "external.hpp"
#include "ThreadPool.hpp"

namespace fs = std::filesystem;

using namespace std;
mutex mtx;
mutex mtx_type;
mutex mtxadj;

std::string base_name(std::string const & path)
{
    return path.substr(path.find_last_of("/\\") + 1);
}



double median(vector<float> vec)
{
        typedef vector<int>::size_type vec_sz;

        vec_sz size = vec.size();
        if (size == 0)
                throw domain_error("median of an empty vector");

        sort(vec.begin(), vec.end());

        vec_sz mid = size/2;

        return size % 2 == 0 ? (vec[mid] + vec[mid-1]) / 2 : vec[mid];
}


string split_region (int th, string chr, int start, int end, int stutter)
{
    int rsize = end - start;
    float subsize = rsize / th;
    string result = "";
    int substart = start;
    for (int a = 1; a <= th; a++)
    {
        int subend = substart + subsize + stutter;
        if (a == th) {subend = end;}
        result.append("\n" + chr + ":" + to_string(substart) + "-" + to_string(subend));
        substart = substart + subsize;
    }
    return result.substr(1);
}


void screen_message (int size, int left, string message, int enter, int quiet)
{
    if(quiet == 1) {return;}
    if ((message.length()+left) > size) {message = message.substr(0,(size-left-1));}
    cout << "\r";
    for (int a = 0; a < left; a++){cout << " ";}
    cout << message;
    for (int a = 0; a < (((size-left)-message.length())); a++){cout << " ";}
    std::cout.flush();
    if (enter == 1) {cout << endl;}
    return;
}

void debug_message (string message)
{
    if (v_debug == 0) {return;}
    cout << "*** debug: " << message << endl;
}


bool ends_with(const std::string &filename, const std::string &ext)
{
    return ext.length() <= filename.length() &&
    std::equal(ext.rbegin(), ext.rend(), filename.rbegin());
}




string GetStdoutFromCommand(string cmd) {
    
    string data;
    FILE * stream;
    const int max_buffer = 256;
    char buffer[max_buffer];
    cmd.append(" 2>&1");
    
    stream = popen(cmd.c_str(), "r");
    if (stream) {
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
        pclose(stream);
    }
    return data;
}





string findfilepath (string v_file)
{
    string v_filename = v_file.substr(0,(v_file.find_last_of("/"))+1);
    return (v_filename);
}





int phred (char v_char) {
    return int(v_char) - 33;
}






string mtrim (string seq, string qual)
{
    int start = -1;
    int end = -1;
    int low = 0;
    int high = 0;
    
    
    
    for (int $a = 0; $a < qual.length(); $a++) {
        double chr_phred = double(phred(qual.at($a)));
        //        double chr_phred = phred(qual.at($a));
        double P = pow(double(10),(chr_phred / double(-10)));
        
        if ((P <= v_mtrim_error) && (start == -1)) {
            start = $a;
            end = $a;
            continue;
        }
        
        if ((P > v_mtrim_error) && (start >= 0)) {
            end = $a - 1;
            if ((end - start) > (high - low)) {
                low = start;
                high = end;
            }
            start = -1;
            end = -1;
            continue;
        }
    }
    
    if (start >= 0) {
        end = qual.length();
        if ((end - start) > (high - low)) {
            low = start;
            high = end;
        }
    }
    
    if (low >= high)
    {
        seq = "";
        qual = "";
    }
    
    if (low < high)
    {
        seq = seq.substr(low,(high-low+1));
        qual = string(seq.size(),'A');
    }
    
    string v_limits = "";
    if (seq.size() >= v_size) {v_limits = seq + "\n+\n" + qual;}

    return (v_limits);
}







double filesize(const char *filename)
{
    FILE *f = fopen(filename,"rb");  /* open the file in read only */
    
    long size = 0;
    if (fseek(f,0,SEEK_END)==0) /* seek was successful */
        size = ftell(f);
    fclose(f);
    return size;
}




string splitcigar (string cigar)
{
    boost::replace_all(cigar, "S", "S,");
    boost::replace_all(cigar, "H", "H,");
    boost::replace_all(cigar, "M", "M,");
    boost::replace_all(cigar, "D", "D,");
    boost::replace_all(cigar, "I", "I,");
    boost::replace_all(cigar, "N", "N,");
    boost::replace_all(cigar, "P", "P,");
    boost::replace_all(cigar, "X", "X,");
    return cigar;
}

string splitNcigar (string cigar)
{
    boost::replace_all(cigar, "N", "N,");
    return cigar;
}




bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}



bool isNewer(std::string const &file, std::vector<std::string> const &otherFiles) {
    auto newest = std::max_element(otherFiles.begin(), otherFiles.end(),
        [&](std::string const &a, std::string const &b) {
            return fs::last_write_time(a) > fs::last_write_time(b);
        });
    return fs::last_write_time(file) >= fs::last_write_time(*newest);
}



void removefile (string v_file, int v_debug)
{
    if (v_debug == 1) {return;}
    if (fileExists(v_file)) {
        const int result = remove(v_file.c_str());
    }
}




int is_num (string str)
{
    if (str == "0") {return 1;}
    if (str == "1") {return 1;}
    if (str == "2") {return 1;}
    if (str == "3") {return 1;}
    if (str == "4") {return 1;}
    if (str == "5") {return 1;}
    if (str == "6") {return 1;}
    if (str == "7") {return 1;}
    if (str == "8") {return 1;}
    if (str == "9") {return 1;}
    return 0;
}




string decompose_mpileup (string cigar)
{
    cigar.erase(std::remove(cigar.begin(), cigar.end(), '$'), cigar.end());
    cigar.erase(std::remove(cigar.begin(), cigar.end(), '^'), cigar.end());
    cigar.erase(std::remove(cigar.begin(), cigar.end(), ']'), cigar.end());
    cigar.erase(std::remove(cigar.begin(), cigar.end(), '~'), cigar.end());
//    cigar.erase(std::remove(cigar.begin(), cigar.end(), 'N'), cigar.end());
//    cigar.erase(std::remove(cigar.begin(), cigar.end(), 'n'), cigar.end());

    string new_cigar = "";
    for (int a = 0; a < cigar.size(); a++)
    {
        string sub = cigar.substr(a,1);
        
        string next;
        if (a < cigar.size()) {next = cigar.substr(a+1,1);}
        
        if ((next == "+") || (next == "-"))
        {
            int start_search = a+2;
            int end_search = a+2;
            for (end_search = start_search; end_search < cigar.size();end_search++) {if (is_num(cigar.substr(end_search,1)) == 0) {break;}}

            int get = end_search - start_search;
            int size_n = stoi(cigar.substr(start_search, get));
            string size_te = cigar.substr(start_search, get);
            

            string value = "";
            value.append(cigar.substr(a,1));
            value.append(next);
            value.append(size_te);
            a = a + 2 + size_te.size();
            value.append(cigar.substr(a,size_n));
            a = a + size_n - 1;
            new_cigar.append(value + ";");
            continue;
        }
        
        if (sub == ".") {new_cigar.append(sub + ";");continue;}
        if (sub == ",") {new_cigar.append(sub + ";");continue;}
        if (sub == "A") {new_cigar.append(sub + ";");continue;}
        if (sub == "T") {new_cigar.append(sub + ";");continue;}
        if (sub == "C") {new_cigar.append(sub + ";");continue;}
        if (sub == "G") {new_cigar.append(sub + ";");continue;}
        if (sub == "a") {new_cigar.append(sub + ";");continue;}
        if (sub == "t") {new_cigar.append(sub + ";");continue;}
        if (sub == "c") {new_cigar.append(sub + ";");continue;}
        if (sub == "g") {new_cigar.append(sub + ";");continue;}
        if (sub == "N") {new_cigar.append(sub + ";");continue;}
        if (sub == "n") {new_cigar.append(sub + ";");continue;}
        if (sub == "*") {new_cigar.append(sub + ";");continue;}
    }
    std::transform(new_cigar.begin(), new_cigar.end(),new_cigar.begin(), ::toupper);
    return new_cigar;
}




size_t uiLevenshteinDistance(const std::string &s1, const std::string &s2)
{
    // from https://rosettacode.org/wiki/Levenshtein_distance#C.2B.2B
    
    const size_t m(s1.size());
    const size_t n(s2.size());
 
    if( m==0 ) return n;
    if( n==0 ) return m;
 
    size_t *costs = new size_t[n + 1];
 
    for( size_t k=0; k<=n; k++ ) costs[k] = k;
 
    size_t i = 0;
    for ( std::string::const_iterator it1 = s1.begin(); it1 != s1.end(); ++it1, ++i )
    {
        costs[0] = i+1;
        size_t corner = i;
 
        size_t j = 0;
        for ( std::string::const_iterator it2 = s2.begin(); it2 != s2.end(); ++it2, ++j )
        {
            size_t upper = costs[j+1];
            if( *it1 == *it2 )
            {
          costs[j+1] = corner;
      }
            else
      {
        size_t t(upper<corner?upper:corner);
                costs[j+1] = (costs[j]<t?costs[j]:t)+1;
      }
 
            corner = upper;
        }
    }
 
    size_t result = costs[n];
    delete [] costs;
 
    return result;
}




string reverse_and_complement (string seq)
{
    string rev = "";
    for (int a = seq.size(); a >= 0;a--)
    {
        string chr = seq.substr(a,1);
        if (chr == "T") {rev.append("A");continue;}
        if (chr == "A") {rev.append("T");continue;}
        if (chr == "G") {rev.append("C");continue;}
        if (chr == "C") {rev.append("G");continue;}
        if (chr == "N") {rev.append("N");continue;}
    }
    return rev;
}

string reverse_string (string seq)
{
    string rev = "";
    for (int a = seq.size(); a >= 0;a--)
    {
        string chr = seq.substr(a,1);
    }
    return rev;
}







string typing_dna (string gene, string fastq1, string fastq2, int maxselect, int maxerror, string maptype, float mincov)
{

        if (filesize(fastq1.c_str()) < 25000) {return "fail:low cov low size";}

    
        int motif_size = 20;
        unordered_map <string,int> motifs;

        ifstream list;
        string reffile = v_db + "/typing/dna/" + gene + "/reference/" + gene + ".fa";
        boost::replace_all(reffile, "\\ ", " ");

        if (! fileExists(reffile.c_str())) {gene_type[gene] = "."; return "disabled";}
        string vlog = v_output + "/log/typing_select.log";
        string vout = v_output + "/" + v_sample + gene + ".typing_search.sam";
        string command = "";
        if (maptype == "paired") {command = v_bwa + " mem -a -t " + v_threads + " '" + reffile + "' " + fastq1 + " " + fastq2 + " > " + vout + " 2>" + vlog;}
        if (maptype == "single") {command = v_bwa + " mem -a -t " + v_threads + " '" + reffile + "' " + fastq1 + " > " + vout + " 2>" + vlog;}
        system(command.c_str());

    
        map <string,float> reference;
        unordered_map <string,float> counterhits;
        unordered_map <string,float> error;
        vector <int> readsizes;
        
        ifstream samsearch;
        samsearch.open (vout.c_str());
        for( std::string item; getline( samsearch, item ); )
        {
                if (item == "") {continue;}
                if (item.substr(0,1) == "[") {continue;}

                if (item.substr(0,3) == "@SQ") {
                                vector <string> samdata;
                                boost::split(samdata,item,boost::is_any_of("\t"));
                                string target = samdata[1].substr(3);
                                string size = samdata[2].substr(3);
                                reference[target] = stoi(size) - 120;
                                continue;
                }
                if (item.substr(0,1) == "@") {continue;}
                
                vector <string> samdata;
                boost::split(samdata,item,boost::is_any_of("\t"));

                if (samdata[2] == "*") {continue;}
                if (samdata[9] != "*") {
                        readsizes.push_back(samdata[9].size());
                        for (int a = 0; a < (samdata[9].size() - motif_size); a++)
                        {
                                string sub = samdata[9].substr(a,motif_size);
                                motifs[sub]++;
                        }
                }
                
                if (samdata.size() > 10) {
                        if (samdata[11].substr(0,2) == "NM") {
                                int nm = stoi(samdata[11].substr(5));
                                if (nm > maxerror) {continue;}
                                counterhits[samdata[2]]++;
                                error[samdata[2]] = error[samdata[2]] + nm;
                        }
                }
        }
        samsearch.close();
        removefile(vout,v_debug);
        

        int sumsizes = 0;
        for (auto & item : readsizes) {sumsizes = sumsizes + item;}
        float meanreadsize = float(sumsizes) / float(readsizes.size());
        
        float highercov = 0;
        unordered_map <string,int> banned;
        for (auto & item : reference)
        {
                string allele = item.first;
                float allelecov = (float(counterhits[allele]) * meanreadsize) / float(item.second);
                if (highercov < allelecov) {highercov = allelecov;}
        }
        
        
        if (highercov < mincov) {return "fail:low cov after BWA";}
        
        for (auto & item : reference)
        {
                string allele = item.first;
                float allelecov = (float(counterhits[allele]) * meanreadsize) / float(item.second);
                if (allelecov < (0.7 * highercov)) {banned[allele] = 1; continue;}
        }


        string motif_file = v_db + "/typing/dna/" + gene + "/motifs/" + gene + ".txt";
        boost::replace_all(motif_file, "\\ ", " ");

        ifstream motif;
        motif.open(motif_file.c_str());
        unordered_map <string,int> accepted;
        for( std::string line; getline( motif, line ); )
        {
                vector <string> item;
                boost::split(item,line,boost::is_any_of("\t"));
                string allele = item[0];
                
                if (banned.find(allele) != banned.end()) {continue;}
                
                for (int a = 1; a < item.size(); a++)
                {
                        string sub = item[a];
                        if (motifs.find(sub) != motifs.end()) {
                                if (motifs[sub] >= (0.2 * highercov)) {
                                        accepted[allele] = 1;
                                        break;
                                }
                        }
                }
        }
        motif.close();
        
        


        ThreadPool select_pool(stoi(v_threads));
        std::vector< std::future<int> > select_step_results;
        map <float,string> scores;
        for (auto & target : counterhits)
        {
                string current = target.first;
                select_step_results.emplace_back(
                                        select_pool.enqueue([&scores,accepted,current,&counterhits,&error,&reference]
                {
                        if (accepted.find(current) == accepted.end()) {return 1;}
                        float cov = (counterhits[current]*100) / reference[current];
                        float noise = error[current] / counterhits[current];
                        float score = cov / noise;
                        mtx.lock();
                        scores[score].append("," + current);
                        mtx.unlock();
                        return 1;
                })
                );
        }
        for(auto && result: select_step_results){result.get();}
        

        
        int count = 0;
        vector <string> selected;
        int best_score_selection = 0;
        for (auto it = scores.rbegin(); it != scores.rend(); it++) {
                                if (best_score_selection == 0) {best_score_selection = it->first;}
                                vector <string> item;
                                boost::split(item,it->second,boost::is_any_of(","));
                                for (auto & allele : item)
                                {
                                        if (allele == "") {continue;}
                                        selected.push_back(allele);
                                        count++;
                                }
                                if (count > maxselect) {break;}
        }

        if (v_debug == 1) {
                for (auto & allele : selected)
                {
                        cout << endl << gene << " allele -> " << allele;
                }
        }
        
        if (selected.size() == 0) {return "fail: no preselect";}



        unordered_map <string,string> head;
        unordered_map <string,int> reads;
        map <pair<string,string>,string> mapdata;
        map <pair<string,string>,int> nmdata;
        
        ThreadPool loading_pool(stoi(v_threads));
        std::vector< std::future<int> > loading_results;
        
        for (auto & allele : selected)
        {
                string current = allele;
                loading_results.emplace_back(
                                        loading_pool.enqueue([current,&mapdata,&nmdata,&head,gene,fastq1,fastq2,&reads,maxerror,maptype]
                {
                        string file = current;
                        replace( file.begin(), file.end(), ':', '_');
                        replace( file.begin(), file.end(), '*', '_');
                        
                        string ref = "'" + v_db + "/typing/dna/" + gene + "/alleles/" + file + ".fa'";
                        boost::replace_all(ref, "\\ ", " ");

                        string cmd = "";
                        if (maptype == "paired"){cmd = v_bwa + " mem -v 1 -a " + ref + " " + fastq1 + " " + fastq2;}
                        if (maptype == "single"){cmd = v_bwa + " mem -v 1 -a " + ref + " " + fastq1;}
                        string out = GetStdoutFromCommand(cmd);
                        
                        vector <string> lines;
                        boost::split(lines,out,boost::is_any_of("\n"));
                        for (auto & item: lines)
                        {
                                if (item == "") {continue;}
                                if (item.substr(0,1) == "[") {continue;}
                                if (item.substr(1,2) == "SQ")
                                {
                                        vector <string> data;
                                        boost::split(data,item,boost::is_any_of("\t"));
                                        string allele = data[1].substr(3);
                                        mtx.lock();
                                        head[allele] = item;
                                        mtx.unlock();
                                        continue;
                                }
                                if (item.substr(0,1) == "@") {continue;}
                                vector <string> data;
                                boost::split(data,item,boost::is_any_of("\t"));
                                if (data[2] == "*") {continue;}
                                if (data.size() >= 11)
                                {
                                        int nm = stoi(data[11].substr(5));
                                        pair <string,string> key = make_pair(data[2],data[0]);
                                        mtx.lock();
                                        mapdata[key].append("\n" + item);
                                        nmdata[key] = nmdata[key] + nm;
                                        reads[data[0]] = 1;
                                        mtx.unlock();
                                }
                        }
                        return 1;
                })
                );
        }
        for(auto && result: loading_results){result.get();}
        if (reads.size() == 0) {return "fail: no loaded data";}




        unordered_map <string,int> used;
        int counter = 0;
        
        string best_hit = "";
        float best_score = 0;
        
        ThreadPool pool(stoi(v_threads));
        std::vector< std::future<int> > results;

        for (auto & alleleA : selected)
        {
                string currentA = alleleA;
                
                for (auto & alleleB : selected)
                {
                        string currentB = alleleB;
                        if (used.find(currentB) != used.end()) {continue;}
                        results.emplace_back(
                                                pool.enqueue([counter, gene, &head, currentA, currentB, &reads, &nmdata, &mapdata, &reference, maxerror,&best_score, &best_hit]
                        {
                                        string samA = "";
                                        string samB = "";
                                        unordered_map <string,int> readsA;
                                        unordered_map <string,int> readsB;
                                        unordered_map <string,int> reads_general;

                                        samA.append(head[currentA]);
                                        samB.append(head[currentB]);
                                
                                        for (auto & readid : reads)
                                        {
                                                int nmA = 1000;
                                                int nmB = 1000;
                                                pair <string,string> keyA = make_pair(currentA,readid.first);
                                                pair <string,string> keyB = make_pair(currentB,readid.first);
                                                if (nmdata.find(keyA) != nmdata.end()) {nmA = nmdata[keyA];}
                                                if (nmdata.find(keyB) != nmdata.end()) {nmB = nmdata[keyB];}
                                                if ((nmA == 1000) && (nmB == 1000)) {continue;}
                                                if ((nmA > 1) && (nmB == 1000)){continue;}
                                                if ((nmB > 1) && (nmA == 1000)){continue;}
                                                if ((nmA > maxerror) && (nmB > maxerror)){continue;}

                                                if (nmA <= nmB) {samA.append(mapdata[keyA]);readsA[readid.first] = 1;}
                                                if (nmA >= nmB) {samB.append(mapdata[keyB]);readsB[readid.first] = 1;}
                                                reads_general[readid.first] = 1;
                                        }
                                        

                                        
                                        string fileA = currentA;
                                        replace( fileA.begin(), fileA.end(), ':', '_');
                                        replace( fileA.begin(), fileA.end(), '*', '_');
                                        string fileB = currentB;
                                        replace( fileB.begin(), fileB.end(), ':', '_');
                                        replace( fileB.begin(), fileB.end(), '*', '_');
                                     
                                        string ref = "'" + v_db + "/typing/dna/" + gene + "/alleles/" + fileA + ".fa'";
                                        boost::replace_all(ref, "\\ ", " ");

                                        
                                        string cmd = "";
                                        string pileA = "";
                                        string outfile = v_output + "/" + v_sample + "_pair_" + fileA + "_" + fileB + ".1.sam";
                                        mtx.lock();
                                        ofstream out;
                                        out.open (outfile.c_str());
                                        out << samA << endl;
                                        out.close();
                                        mtx.unlock();
                                        cmd = v_samtools + " sort " + outfile + " | " + v_samtools + " mpileup --verbosity 0 -a --reference " + ref + " -";
                                        pileA = GetStdoutFromCommand(cmd);
                                        mtx.lock(); removefile(outfile,0); mtx.unlock();
                                        
                                        string pileB = "";
                                        if (currentA != currentB) {
                                                ref = "'" + v_db + "/typing/dna/" + gene + "/alleles/" + fileB + ".fa'";
                                                boost::replace_all(ref, "\\ ", " ");

                                                cmd = "";
                                                outfile = v_output + "/" + v_sample + "_pair_" + fileA + "_" + fileB + ".2.sam";
                                                mtx.lock();
                                                out.open (outfile.c_str());
                                                out << samB << endl;
                                                out.close();
                                                mtx.unlock();
                                                cmd = v_samtools + " sort " + outfile + " | " + v_samtools + " mpileup --verbosity 0 -a --reference " + ref + " -";
                                                pileB = GetStdoutFromCommand(cmd);
                                                mtx.lock(); removefile(outfile,0); mtx.unlock();
                                        }
                                
                                        vector <string> pileAdata;
                                        vector <string> pileBdata;
                                        boost::split(pileAdata,pileA,boost::is_any_of("\n"));
                                        boost::split(pileBdata,pileB,boost::is_any_of("\n"));
                                        
                                        float covA_sum;
                                        float covA_count;
                                        float errorA_sum = 1;
                                        float errorA_count = 1;
                                        
                                        for (auto item : pileAdata)
                                        {
                                                if (item == "") {continue;}
                                                if (item.substr(0,1) == "[") {continue;}
                                                vector <string> data;
                                                boost::split(data,item,boost::is_any_of("\t"));
                                                if (((data[2] == "N") || (data[2] == "0")) || (data[2] == "*")) {continue;}
                                                float cov = stof(data[3]);
                                                if (cov == 0) {continue;}
                                                string cigar = data[4];
                                                float ref = 0;
                                                for (int a = 0; a < cigar.size(); a++)
                                                {
                                                    if ((cigar.substr(a,1) == ".") || (cigar.substr(a,1) == ",")) {ref++;}
                                                }
                                                
                                                float indel = 0;
                                                for (int a = 0; a < cigar.size(); a++)
                                                {
                                                        if ((cigar.substr(a,1) == "+") || (cigar.substr(a,1) == "-")) {indel++;}
                                                }
                                                errorA_sum = errorA_sum + ((cov - ref) + indel);
                                                if (((cov - ref) + indel) > 0) {errorA_count++;}
                                        }
                                        

                                        float covB_sum;
                                        float covB_count;
                                        float errorB_sum = 1;
                                        float errorB_count = 1;
                                        
                                        if (currentA != currentB) {
                                            for (auto item : pileBdata)
                                            {
                                                    if (item == "") {continue;}
                                                    if (item.substr(0,1) == "[") {continue;}
                                                    vector <string> data;
                                                    boost::split(data,item,boost::is_any_of("\t"));
                                                    if (((data[2] == "N") || (data[2] == "0")) || (data[2] == "*")) {continue;}
                                                    float cov = stof(data[3]);
                                                    if (cov == 0) {continue;}
                                                    string cigar = data[4];
                                                    float ref = 0;
                                                    for (int a = 0; a < cigar.size(); a++)
                                                    {
                                                        if ((cigar.substr(a,1) == ".") || (cigar.substr(a,1) == ",")) {ref++;}
                                                    }
                                                    
                                                    float indel = 0;
                                                    for (int a = 0; a < cigar.size(); a++)
                                                    {
                                                            if ((cigar.substr(a,1) == "+") || (cigar.substr(a,1) == "-")) {indel++;}
                                                    }
                                                    errorB_sum = errorB_sum + ((cov - ref) + indel);
                                                    if (((cov - ref) + indel) > 0) {errorB_count++;}
                                            }
                                        }
                                        
                                        
                                        if (currentA == currentB) {
                                            errorB_sum = errorA_sum;
                                            errorB_count = errorA_count;
                                        }
                                        
                                        float errorA = errorA_sum / errorA_count;
                                        float errorB = errorB_sum / errorB_count;
                                        float errorsum = errorA + errorB;

                                        float size = (max(reference[currentA], reference[currentB])) / 1000;
                                        float sizekb = ((float)((int)(size * 100))) / 100;
                                        float abund = float(reads_general.size()) / sizekb;
                                        
                                        float score = abund / errorsum;
                                            
                                        
                                        mtx.lock();
                                        if (v_debug == 1) {
                                            cout << currentA << "\t" << currentB << "\t" << errorA << "\t" << errorB << "\t" << errorsum << "\t" << abund << "\t" << score << endl;
                                        }
                                        if (best_score == score) {best_hit.append("," + currentA + "\t" + currentB);}
                                        if (best_score < score) {best_score = score; best_hit = currentA + "\t" + currentB;}
                                        mtx.unlock();
                                        return counter;
                        })
                        );
                                counter++;
                                used[currentA] = 1;
                }
        }
        for(auto && result: results){result.get();}

        if (best_hit == "") {return "fail: no best hit";}
        
        vector <string> typed;
         if (best_hit != "") {
            vector <string> subbest;
            boost::split(subbest,best_hit,boost::is_any_of(","));
            best_hit = subbest[0];
            boost::split(typed,best_hit,boost::is_any_of("\t"));
         }
        
        
        if (best_hit != "")
        {
            string fq;
            if (v_map_type == "single") {fq = " " + v_output + v_sample + gene + "_R0.fastq";}
            if (v_map_type == "paired") {fq = " " + v_output + v_sample + gene + "_R1.fastq " + v_output + v_sample + gene + "_R2.fastq";}
            
            vector <string> typed;
            boost::split(typed,best_hit,boost::is_any_of("\t"));
            string file = typed[0];
            replace( file.begin(), file.end(), ':', '_');
            replace( file.begin(), file.end(), '*', '_');
            string ref = " '" + v_db + "/typing/dna/" + gene + "/alleles/" + file + ".fa' ";
            boost::replace_all(ref, "\\ ", " ");

            v_command = v_bwa + " mem -v 1 -t " + v_threads + ref + fq;
            string out = GetStdoutFromCommand(v_command);
            
            
            unordered_map <string,int> nm_allele_A;
            unordered_map <string,int> nm_allele_B;
            
            vector <string> samdata;
            boost::split(samdata,out,boost::is_any_of("\n"));
            for (auto & item : samdata)
            {
                if (item == "") {continue;}
                if (item.substr(0,1) == "[") {continue;}
                vector <string> data;
                boost::split(data,item,boost::is_any_of("\t"));
                    
                if (data[2] == "*") {continue;}
                if (data.size() >= 11)
                {
                    int nm = stoi(data[11].substr(5));
                    nm_allele_A[data[0]] = nm_allele_A[data[0]] + nm;
                }
                

                
            }
            for (auto &item : nm_allele_A)
            {
                if (item.second <= 1)
                {
                    pair <string,string> key = make_pair(gene,item.first);
                    if (reads_possibly_mapped.find(key) != reads_possibly_mapped.end())
                    {
                        reads_possibly_mapped[key] = item.second;
                    }
                    else
                    {
                        if (item.second < reads_possibly_mapped[key]) {
                            reads_possibly_mapped[key] = item.second;
                        }
                    }
                }
            }
 
            if (typed[0] != typed[1]) {
                file = typed[1];
                replace( file.begin(), file.end(), ':', '_');
                replace( file.begin(), file.end(), '*', '_');
                ref = " '" + v_db + "/typing/dna/" + gene + "/alleles/" + file + ".fa' ";
                boost::replace_all(ref, "\\ ", " ");

                v_command = v_bwa + " mem -v 1 -t " + v_threads + ref + fq;
                    
                out = GetStdoutFromCommand(v_command);
                
                samdata.clear();
                boost::split(samdata,out,boost::is_any_of("\n"));
                for (auto & item : samdata)
                {
                    if (item == "") {continue;}
                    if (item.substr(0,1) == "[") {continue;}
                            
                    vector <string> data;
                    boost::split(data,item,boost::is_any_of("\t"));
                    
                            
                    if (data[2] == "*") {continue;}
                    if (data.size() >= 11)
                    {
                        int nm = stoi(data[11].substr(5));
                        nm_allele_B[data[0]] = nm_allele_B[data[0]] + nm;
                    }
                }
                    
                for (auto &item : nm_allele_B)
                {
                    if (item.second <= 1 )
                    {
                        pair <string,string> key = make_pair(gene,item.first);
                        if (reads_possibly_mapped.find(key) != reads_possibly_mapped.end())
                        {
                            if (reads_possibly_mapped[key] > item.second)
                            {
                                reads_possibly_mapped[key] = item.second;
                            }
                        }
                        else {reads_possibly_mapped[key] = item.second;}
                    }
                }
            }
        }
        

        if (best_hit != "") {return best_hit;}
        return "fail: no best hit";
}






unsigned long long getTotalSystemMemory()
{
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}


string treat_genotype (string genotype, int copyn)
{
    
	
	
	vector <string> fields;
    boost::split(fields,genotype,boost::is_any_of(":"));
    if (fields[0].find(".") != std::string::npos) {return genotype;}
    vector <string> alleles;
    
    if (copyn != 1) {
        boost::split(alleles,fields[0],boost::is_any_of("/"));
    }
    if (copyn == 1) {alleles.push_back(fields[0]);}
    
    float cov = stof(fields[1]);
    string depth = fields[2];
    vector <string> subdepth;
    boost::split(subdepth,depth,boost::is_any_of(","));
    
    
    int change = 0;
    for (int a = 0; a < alleles.size();a++)
    {
        float alleledepth = stof(subdepth[stoi(alleles[a])]);
        if (alleledepth <= 2) {alleles[a] = ".";change++;}
    }
    
    if (change != 0)
    {
        string newgeno = boost::algorithm::join(alleles, "/");
        return newgeno + ":" + fields[1] + ":" + fields[2] + ":" + fields[3];
    }

    
    map <string,int> allelecount;
    for(auto item : alleles) {allelecount[item]++;}
    string type = "";
    if (allelecount.size() == 1) {type = "homoz";}
    if (allelecount.size() > 1) {type = "heteroz";}

    

    if (type == "homoz")
    {
        if (cov < 6)
        {
            for (int a = 1; a < alleles.size();a++)
            {
                alleles[a] = ".";
            }
        }
        

        if (cov >= 6)
        {
            for (int a = 0; a < alleles.size();a++)
            {
                float alleledepth = stof(subdepth[stoi(alleles[a])]);
                float ratio = alleledepth / cov;
                if (copyn == 1){
                    if (ratio <= 0.7) {alleles[a] = ".";break;}
                }
                if (copyn > 1){
                    if (ratio <= 0.85) {alleles[a] = ".";break;}
                }
            }
        }

        string newgeno = boost::algorithm::join(alleles, "/");
        return newgeno + ":" + fields[1] + ":" + fields[2] + ":" + fields[3];
    }
    
    if (type == "heteroz")
    {
        // is there a major allele?
        string majorallele = "";
        
        for (auto item : alleles)
        {
            float alleledepth = stof(subdepth[stoi(item)]);
            float ratio = alleledepth / cov;
            if (ratio >= 0.9) {majorallele = item;}
        }
        
        if (majorallele != "")
        {
            for (int a = 0; a < alleles.size();a++)
            {
                alleles[a] = majorallele;
            }
            string newgeno = boost::algorithm::join(alleles, "/");
            return newgeno + ":" + fields[1] + ":" + fields[2] + ":" + fields[3];
        }
        
        if (majorallele == "")
        {
            for (int a = 0; a < alleles.size();a++)
            {
                float alleledepth = stof(subdepth[stoi(alleles[a])]);
                float ratio = alleledepth / cov;
                if (ratio <= 0.25) {alleles[a] = ".";}
            }
            string newgeno = boost::algorithm::join(alleles, "/");
            return newgeno + ":" + fields[1] + ":" + fields[2] + ":" + fields[3];
        }
    }
    return genotype;
}



string adjust_freeb (string gene, string fq1, string fq2, string out, string threads, string type, string db)
{
    int size = 20;
    int debug = v_debug;
    string fas = db + "/adjustment/dna/fasta/" + type + "/" + gene;
    int minnm = 10;
    int nselect = 100;
    
    string fails = "failed";
    
    vector <string> files_to_remove;
    
    int match_level = 6;
    if (type == "cds") {match_level = 10;}
    
    
    if (gene == "") {return fails;}
    if (fq1 == "") {return fails;}
    if (! fileExists(fq1)) {return fails;}
    if ((fq2 != "") && (! fileExists(fq2))) {return fails;}
    if (out == "") {return fails;}
    if (! fileExists(out)) {return fails;}
    if (threads == "") {threads = "1";}
    if (type == "") {type = "cds";}
    if (fas == "") {return fails;}
    if (! fileExists(fas)) {return fails;}
    

    
    map <string,int> banned;
    banned["KIR2DS2.016"];
    banned["KIR2DL1.026"];
    banned["KIR2DL1.026N"];

/*

    debug_message("Selecting alleles ... start");
    map <string,int> motif;
    
    ifstream FQ;
    FQ.open (fq1.c_str());
    for( std::string item; getline( FQ, item ); )
    {
        string id = item;
        getline( FQ, item );
        string seq = item;
        getline( FQ, item );
        string info = item;
        getline( FQ, item );
        string qual = item;
        
        string rev = reverse_and_complement(seq);
        if (seq.size() < size) {continue;}
        for (int a = 0; a < (seq.size() - size); a++)
        {
            string sub = seq.substr(a,size);
            motif[sub]++;
        }
        for (int a = 0; a < (rev.size() - size); a++)
        {
            string sub = rev.substr(a,size);
            motif[sub]++;
        }
    }
    FQ.close();
    
    if (fq2 != "") {
        ifstream FQ;
        FQ.open (fq2.c_str());
        for( std::string item; getline( FQ, item ); )
        {
            string id = item;
            getline( FQ, item );
            string seq = item;
            getline( FQ, item );
            string info = item;
            getline( FQ, item );
            string qual = item;
            
            string rev = reverse_and_complement(seq);
            for (int a = 0; a < (seq.size() - size); a++)
            {
                string sub = seq.substr(a,size);
                motif[sub]++;
            }
            for (int a = 0; a < (rev.size() - size); a++)
            {
                string sub = rev.substr(a,size);
                motif[sub]++;
            }
        }
        FQ.close();
    }
   */ 
    
    
    vector <string> fasfiles;
    for (const auto & entry : fs::directory_iterator(fas))
    {
        if (entry.path().extension() != ".fas") {continue;}
        fasfiles.push_back(entry.path());
    }
    if (fasfiles.size() == 0) {return fails;}
    
    
    vector <string> selected_alleles_list;
    for (auto item : fasfiles)
    {
        string file = item;
        string allele = base_name(file);
        boost::replace_all(allele, ".fas", "");
        selected_alleles_list.push_back(allele);   
    }
  /*  
    map <float,string> allele_motif;
    for (auto item : fasfiles)
        {
            string file = item;
            string allele = base_name(file);
            
            boost::replace_all(allele, ".fas", "");
            
            float tested = 0;
            float valid = 0;
            
            ifstream FAS;
            FAS.open (file.c_str());
            for( std::string item; getline( FAS, item ); )
            {
                string id = item;
                getline( FAS, item );
                string seq = item;
                for (int a = 0; a < (seq.size() - size); a++)
                {
                    string sub = seq.substr(a,size);
                    if (sub.find("N") != std::string::npos) {continue;}
                    tested++;
                    if (motif.find(sub) != motif.end()) {
                        if (motif[sub] >= match_level) {valid++;}
                    }
                }
            }
            FAS.close();
            float ratio = valid / tested;
            allele_motif[ratio].append("," + allele);
        }
    
    
    vector <string> selected_alleles_list;
    for (auto iter = allele_motif.rbegin(); iter != allele_motif.rend(); ++iter) {
        
        vector <string> alleles;
        string all = iter->second.substr(1);
        float score = iter->first;
        boost::split(alleles,all,boost::is_any_of(","));
        
        for (auto item: alleles)
        {
            if (banned.find(item) != banned.end()) {continue;}
            selected_alleles_list.push_back(item);
        }
        //if (selected_alleles_list.size() >= nselect) {break;}
    }
    
    if (selected_alleles_list.size() == 0) {return fails;}
    */

    debug_message("List of alleles ...");


    
    map <string,int> select_alleles;
    for (auto item : selected_alleles_list)
    {
        select_alleles[item] = 1;
        debug_message(item);
        
    }
    
    string file = out + "/list_of_allele.preselection.txt";
    ofstream alleleout (file.c_str());
    selected_alleles_list.clear();
    for (auto item : select_alleles)
    {
        selected_alleles_list.push_back(item.first);
        alleleout << item.first << endl;
    }
    
    

    




    
    
    
    
    
    
    
    
    
    debug_message("Quick genotyping ... start");
    
    vector <string> bamlist;
    for (auto item : selected_alleles_list)
    {
        string bam = db + "/adjustment/dna/bam/" + gene + "/" + item + ".bam";
        if (fileExists(bam)) {bamlist.push_back(bam);}
    }
    
    string v_reference = "'" + db + "/reference/dna/loci/" + gene + ".fas'";
    boost::replace_all(v_reference, "\\ ", " ");
    string samtmp = out + "/" + gene + ".tmp.sam";
    string sam = out + "/" + gene + ".sam";
    string bam = out + "/" + gene + ".bam";
    string cmd = v_bwa + " mem -t " + threads + " -B 2 -O 3,3 -L 3,3 -o " + samtmp + " -R '@RG\\tID:Test\\tLB:Test\\tSM:Test\\tPL:illumina\\tPU:Test' " + v_reference + " " + fq1;
    if (fq2 != "") {cmd.append(" " + fq2);}
    GetStdoutFromCommand(cmd);

    files_to_remove.push_back(sam);
    files_to_remove.push_back(samtmp);

    
        
    ifstream samdata (samtmp.c_str());
    ofstream samout (sam.c_str());
    int v_pos_correct = position_hg38[gene] - 1;
    samout << "@SQ\tSN:chr19\tLN:58617616" << endl;
    samout << "@SQ\tSN:chr19_KI270890v1_alt\tLN:184499" << endl;
    samout << "@SQ\tSN:chr19_KI270921v1_alt\tLN:282224" << endl;
    samout << "@SQ\tSN:chr19_KI270938v1_alt\tLN:1066800" << endl;
    samout << "@SQ\tSN:chr6\tLN:170805979" << endl;
    
    for( std::string line; getline( samdata, line ); )
    {
        if (line == "") {continue;}
        if (line.substr(0,3) == "@RG") {samout << line << endl; continue;}
        if (line.substr(0,1) == "@") {continue;}
        vector<string> sam_line;
        boost::split(sam_line,line,boost::is_any_of("\t"));
        
        sam_line[3] = to_string((stoi(sam_line[3]) + v_pos_correct));
        sam_line[7] = to_string((stoi(sam_line[7]) + v_pos_correct));
        sam_line[2] = chr_hg38[gene];
        
        string newline = "";
        for (auto item : sam_line) {newline.append("\t" + item);}
        samout << newline.substr(1) << endl;
    }
    samdata.close();
    samout.close();
    
    
            
    cmd = v_samtools + " sort -@ " + threads + " -o " + bam + " " + sam;
    GetStdoutFromCommand(cmd);
    cmd = v_samtools + " index " + bam;
    GetStdoutFromCommand(cmd);


    
    string cnref = out + "/" + gene + ".cn.txt";
    files_to_remove.push_back(cnref);
    ofstream cn;
    cn.open (cnref.c_str());
    cn << "Test\t2" << endl;
    for (auto item : selected_alleles_list)
    {
        cn << "ref." << item << "\t1" << endl;
    }
    cn.close();
    
    
    
    
    vector <string> subregions;

    string bed = "";
    if (type == "cds"){bed = db + "/genotype/bed/" + gene + ".exons.bed"; }
    if (type == "full"){bed = db + "/genotype/bed/" + gene + ".full.bed"; }

    if (type == "full") {
        ifstream exon(bed.c_str());
        for( std::string line; getline( exon, line ); )
        {
            if (line == "") {continue;}
            vector <string> data;
            boost::split(data,line,boost::is_any_of("\t"));
            int start = stoi(data[1]) - 5;
            int end = stoi(data[2]) + 5;
            string sub = data[0] + ":" + to_string(start) + "-" + to_string(end);
            string subregions_txt = split_region(stoi(v_threads),data[0],start,end,20);
            boost::split(subregions,subregions_txt,boost::is_any_of("\n"));
    //        subregions.push_back(sub);
        }
        exon.close();
    }
    if (type == "cds") {
        ifstream exon(bed.c_str());
        for( std::string line; getline( exon, line ); )
        {
            if (line == "") {continue;}
            vector <string> data;
            boost::split(data,line,boost::is_any_of("\t"));
            int start = stoi(data[1]) - 5;
            int end = stoi(data[2]) + 5;
            string sub = data[0] + ":" + to_string(start) + "-" + to_string(end);
            subregions.push_back(sub);
        }
        exon.close();
    }
    
    ThreadPool poolfree(stoi(v_threads));
    std::vector< std::future<int> > results_free;
    
    int call_count = 0;
    for (auto region : subregions)
    {
        call_count++;
        results_free.emplace_back(
            poolfree.enqueue([call_count,region,gene,out,&files_to_remove,db,bam,cnref,&bamlist,type]
                {
                    string current_region = "";
                    current_region = region;
                    boost::replace_all(current_region, "chrchr", "chr");
                    
                    string region_out = current_region;
                    boost::replace_all(region_out, ":", ".");
                    boost::replace_all(region_out, "-", ".");

                    string vcfout = out + "/" + gene + "." + region_out + ".vcf";
                    mtxadj.lock(); files_to_remove.push_back(vcfout); mtxadj.unlock();
                    string hg38ref = db + "/reference/hg38/reference.fasta";

                    string cmd = v_freebayes + " -A " + cnref + " -f " + hg38ref + " -b " + bam + " -r " + region;
                    for (auto item : bamlist) {cmd.append(" -b " + item);}
                    cmd.append(" -v " + vcfout);
                    GetStdoutFromCommand(cmd);
                    
                    if ((! fileExists(vcfout)) || (filesize(vcfout.c_str()) == 0)) {
                        string out = GetStdoutFromCommand(cmd.c_str());
                        int failed = 0;
                        if (filesize(vcfout.c_str()) == 0) {failed = 1;}
                        if (out.find("Failed") != std::string::npos) {failed = 1;}
                        if (out.find("Too many") != std::string::npos) {failed = 2;}
                        if (failed == 1) {
                            mtxadj.lock();
                            cout << endl << "Freebayes failed. Quitting..." << endl;
                            mtxadj.unlock();
                            exit(0);
                            return 1;
                        }
                        if (failed == 2) {
                            mtxadj.lock();
                            cout << endl << "Freebayes failed because of too many files." << endl;
                            cout << "You should adjust 'ulimit -n' to a higher value. Quitting..." << endl;
                            mtxadj.unlock();
                            exit(0);
                            return 1;
                        }
                    }
                    return 1;
                })
        );

    }
    
    for(auto && result: results_free){result.get();} // waiting for all threads
    debug_message("Quick genotyping ... freebayes ...end");


    debug_message("Loading vcf ... start");

      int count = 0;
      vector <string> head;
      vector <string> samples;
      map <int,string> snp_data;
      map <pair<int,string>,string> vcf_data;
      string sampleline = "";

      for (auto region : subregions)
      {
          string current_region = "";
          current_region = region;
          boost::replace_all(current_region, "chrchr", "chr");
          
          string region_out = current_region;
          boost::replace_all(region_out, ":", ".");
          boost::replace_all(region_out, "-", ".");

          string vcfout = out + "/" + gene + "." + region_out + ".vcf";

          ifstream vcffile (vcfout.c_str());
          for( std::string line; getline( vcffile, line ); )
          {
              if (line == "") {continue;}
              if (line.substr(0,2) == "##")
              {
                  if (count == 0) { head.push_back(line);}
                  continue;
              }
              if (line.substr(0,2) == "#C")
              {
                  boost::split(samples,line,boost::is_any_of("\t"));
                  sampleline = line;
                  continue;
              }
              vector <string> data;
              boost::split(data,line,boost::is_any_of("\t"));
              
              snp_data[stoi(data[1])] = data[0] + "\t" + data[1] + "\t" + data[2] + "\t" + data[3] + "\t" + data[4] + "\t" + data[5] + "\t" + data[6] + "\t" + data[7] + "\t" + "GT:DP:AD:PS";
              
              for (int a = 9; a < data.size(); a++)
              {
                  pair <int,string> k;
                  k = make_pair(stoi(data[1]),samples[a]);
                  if (data[a] == ".") {
                      vcf_data[k] = ".:.:.:.";
                  }
                  if (data[a] != ".") {
                      vector <string> fields;
                      boost::split(fields,data[a],boost::is_any_of(":"));
                      vcf_data[k] = fields[0] + ":" + fields[1] + ":" + fields[2] + ":.";
                  }
              }
          }
          vcffile.close();
          count++;
      }
      debug_message("Loading vcf ... end");
      
      debug_message("Combining vcf ... start");
      string outvcf = out + "/" + gene + ".combined.vcf";
      files_to_remove.push_back(outvcf);
      ofstream COMBVCF;
      COMBVCF.open (outvcf.c_str());
      
      for (auto item : head)
      {
          COMBVCF << item << endl;
      }
      COMBVCF << "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">" << endl;
      
      COMBVCF << sampleline << endl;

      for (auto item : snp_data)
      {
          int pos = item.first;
          string info = item.second;
          COMBVCF << info;
          
          for (int a = 9; a < samples.size(); a++)
          {
              
              string sample = samples[a];
              pair <int,string> k;
              k = make_pair(pos,sample);
              string value = vcf_data[k];
              
              if (sample.find("ref.") != std::string::npos) {
                  COMBVCF << "\t" << value;
                  continue;
              }
              pair <string,string> j;
              j = make_pair(sample,gene);
              string newgenotype = treat_genotype(value,2);
              COMBVCF << "\t" << newgenotype;
              continue;
          }
          COMBVCF << endl;
      }
      COMBVCF.close();

      debug_message("Combining vcf ... end");

      debug_message("Trimming vcf ... start");
      string outvcftrim =  out + "/" + gene + ".combined.trim.vcf";
      files_to_remove.push_back(outvcftrim);
      cmd = v_bcftools + " view --trim-alt-alleles " + outvcf + " -o " + outvcftrim;
      GetStdoutFromCommand(cmd.c_str());
      debug_message("Trimming vcf ... end");

      debug_message("Treating VCF ... start");
      string outvcftrimtreated = out + "/" + gene + ".combined.trim.treated.vcf";
      ifstream vcf (outvcftrim.c_str());
      ofstream TREATEDVCF;
      TREATEDVCF.open (outvcftrimtreated.c_str());
      for( std::string line; getline(vcf, line ); )
      {
          if (line == "") {continue;}
          if (line.substr(0,1) == "#")
          {
              TREATEDVCF << line << endl;
              continue;
          }
          vector <string> data;
          boost::split(data,line,boost::is_any_of("\t"));

          if (data[4] == ".") {continue;}
          if (data[4] == "*") {continue;}
          if (data[5] == "0") {continue;}
          TREATEDVCF << line << endl;
      }
      
      TREATEDVCF.close();
      vcf.close();
      debug_message("Treating VCF ... end");
      

    
    
    
    
    
    
        
    // reloading
    debug_message("Reloading VCF ... start");
    head.clear();
    samples.clear();
    snp_data.clear();
    vcf_data.clear();
    sampleline = "";
    float countmiss = 0;
    float countsnps = 0;
    map <string,float> refcountmiss;
    map <string,float> refcountvar;
    
    ifstream vcffile (outvcftrimtreated.c_str());
    for( std::string line; getline( vcffile, line ); )
    {
        if (line == "") {continue;}
        if (line.substr(0,2) == "##")
        {
            head.push_back(line);
            continue;
        }
        if (line.substr(0,2) == "#C")
        {
            boost::split(samples,line,boost::is_any_of("\t"));
            sampleline = line;
            continue;
        }
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));
        
        snp_data[stoi(data[1])] = data[0] + "\t" + data[1] + "\t" + data[2] + "\t" + data[3] + "\t" + data[4] + "\t" + data[5] + "\t" + data[6] + "\t" + data[7] + "\t" + data[8];
        
        for (int a = 9; a < data.size(); a++)
        {
            pair <int,string> k;
            k = make_pair(stoi(data[1]),samples[a]);
            vcf_data[k] = data[a];
            
            if (samples[a] == "Test")
            {
                countsnps++;
                vector <string> fields;
                boost::split(fields,data[a],boost::is_any_of(":"));
                if (fields[0].find(".") != std::string::npos) {countmiss++;}
            }

            if (samples[a] != "Test")
            {
                refcountvar[samples[a]]++;
                vector <string> fields;
                boost::split(fields,data[a],boost::is_any_of(":"));
                if (fields[0].find(".") != std::string::npos) {refcountmiss[samples[a]]++;}
            }
        }
    }
    vcffile.close();
    
    float missprop = countmiss / countsnps;
    if (missprop > 0.3)
    {
        
        for (auto item : files_to_remove) {removefile(item, v_debug);}
        return "failed!";
    }
    


            
    map <pair<string,int>,string> allele_data;
    map <string,int> allele_list;
    map <pair<string,string>,float> main_results;
    map <float,string> score_track;
            
    for (auto item : snp_data)
    {
        int pos = item.first;
        
        for (int a = 9; a < samples.size(); a++)
        {
            string sample = samples[a];
            if (sample.find("ref.") == std::string::npos) {continue;}

            float missratio = refcountmiss[sample] / refcountvar[sample];
            debug_message("Reference " + sample + " has " + to_string(missratio) + "  missing alleles.");
            if (type != "cds") {   if (missratio > 0.10) { continue; } }
            
            pair <int,string> k;
            k = make_pair(pos,sample);
            
            string info = vcf_data[k];
            
            vector <string> data;
            boost::split(data,info,boost::is_any_of(":"));
            if ((data[0] == ".") || (data[0] == "./.")) {continue;}
            
            pair <string,int> j;
            j = make_pair(sample,pos);
            allele_list[sample] = 1;
            allele_data[j] = data[0];
        }
    }
    debug_message("Reloading VCF ... end");        
            
            
    
    debug_message("Pairing ... start");
    ThreadPool poolsamples(stoi(threads));
    std::vector< std::future<int> > results_samples;
            
            
    int turn = 1;
    for (int a = 9; a < samples.size(); a++)
    {
        string sample = samples[a];
        if (sample.find("ref.") != std::string::npos) {continue;}
        
        pair <string,string> k;
        k = make_pair(sample,gene);
        turn++;
        
        
        results_samples.emplace_back(
            poolsamples.enqueue([turn,gene,sample,snp_data,&vcf_data,allele_list,&allele_data,&main_results,out,&score_track]
                {

                string outgenotypelist = "";
                outgenotypelist = out + "/" + gene + "." + sample + ".report.txt";
                ofstream REPORT;
                
                REPORT.open (outgenotypelist.c_str());
                REPORT << "SAMPLE\tGENE\tCOPY_NUMBER\tALLELE_1\tALLELE_2\tCOMPATIBILITY\tPOSITIVE_SNPS\tTOTAL_SNPS\tERRORS\tMISSINGS\n";

                
                map <int,string> sample_data;
                map <int,string> sample_position_ps;
                map <string,string> sample_ps_positions;
                for (auto item : snp_data)
                {
                    int pos = item.first;
                    pair <int,string> k;
                    k = make_pair(pos,sample);
                    string info = vcf_data[k];
                    vector <string> data;
                    boost::split(data,info,boost::is_any_of(":"));
                    if (data[0] == "") {continue;}
                    if (data[0] == ".") {continue;}
                    if (data[0] == "./.") {continue;}
                    if (data[0] == "././.") {continue;}
                    if (data[0] == "./././.") {continue;}
                    sample_data[pos] = data[0];
                    string ps = data[data.size()-1];
                    sample_position_ps[pos].append("," + ps);
                    sample_ps_positions[ps].append("," + to_string(pos));
                }
                    

                map <string,int> used;
                map <pair<float,string>,string> paired_results;
                
                for (auto alleleA : allele_list)
                {
                    for (auto alleleB : allele_list)
                    {
                        if (used.find(alleleB.first) != used.end()) {continue;}
                        
                        string paired = alleleA.first + "\t" + alleleB.first;
                        
                        int ncomparisons = 0;
                        int npositive = 0;
                        string errors = "";
                        string missings = "";
                        
                        for (auto item : sample_ps_positions)
                        {
                            string ps = item.first;
                            string poslist = item.second.substr(1);
                            vector <string> pos;
                            boost::split(pos,poslist,boost::is_any_of(","));
                            

                            for (auto position : pos)
                            {
                                pair <string,int> j;
                                j = make_pair(alleleA.first,stoi(position));
                                string nucA = ".";
                                if (allele_data.find(j) != allele_data.end())
                                {
                                    nucA = allele_data[j];
                                }
                                pair <string,int> k;
                                k = make_pair(alleleB.first,stoi(position));
                                string nucB = ".";
                                if (allele_data.find(k) != allele_data.end())
                                {
                                    nucB = allele_data[k];
                                }
                                
                                if ((nucA == ".") && (nucB == ".")) {continue;}
                                
                                string sample_genotype = "./.";
                                if (sample_data.find(stoi(position)) != sample_data.end())
                                {
                                    sample_genotype = sample_data[stoi(position)];
                                }
                                boost::replace_all(sample_genotype, "|", "/");
                                vector <string> sample_allele;
                                boost::split(sample_allele,sample_genotype,boost::is_any_of("/"));
                                
                                if (sample_genotype.find(".") != std::string::npos)
                                {
                                    missings.append("," + position);
                                }
                                
                                if ((sample_genotype == "./.") || (sample_genotype == "."))
                                {
                                    continue;
                                }
                                
                                unordered_map <string,int> valid_alleles;
                                valid_alleles[nucA]++;
                                valid_alleles[nucB]++;
                                
                                ncomparisons++;
                                
                                if ((nucA == ".") || (nucB == "."))
                                {
                                    if ((sample_allele[0] == ".") || (sample_allele[1] == "."))
                                    {
                                        npositive++;
                                        continue;
                                    }
                                }
                                
                                
                                if (sample_allele[0] == ".") {
                                    if (valid_alleles.find(sample_allele[1]) != valid_alleles.end())
                                    {
                                        npositive++;
                                    }
                                    else {
                                        errors.append("," + position);
                                    }
                                    continue;
                                }
                                
                                if (sample_allele[1] == ".") {
                                    if (valid_alleles.find(sample_allele[0]) != valid_alleles.end())
                                    {
                                        npositive++;
                                    }
                                    else {
                                        errors.append("," + position);
                                    }
                                    continue;
                                }
                                
                                if ((sample_allele[0] != ".") && (sample_allele[1] != "."))
                                {
                                    if ((nucA != ".") && (nucB != "."))
                                        {
                                            string testA = nucA + "/" + nucB;
                                            string testB = nucB + "/" + nucA;
                                            if ((sample_genotype == testA) || (sample_genotype == testB))
                                                {
                                                    npositive++;
                                            }
                                            else {
                                                errors.append("," + position);
                                            }
                                            continue;
                                    }
                                    
                                    if (nucA != ".")
                                        {
                                            if ((sample_allele[0] == nucA) || (sample_allele[1] == nucA))
                                                {
                                                    npositive++;
                                            }
                                            else {
                                                errors.append("," + position);
                                            }
                                            continue;
                                    }
                                    if (nucB != ".")
                                        {
                                            if ((sample_allele[0] == nucB) || (sample_allele[1] == nucB))
                                                {
                                                    npositive++;
                                            }
                                            else {
                                                errors.append("," + position);
                                            }
                                            continue;
                                    }
                                    
                            }
                        }

                        } // end loop positions
                        
                        
                        float ratio = float(npositive) / float(ncomparisons);
                        pair <float,string> j;
                        j = make_pair(ratio,paired);
                        paired_results[j] = to_string(npositive) + "\t" + to_string(ncomparisons) + "\t" + errors + "\t" + missings;
                    } // end pair
                    used[alleleA.first] = 1;
                } //end pair
                
                float higher_score = 0;
                for (auto iter = paired_results.rbegin(); iter != paired_results.rend(); ++iter) {
                    float score = iter->first.first;
                    string paired = iter->first.second;
                    string info = iter->second;
                    
                    if (higher_score < score) {higher_score = score;}
                    
                    
                    string subpaired = paired;
                    boost::replace_all(subpaired, "ref.", "");
                    
                    vector <string> data;
                    boost::split(data,info,boost::is_any_of("\t"));
                    string positives = data[0];
                    string comp = data[1];
                    string errors = data[2];
                    string miss = data[3];
                    if (errors != "") {errors = errors.substr(1);}
                    if (miss != "") {miss = miss.substr(1);}
                    if (errors == "") {errors = "none";}
                    if (miss == "") {miss= "none";}
                    REPORT << sample << "\t" << gene << "\t2\t" << subpaired << "\t" << score << "\t" << positives << "\t" << comp << "\t" << errors << "\t" << miss << endl;
                    
                    score_track[score].append(";" + subpaired);
                    
                    if (score == higher_score)
                    {
                        vector <string> suballeles;
                        boost::split(suballeles,paired,boost::is_any_of("\t"));
                        
                        string convertA = suballeles[0];
                        string convertB = suballeles[1];
                        boost::replace_all(convertA, "ref.", "");
                        boost::replace_all(convertB, "ref.", "");
                        
                        vector <string> split;
                        boost::split(split,convertA,boost::is_any_of("."));
                        convertA = split[0] + "*" + split[1].substr(0,5);
                        
                        mtxadj.lock();
                        string gen = convertA + "+" + convertB;
                        pair <string,string> k;
                        k = make_pair(sample,gen);
                        main_results[k] = score;
                        mtxadj.unlock();
                    }
                    
                    
                    
                }
                REPORT.close();
                return 1;

                })
        );
    } // final samples
            
     for(auto && result: results_samples){result.get();} // waiting for all threads
    debug_message("Pairing ... end");        
            
            ofstream final_report;
            string outfinal = "";
            outfinal = out + "/" + gene + ".calls.txt";
            
            final_report.open (outfinal.c_str());
            
            map <string,string> final_calls;
            map <string,float> final_scores;
            for (auto item : main_results)
            {
                string sample = item.first.first;
                string gen = item.first.second;
                float score = item.second;
                final_calls[sample].append(";" + gen);
                final_scores[sample] = score;
            }
            
            final_report << "SAMPLE\tALLELE_CALLS\tSNP_COMPATIBILITY" << endl;
            for (auto item : final_calls)
            {
                string sample = item.first;
                string gen = item.second;
                boost::replace_all(gen, ".", "*");
                
                pair <string,string> k;
                k = make_pair(sample,gene);
                final_report << sample << "\t" << "\t" << gen.substr(1) << "\t" << final_scores[sample] << endl;
            }
            
            final_report.close();
            

    select_alleles.clear();
    selected_alleles_list.clear();
    vector <string> selected_diplotypes;
    count = 0;
    
    for (auto iter = score_track.rbegin(); iter != score_track.rend(); ++iter)
    {
        float score = iter->first;
        string pair = iter->second;
        vector <string> data;
        string substr_pair = pair.substr(1);
        boost::split(data,substr_pair,boost::is_any_of(";"));
        for (auto sub : data)
        {
            selected_diplotypes.push_back(sub);
            vector <string> alleles;
            boost::split(alleles,sub,boost::is_any_of("\t"));
            select_alleles[alleles[0]] = 1;
            select_alleles[alleles[1]] = 1;
        }
        count++;
        if (count == 1) {break;}
    }
    
    
    debug_message("Selected dyplotypes after quick genotyping ...");

    if (debug == 1) {
        for (auto item : selected_diplotypes)
        {
            debug_message(item);
        }
    }
    
    if (gene != "KIR2DL5AB"){
        if (selected_diplotypes.size() >= 10) {return "failed";}
    }


    
    
    debug_message("Testing depth ... start");
    
    
    map <pair<string,int>,float> depth_track;
    map <pair<string,string>,string> read_track;
    map <string,string> head_track;
    map <string,int> reads;
    map <pair<string,string>,int> reads_nm_per_allele;
    
    ThreadPool cov_pool(stoi(threads));
    std::vector< std::future<int> > cov_results;
    for (auto item : fasfiles)
        {
            
            cov_results.emplace_back(
                cov_pool.enqueue([item,select_alleles,out,fq1,fq2,&files_to_remove,&depth_track,&head_track,&read_track,&reads,&reads_nm_per_allele]
                    {
                        
                        string file = item;
                        string allele = base_name(file);
                        boost::replace_all(allele, ".fas", "");
                        if (select_alleles.find(allele) == select_alleles.end()) {return 1;}
                        
                        string outsam = out + "/" + allele + ".sam";
                        string outsam2 = out + "/" + allele + ".2.sam";
                        
                        string cmd = "";
                        cmd = v_bwa + " mem -t 1 " + " -o " + outsam + " " +  file + " " + fq1;
                        GetStdoutFromCommand(cmd);
                        if (fq2 != "") {
                            string cmd = "";
                            cmd = v_bwa + " mem -t 1 " + " -o " + outsam2 + " " +  file + " " + fq2;
                            GetStdoutFromCommand(cmd);
                        }
                        
                        if (fq2 != "") {
                            ifstream FQ;
                            
                            FQ.open (outsam2 .c_str());
                            ofstream writer(outsam, ios::app);
                            for( std::string item; getline( FQ, item ); )
                            {
                                if (item == "") {continue;}
                                if (item.substr(0,1) == "@") {continue;}
                                vector <string> data;
                                boost::split(data,item,boost::is_any_of("\t"));
                                data[0].append(":2b");
                                string newline = "";
                                for (auto value : data) {newline.append("\t" + value);}
                                newline = newline.substr(1);
                                writer << newline << endl;
                            }
                            FQ.close();
                            writer.close();
                            mtxadj.lock();
                            files_to_remove.push_back(outsam);
                            files_to_remove.push_back(outsam2);
                            mtxadj.unlock();
                        }
                        
                        string outbam =  out + "/" + allele + ".bam";
                        string sam = out + "/" + allele + ".sam";
                        cmd = v_samtools +  " sort " + " -o " + outbam + " " + sam;
                        GetStdoutFromCommand(cmd);
                        
                        mtxadj.lock();
                        files_to_remove.push_back(outbam);
                        mtxadj.unlock();
                        
                        cmd = v_samtools + " depth " + outbam;
                        string depth = GetStdoutFromCommand(cmd);
                        vector <string> depthdata;
                        boost::split(depthdata,depth,boost::is_any_of("\n"));
                        
                        for (auto item : depthdata)
                        {
                            if (item == "") {continue;}
                            if (item.substr(0,1) == "[") {continue;}
                            vector <string> data;
                            boost::split(data,item,boost::is_any_of("\t"));
                            pair <string,int> key = make_pair(allele,stoi(data[1]));
                            mtxadj.lock();
                            depth_track[key] = stof(data[2]);
                            mtxadj.unlock();
                        }
                        
                        ifstream SAM;
                        SAM.open (outsam .c_str());
                        for( std::string item; getline( SAM, item ); )
                        {
                            if (item == "") {continue;}
                            if (item.substr(0,1) == "[") {continue;}
                            if (item.substr(0,2) == "@P") {continue;}
                            if (item.substr(0,2) == "@S") {
                                mtxadj.lock();
                                head_track[allele].append(item + "\n");
                                mtxadj.unlock();
                                continue;
                            }
                            vector <string> data;
                            boost::split(data,item,boost::is_any_of("\t"));
                            if (data[2] == "*") {continue;}
                            pair <string,string> key = make_pair(allele,data[0]);
                            mtxadj.lock();
                            read_track[key] = item;
                            mtxadj.unlock();
                            int nm = 100;
                            if (data.size() >= 11)
                                {
                                    string nmfield = data[11];
                                    boost::replace_all(nmfield, "NM:i:", "");
                                    nm = stoi(nmfield);
                                }
                            mtxadj.lock();
                            reads_nm_per_allele[key] = nm;
                            reads[data[0]] = 1;
                            mtxadj.unlock();
                        }
                        SAM.close();
                        
                        return 1;
                    })
            );
        }
    for(auto && result: cov_results){result.get();}
    debug_message("Testing depth ... end");
    
    
    
    
       
    debug_message("Pairing");
    

    map <pair<string,string>,float> scores;
    map <string,int> used;

//    ThreadPool pair_pool(1);
    ThreadPool pair_pool(stoi(v_threads));
    std::vector< std::future<int> > pair_results;
    
    
    for (auto diplotype : selected_diplotypes)
    {
        vector <string> alleles;
        boost::split(alleles,diplotype,boost::is_any_of("\t"));
        string alleleA = alleles[0];
        string alleleB = alleles[1];
        string pair_of_alleles = diplotype;
        

                
        pair_results.emplace_back(
            pair_pool.enqueue([alleleA, alleleB, out, &files_to_remove, &head_track,&reads, &read_track, minnm, db, gene, type, &depth_track, pair_of_alleles, &scores]
                {
                    
            string outsam = out + "/" + alleleA + "_and_" + alleleB + ".tmp.sam";
            string outsam2 = out + "/" + alleleA + "_and_" + alleleB + ".sam";
            mtxadj.lock();
            files_to_remove.push_back(outsam2);
            files_to_remove.push_back(outsam);
            mtxadj.unlock();
                    
                    float depth_sum = 0;
                    float depth_count = 0;
                    for (auto item : depth_track)
                    {
                        string allele = item.first.first;
                        int pos = item.first.second;
                        float depth = item.second;
                        if (depth > 0) {
                            if ((allele == alleleA) || (allele == alleleB))
                                {
                                    depth_sum = depth_sum + depth;
                                    depth_count++;
                                }
                        }
                    }
                    
                    float coverage = depth_sum / depth_count;
                    float match_level_adapted = coverage / 4;
                    
                    ofstream outfile;
                    outfile.open (outsam.c_str());
                    outfile << head_track[alleleA];
                    if (alleleB != alleleA) {
                        outfile << head_track[alleleB];
                    }
                    
                    float included_reads = 0;
                    float duplicated_reads = 0;
                    
                    
                    for (auto item : reads)
                    {
                        string read = item.first;
                        pair <string,string> keyA = make_pair(alleleA,read);
                        string lineA = read_track[keyA];
                        pair <string,string> keyB = make_pair(alleleB,read);
                        string lineB = read_track[keyB];
                        
                        
                        vector <string> dataA;
                        vector <string> dataB;
                        boost::split(dataA,lineA,boost::is_any_of("\t"));
                        boost::split(dataB,lineB,boost::is_any_of("\t"));
                        
                        int nmA = 100;
                        if (dataA.size() >= 11)
                        {
                            string nmfield = dataA[11];
                            boost::replace_all(nmfield, "NM:i:", "");
                            nmA = stoi(nmfield);
                        }
                        int nmB = 100;
                        if (dataB.size() >= 11)
                        {
                            string nmfield = dataB[11];
                            boost::replace_all(nmfield, "NM:i:", "");
                            nmB = stoi(nmfield);
                        }
                        if (alleleA == alleleB) {
                            outfile << lineA << endl;
                            included_reads++;
                        }
                        
                        if (alleleA != alleleB) {
                            if ((nmA <= minnm) && (nmA == nmB)) {
                                dataA[0].append("-1");
                                lineA = "";
                                for (auto value : dataA) {lineA.append("\t" + value);}
                                lineA = lineA.substr(1);
                                outfile << lineA << endl;
                                outfile << lineB << endl;
                                included_reads++;
                                duplicated_reads++;
                            }
                            
                            if ((nmA <= minnm) && (nmA < nmB)) {
                                outfile << lineA << endl;
                                included_reads++;
                            }
                            if ((nmB <= minnm) && (nmB < nmA)) {
                                outfile << lineB << endl;
                                included_reads++;
                            }
                        }
                    }
                    outfile.close();
                    
                    
                    ifstream t(outsam);
                    ofstream o(outsam2);
                    for( std::string item; getline( t, item ); )
                    {
                        if (item == "") {continue;}
                        o << item << endl;
                    }
                    t.close();
                    o.close();
                    
                    
                    string outbam = out + "/" + alleleA + "_and_" + alleleB + ".bam";
                    string cmd = v_samtools + " sort -o " + outbam + " " + outsam2;
                    GetStdoutFromCommand(cmd);
                    cmd = v_samtools + " index " + outbam;
                    GetStdoutFromCommand(cmd);
                    string fasfileA = db + "/adjustment/dna/fasta/" + type + "/" + gene + "/" + alleleA + ".fas";
                    string fasfileB = db + "/adjustment/dna/fasta/" + type + "/" + gene + "/" + alleleB + ".fas";
                    
                    mtxadj.lock();
                    files_to_remove.push_back(outbam);
                    files_to_remove.push_back(outbam + ".bai");
                    files_to_remove.push_back(outsam);
                    mtxadj.unlock();
                    
                    cmd = v_samtools + " mpileup -a -r " + alleleA + " -f " + fasfileA + " " + outbam;
                    string pileA = GetStdoutFromCommand(cmd);
                    cmd = v_samtools + " mpileup -a -r " + alleleB + " -f " + fasfileB + " " + outbam;
                    string pileB = GetStdoutFromCommand(cmd);
                    
                    
                    vector <string> pileAdata;
                    boost::split(pileAdata,pileA,boost::is_any_of("\n"));
                    vector <string> pileBdata;
                    boost::split(pileBdata,pileB,boost::is_any_of("\n"));
                    
                    float testedA = 0;
                    float validA = 0;
                    float sumcovA = 0;
                    float noise = 0;
                    
                    for (auto item : pileAdata)
                    {
                        if (item == "") {continue;}
                        if (item.substr(0,1) == "[") {continue;}
                        vector <string> data;
                        boost::split(data,item,boost::is_any_of("\t"));
                        if (data.size() < 5) {continue;}
                        string bases = data[4];
                        if (bases == "*") {continue;}
                        if (data[2] == "N") {continue;}
                        if (data[1] == "") {continue;}
                        if (data[1] == "*") {continue;}
                        int pos = stoi(data[1]);
                        float cov = stof(data[3]);
                        
                        float refcount = 0;
                        for (int a = 0; a < bases.size(); a++)
                        {
                        if ((bases.substr(a,1) == ".") || (bases.substr(a,1) == ","))
                            {
                                refcount++;
                            }
                        }
                    
                        noise = noise + (cov - refcount);
                        pair <string,int> key = make_pair(alleleA,pos);
                        
                        if (depth_track[key] >= match_level_adapted)
                        {
                            testedA++;
                            if (refcount >= match_level_adapted)
                            {
                                validA++;
                                sumcovA = sumcovA + refcount;
                            }
                        }
                    }
                    
                    float noiseA = 100;
                    if (testedA > 0) {noiseA = noise / testedA;}
                    
                    float compA = 0;
                    if (testedA > 0) {compA = validA / testedA;}
                    
                    float meancovA = 0;
                    if (testedA > 0) {meancovA = sumcovA / testedA;}
                    
                    
                    
                    float testedB = 0;
                    float validB = 0;
                    float sumcovB = 0;
                    noise = 0;
                    
                    for (auto item : pileBdata)
                    {
                        if (item == "") {continue;}
                        if (item.substr(0,1) == "[") {continue;}
                        vector <string> data;
                        boost::split(data,item,boost::is_any_of("\t"));
                        if (data.size() < 5) {continue;}
                        string bases = data[4];
                        if (bases == "*") {continue;}
                        if (data[2] == "N") {continue;}
                        if (data[1] == "") {continue;}
                        if (data[1] == "*") {continue;}
                        
                        int pos = stoi(data[1]);
                        float cov = stof(data[3]);
                        
                        float refcount = 0;
                        for (int a = 0; a < bases.size(); a++)
                        {
                            if ((bases.substr(a,1) == ".") || (bases.substr(a,1) == ","))
                            {
                                refcount++;
                            }
                        }
                        
                        noise = noise + (cov - refcount);
                        pair <string,int> key = make_pair(alleleB,pos);
                        
                        if (depth_track[key] >= match_level_adapted)
                        {
                            testedB++;
                            if (refcount >= match_level_adapted)
                            {
                                validB++;
                                sumcovB = sumcovB + refcount;
                            }
                        }
                    }
                    
                    float noiseB = 100;
                    if (testedB > 0) {noiseB = noise / testedB;}
                    
                    float compB = 0;
                    if (testedB > 0) {compB = validB / testedB;}
                    
                    float meancovB = 0;
                    if (testedB > 0) {meancovB = sumcovB / testedB;}
                    
                    float duplication_level = duplicated_reads / included_reads;
                    if (alleleA == alleleB) {duplication_level = 1;}
                    
                    
                    mtxadj.lock();
                    pair <string,string> key = make_pair(pair_of_alleles,"compA");
                    scores[key] = compA;
                    key = make_pair(pair_of_alleles,"compB");
                    scores[key] = compB;
                    key = make_pair(pair_of_alleles,"testA");
                    scores[key] = testedA;
                    key = make_pair(pair_of_alleles,"testB");
                    scores[key] = testedB;
                    key = make_pair(pair_of_alleles,"covA");
                    scores[key] = meancovA;
                    key = make_pair(pair_of_alleles,"covB");
                    scores[key] = meancovB;
                    key = make_pair(pair_of_alleles,"noiseA");
                    scores[key] = noiseA;
                    key = make_pair(pair_of_alleles,"noiseB");
                    scores[key] = noiseB;
                    key = make_pair(pair_of_alleles,"dup");
                    scores[key] = duplication_level;
                    mtxadj.unlock();
                    
                    
                    return 1;
                })
        );
            
        }
    for(auto && result: pair_results){result.get();}
    
    debug_message("Pairing ... end");


    map <float,string> final_data;
    
    if (debug == 1) {
        cout << "\nPairing results:\n";
        cout << "ALLELE_A\tALLELEB\tCOMPATIBILITY_A\tCOMPATIBILITY_B\tCOV_A\tCOV_B\tDUP_LEVEL\tNOISE_A\tNOISE_B\tSCORE\tCOMP\tCOV\tNOISE\n";
    }
    
    
    map <string,int> skip;
    for (auto item : scores)
        {
            string pair_of_alleles = item.first.first;
            if (skip.find(pair_of_alleles) != skip.end()) {continue;}
            skip[pair_of_alleles] = 1;
            vector <string> alleles;
            boost::split(alleles,pair_of_alleles,boost::is_any_of("\t"));
            
            float score = 0;
            float comp = 0;
            float cov = 0;
            
            pair <string,string> key = make_pair(pair_of_alleles,"noiseA");
            float noiseA = scores[key];
            key = make_pair(pair_of_alleles,"noiseB");
            float noiseB = scores[key];
            key = make_pair(pair_of_alleles,"compA");
            float compA = scores[key];
            key = make_pair(pair_of_alleles,"compB");
            float compB = scores[key];
            key = make_pair(pair_of_alleles,"covA");
            float covA = scores[key];
            key = make_pair(pair_of_alleles,"covB");
            float covB = scores[key];
            key = make_pair(pair_of_alleles,"dup");
            float dup = scores[key];
            
            
            float noise = 1 - (((noiseA + noiseB) / 2) / 10);
            
            comp = ((compA + compB)/2);
            cov = (covA + covB) - ((covA + covB) * dup) / 2;
            
            score = cov * comp * noise;
            
            final_data[score].append("," + pair_of_alleles);
            
            if (debug == 1)
                {
                    cout << pair_of_alleles;
                    cout << "\t" << compA;
                    cout << "\t" << compB;
                    cout << "\t" << covA;
                    cout << "\t" << covB;
                    cout << "\t" << dup;
                    cout << "\t" << noiseA;
                    cout << "\t" << noiseB;
                    cout << "\t" << score;
                    cout << "\t" << comp;
                    cout << "\t" << cov;
                    cout << "\t" << noise;
                    cout << endl;
                    
                }
        }
    
    string valid_pair = "";
    for (auto iter = final_data.rbegin(); iter != final_data.rend(); ++iter) {
        
        string all = iter->second.substr(1);
        vector <string> list;
        boost::split(list,all,boost::is_any_of(","));
        valid_pair = list[0];
        break;
    }
    
    pair <string,string> key = make_pair(valid_pair,"compA");
    float compA = scores[key];
    key = make_pair(valid_pair,"compB");
    float compB = scores[key];
    float comp = ((compA + compB)/2);
    if (debug == 1)
    {
        cout << endl << "Valid pair: " << valid_pair << endl;
    }
    vector <string> alleles;
    boost::split(alleles,valid_pair,boost::is_any_of("\t"));
    
    map <pair<string,string>,float> read_db_for;
    map <pair<string,string>,float> read_db_rev;
    
    for( auto item : reads_nm_per_allele)
        {
            string allele = item.first.first;
            string read = item.first.second;
            string subread = read;
            int nm = item.second;
            if ((allele == alleles[0]) || (allele == alleles[1]))
                {
                    boost::replace_all(subread, ":2b-1", "");
                    boost::replace_all(subread, ":2b", "");
                    boost::replace_all(subread, "-1", "");
                    pair <string,string> key = make_pair(subread,allele);
                    if (read.find(":2b") == std::string::npos)
                        {
                            read_db_for[key] = nm;
                        }
                    else {
                        read_db_rev[key] = nm;
                    }
                }
        }
    
    map <string,string> final_db;
    map <string,int> used_read;
    
    ofstream outfile;
    string outreport = out + "/" + gene + ".readlist.txt";
    outfile.open (outreport.c_str());
    used.clear();
    map <string,int> track_read_nm;
    for (auto item : read_db_for)
    {
        string read = item.first.first;
        if (used_read.find(read) != used_read.end()) {continue;}
        
        pair <string,string> key = make_pair(read,alleles[0]);
        int subAf = 100; if (read_db_for.find(key) != read_db_for.end()) {subAf = read_db_for[key];}
        int subAr = 100; if (read_db_rev.find(key) != read_db_rev.end()) {subAr = read_db_rev[key];}
        key = make_pair(read,alleles[1]);
        int subBf = 100; if (read_db_for.find(key) != read_db_for.end()) {subBf = read_db_for[key];}
        int subBr = 100; if (read_db_rev.find(key) != read_db_rev.end()) {subBr = read_db_rev[key];}
        
        if ((subAf == 100) && (subAr < 100)) {subAf = subAr;}
        if ((subAr == 100) && (subAf < 100)) {subAr = subAf;}
        
        if ((subBf == 100) && (subBr < 100)) {subBf = subBr;}
        if ((subBr == 100) && (subBf < 100)) {subBr = subBf;}
        
        int valid_for = 100;
        if (subAf <= subBf) {valid_for = subAf;}
        if (subAf > subBf) {valid_for = subBf;}
        
        int valid_rev = 100;
        if (subAr <= subBr) {valid_rev = subAr;}
        if (subAr > subBr) {valid_rev = subBr;}
        
        outfile << read << "\tfor\t" << valid_for << endl;
        outfile << read << "\trev\t" << valid_rev << endl;
        
        int sum = valid_for + valid_rev;
        track_read_nm[read] = sum;
        
        used[read] = 1;
    }
    outfile.close();
    
    mtxadj.lock();
    files_to_remove.push_back(samtmp);
    mtxadj.unlock();
    
    for (auto item : files_to_remove)
    {
        removefile(item,v_debug);
    }
    
    string readfile = out + "/" + gene + ".readlist.txt";
    if (! fileExists(readfile)) {return "failed";}
    
    
    map <string,int> max_mismatch;

    max_mismatch["KIR2DL1"] = 1;
    max_mismatch["KIR2DL2"] = 0;
    max_mismatch["KIR2DL3"] = 0;
    max_mismatch["KIR2DL4"] = 1;
    max_mismatch["KIR2DL5AB"] = 1;
    max_mismatch["KIR2DP1"] = 1;

    max_mismatch["KIR2DS1"] = 0;
    max_mismatch["KIR2DS2"] = 0;
    max_mismatch["KIR2DS3"] = 0;
    max_mismatch["KIR2DS4"] = 0;
    max_mismatch["KIR2DS5"] = 0;

    max_mismatch["KIR3DL1"] = 0;
    max_mismatch["KIR3DL2"] = 1;
    max_mismatch["KIR3DL3"] = 1;
    max_mismatch["KIR3DS1"] = 0;
 


    for (auto item : track_read_nm) {
        int max = max_mismatch[gene];
        if (item.second <= max) {
            pair <string,string> key = make_pair(gene,item.first);
            reads_possibly_mapped[key] = item.second;
        }
    }

    return alleles[0] + "\t" + alleles[1];

}



string adjust(string gene, string r1, string r2, string out, string nthreads)
{
    string return_res = "failed";
    
    string script = v_db + "/adjustment/dna/adjust.pl";
    
    string fas = v_db + "/adjustment/dna/fasta/";

    if (v_exome == 1) {
        fas.append("/cds/" + gene);
    }
    else {
        fas.append("/full/" + gene);
    }
    
    string command = "perl " + script;
    command.append(" -g " + gene);
    command.append(" -f " + r1);
    if (r2 != "") {
        command.append(" -r " + r2);
    }
    command.append(" -o " + out);
    command.append(" -t " + nthreads);

    if (v_debug == 1)
    {
        command.append(" -d 1");
    }
    
    if (v_exome == 1)
    {
        command.append(" -x 'cds'");
    }
    if (v_exome == 0)
    {
        command.append(" -x 'full'");
    }
    command.append(" -y " + fas);
    
//    cout << endl << command << endl;
    
    string outscript = GetStdoutFromCommand(command);
    if (outscript.find("ailed!") != std::string::npos) {
        return "failed";
    }
    if (outscript == "") {return "failed";}
    
    
    
    string readfile = out + "/" + gene + ".readlist.txt";
    if (! fileExists(readfile)) {return "failed";}
    
    ifstream readdb;
    readdb.open (readfile.c_str());
    for( std::string item; getline( readdb, item ); )
    {
        if (item == "") {continue;}
        vector <string> data;
        boost::split(data,item,boost::is_any_of("\t"));
        string read = data[0];
        string nm = data[1];
        
        mtx_type.lock();
        pair <string,string> key = make_pair(gene,data[0]);
        sequence_list_r1[key] = stoi(nm) + 1;
        sequence_list_r2[key] = stoi(nm) + 1;
        mtx_type.unlock();
        
        
        if (stoi(nm) == 0) {
            pair <string,string> key = make_pair(gene,read);
            reads_possibly_mapped[key] = stoi(nm);
        }
    }
    
    return outscript;
}


std::map<std::string, std::string> TypeCDS( string gene, string r1, string r2)
{
    
    map <string,string> return_res;
    
    return_res["alleleA"] = "NA";
    return_res["alleleB"] = "NA";
    
    vector <string> valid_alleles;
    
    map <pair<string,string>,int> motifs;
    map <string,float> motif_count;
    map <string,float> allele_list_motif;
    ifstream motif;
    string file = v_db + "/adjustment/cds_motifs/" + gene + ".txt";
    motif.open (file.c_str());
    for( std::string line; getline( motif, line ); )
    {
        if (line == "") {continue;}
        vector <string> fields;
        boost::split(fields,line,boost::is_any_of("\t"));
        vector <string> list;
        boost::split(list,fields[1],boost::is_any_of(","));
        for (auto item : list)
        {
            pair <string,string> key = make_pair(fields[0],item);
            motifs[key] = 1;
            motif_count[item] = 0;
            allele_list_motif[fields[0]] = 1;
        }
    }
    motif.close();
    
    ifstream r1file;
    r1file.open (r1.c_str());
    for( std::string line; getline( r1file, line ); )
    {
        string id = line;
        getline( r1file, line );
        string seq = line;
        getline( r1file, line );
        string info = line;
        getline( r1file, line );
        string qual = line;

        for (int a = 0; a < (seq.size() - 15); a++)
        {
            string sub = seq.substr(a,15);
            if (motif_count.find(sub) != motif_count.end())
            {
                motif_count[sub]++;
            }
        }
    }
    r1file.close();
    
    for (auto item :  allele_list_motif)
    {
        string allele = item.first;
        float valid = 0;
        float comp = 0;
        for (auto list : motif_count)
        {
            string motif = list.first;
            float count = list.second;
            pair <string,string> key = make_pair(allele,motif);
            if (motifs.find(key) != motifs.end())
            {
                if (count > 2) {valid++;}
                comp++;
            }
        }
        float score = valid / comp;
        if (score > 0.7) {valid_alleles.push_back(allele);}
    }
    
    if (valid_alleles.size() == 0) {return return_res;}
    
    
    
    
    string cmd = "mkdir " + v_output + "/adjustment/";
    GetStdoutFromCommand(cmd.c_str());
    cmd = "mkdir " + v_output + "/adjustment/" + gene;
    GetStdoutFromCommand(cmd.c_str());
    
    

    
    // CHECKING TARGETS
    
    map <string,float> sequence_size;
    unordered_map <string,string> sequence_content;
    unordered_map <string,float> sequence_frags;
    string reference = v_db + "/adjustment/cds/" + gene + ".fas";
    ifstream fasta;
    fasta.open (reference.c_str());
    for( std::string line; getline( fasta, line ); )
    {
        string id = line;
        getline( fasta, line );
        string seq = line;
        boost::replace_all(id, ">", "");
        sequence_content[id] = seq;
        vector <string> seqsplit;
        boost::split(seqsplit,seq,boost::is_any_of("N"));
        float count = 0;
        for (auto item : seqsplit)
        {
            if (item == "") {continue;}
            count++;
        }
        sequence_frags[id] = count;
        boost::replace_all(seq, "N", "");
        sequence_size[id] = float(seq.length());
    }
    fasta.close();
    
    if (sequence_size.size() == 0) {return return_res;}
    
    



    /*
    
    
    // SELECTING
    
    
    ThreadPool pool_first_round(stoi(v_threads));
    std::vector< std::future<int> > results_first_round;
    map <float,string> allele_first_round;
    int loop = 1;

    for (auto item : sequence_size)
    {
        string allele = item.first;
        
        results_first_round.emplace_back(
            pool_first_round.enqueue([loop, allele, gene, &allele_first_round,r1,r2]
                {
                    
                    string tmp_ref = v_db + "/adjustment/cds/" + gene + "/" + allele + ".fas";
                    string log = v_output + "/adjustment/" + gene + "/tmp_" + allele + ".log";;
                    
                    string cmd = v_bwa + " mem -B2 -L 1 " + tmp_ref + " " + r1 + " " + " 2> " + log + " | " + v_samtools + " sort - | " + v_samtools + " mpileup -f " + tmp_ref + " -a - 2>" + log;
                    string pileA = GetStdoutFromCommand(cmd.c_str());
                    cmd = v_bwa + " mem -B2 -L 1 " + tmp_ref + " " + r2 + " " + " 2> " + log + " | " + v_samtools + " sort - | " + v_samtools + " mpileup -f " + tmp_ref + " -a - 2>" + log;
                    string pileB = GetStdoutFromCommand(cmd.c_str());
                    
                    
                    vector <string> piledataA;
                    boost::split(piledataA,pileA,boost::is_any_of("\n"));
                    vector <string> piledataB;
                    boost::split(piledataB,pileB,boost::is_any_of("\n"));

    
                    float valid = 0;
                    float comparisons = 0;
                    
                    for (auto line : piledataA)
                    {
                        if (line == "") {continue;}
                        if (line.substr(0,1) == "[") {continue;}
                        vector <string> fields;
                        boost::split(fields,line,boost::is_any_of("\t"));
                        if (fields[2] == "N") {continue;}
                        if (fields[3] == "0") {continue;}
                        string allele = fields[0];
                        string bases = fields[4];
                        int x = boost::count(bases, ',');
                        int y = boost::count(bases, '.');
                        float refcov = float(x) + float(y);
                        if ((refcov >= 4))
                        {
                            valid++;
                        }
                        comparisons++;
                    }
            
                    for (auto line : piledataB)
                    {
                        if (line == "") {continue;}
                        if (line.substr(0,1) == "[") {continue;}
                        vector <string> fields;
                        boost::split(fields,line,boost::is_any_of("\t"));
                        if (fields[2] == "N") {continue;}
                        if (fields[3] == "0") {continue;}
                        string allele = fields[0];
                        string bases = fields[4];
                        int x = boost::count(bases, ',');
                        int y = boost::count(bases, '.');
                        float refcov = float(x) + float(y);
                        if ((refcov >= 4))
                        {
                            valid++;
                        }
                        comparisons++;
                    }
                    
                    float score = valid / comparisons;
                    mtx_type.lock();
                    allele_first_round[score].append("," + allele);
                    mtx_type.unlock();
                    removefile(log, 0);
                    return 1;
                })
        );
        loop++;
    }

    for(auto && result: results_first_round){result.get();} // waiting for all threads
    
    
    vector <string> valid_alleles;
    int round = 0;
    for (auto iter = allele_first_round.rbegin(); iter != allele_first_round.rend(); ++iter) {
        vector <string> alleles;
        if (iter->first < 0.95) {continue;}
        boost::split(alleles,iter->second,boost::is_any_of(","));
        for (auto item : alleles)
        {
            if (item == "") {continue;}
            valid_alleles.push_back(item);
        }
        round++;
        if (valid_alleles.size() > 10) {break;}
        if (round == 2) {break;}
    }

    if (valid_alleles.size() < 2) {return return_res;}
    
  */
    
    
    
// PAIRING
    
    
    map <float,string> final_scores;
    string comparison_data = "";
    
    ThreadPool pool_pairs(stoi(v_threads));
    std::vector< std::future<int> > results_pairs;
    int loop = 1;
    map <string,int> used;
    for (auto item : valid_alleles)
    {
        string alleleA = item;
        for (auto sub : valid_alleles)
        {
            string alleleB = sub;
            if (used.find(alleleB) != used.end()) {continue;}
            
            results_pairs.emplace_back(
                pool_pairs.enqueue([loop, gene, alleleA, alleleB, &final_scores, &comparison_data,&sequence_size,&sequence_frags,&sequence_content,r1,r2]
                    {
                        
                        ofstream tmp_ref_out;
                        string rev = v_output + "/adjustment/" + gene + "/tmp_" + alleleA + "_" + alleleB + ".fas";
                        tmp_ref_out.open (rev.c_str());
                        tmp_ref_out << ">" << alleleA << endl;
                        tmp_ref_out << sequence_content[alleleA] << endl;
                        if (alleleA != alleleB)
                        {
                            tmp_ref_out << ">" << alleleB << endl;
                            tmp_ref_out << sequence_content[alleleB] << endl;
                        }
                        tmp_ref_out.close();
                        
                        string cmd = v_bwa + " index " + rev;
                        GetStdoutFromCommand(cmd.c_str());
                        
                        string bam = v_output + "/adjustment/" + gene + "/tmp_" + alleleA + "_" + alleleB + ".bam";;
                        string log = v_output + "/adjustment/" + gene + "/tmp_" + alleleA + "_" + alleleB + ".log";;
                        cmd = v_bwa + " mem -B2 -L 1 " + rev + " " + r1 + " " + " 2> " + log + " | " + v_samtools + " sort - | " + v_samtools + " mpileup -f " + rev + " -a - 2>" + log;
                        string pileA = GetStdoutFromCommand(cmd.c_str());
                        cmd = v_bwa + " mem -B2 -L 1 " + rev + " " + r2 + " " + " 2> " + log + " | " + v_samtools + " sort - | " + v_samtools + " mpileup -f " + rev + " -a - 2>" + log;
                        string pileB = GetStdoutFromCommand(cmd.c_str());
                        vector <string> piledataA;
                        boost::split(piledataA,pileA,boost::is_any_of("\n"));
                        vector <string> piledataB;
                        boost::split(piledataB,pileB,boost::is_any_of("\n"));
                        
                        
                        map <string,float> valid;
                        map <string,float> comparisons;
                        map <string,float> covsum;
                        map <string,float> noise;
                        
                        for (auto line : piledataA)
                        {
                            if (line == "") {continue;}
                            if (line.substr(0,1) == "[") {continue;}
                            vector <string> fields;
                            boost::split(fields,line,boost::is_any_of("\t"));
                            if (fields[2] == "N") {continue;}
                            string bases = fields[4];
                            float cov = stof(fields[3]);
                            string allele = fields[0];
                            
                            if (cov < 5) {continue;}
                            covsum[allele] = covsum[allele] + cov;
                            int x = boost::count(bases, ',');
                            int y = boost::count(bases, '.');
                            float refcov = float(x) + float(y);
                            if ((refcov >= 5))
                            {
                                valid[allele]++;
                            }
                            noise[allele] = noise[allele] + (cov - refcov);
                            comparisons[allele]++;
                        }
                        
                        for (auto line : piledataB)
                        {
                            if (line == "") {continue;}
                            if (line.substr(0,1) == "[") {continue;}
                            vector <string> fields;
                            boost::split(fields,line,boost::is_any_of("\t"));
                            if (fields[2] == "N") {continue;}
                            string bases = fields[4];
                            float cov = stof(fields[3]);
                            string allele = fields[0];
                            
                            if (cov < 5) {continue;}
                            covsum[allele] = covsum[allele] + cov;
                            int x = boost::count(bases, ',');
                            int y = boost::count(bases, '.');
                            float refcov = float(x) + float(y);
                            if ((refcov >= 5))
                            {
                                valid[allele]++;
                            }
                            noise[allele] = noise[allele] + (cov - refcov);
                            comparisons[allele]++;
                        }
                    
                    
                    
                    
                        float scoreA = valid[alleleA] / comparisons[alleleA];
                        float errorA = noise[alleleA] / covsum[alleleA];
                        float scoreB = valid[alleleB] / comparisons[alleleB];
                        float errorB = noise[alleleB] / covsum[alleleB];
                        float covA = covsum[alleleA] / comparisons[alleleA];
                        float covB = covsum[alleleB] / comparisons[alleleB];
                        
                        float comp = ((scoreA + scoreB) / 2);
                        float cov = ((covA + covB) / 2);
                        if (alleleA == alleleB) {cov = cov / 2;}
                        float compadj = pow(comp,10);
                        float error = 1-((errorA + errorB) / 2);
                        float erroradj = pow(error,10);
//                        float score = compadj * cov * erroradj;
                        float score = compadj * erroradj;
                        
                        final_scores[score].append("," + alleleA + "\t" + alleleB);
                        comparison_data.append(alleleA + "\t" + alleleB + "\t" + to_string(scoreA) + "\t" + to_string(scoreB) + "\t" + to_string(cov) + "\t" + to_string(score) + "\t" + to_string(errorA) + "\t" + to_string(errorB) + "\n");
                        mtx_type.unlock();
                        
                        removefile(rev, 0);
                        removefile(log, 0);
                        removefile(rev + ".sa", 0);
                        removefile(rev + ".pac", 0);
                        removefile(rev + ".fai", 0);
                        removefile(rev + ".bwt", 0);
                        removefile(rev + ".ann", 0);
                        removefile(rev + ".amb", 0);
                        
                        return 1;
                    })
            );
            loop++;
        }
        used[alleleA] = 1;
    }
    for(auto && result: results_pairs){result.get();} // waiting for all threads
    
    
    
    string paired = "";
    for (auto iter = final_scores.rbegin(); iter != final_scores.rend(); ++iter) {
        string pairs = iter->second.substr(1);
        paired = pairs;
        return_res["score"] = to_string(iter->first);
        return_res["pairs"] = paired;
        break;
    }

    vector <string> list_of_pairs;
    boost::split(list_of_pairs,paired,boost::is_any_of(","));
    vector <string> alleles;
    boost::split(alleles,list_of_pairs[0],boost::is_any_of("\t"));
    return_res["alleleA"] = alleles[0];
    return_res["alleleB"] = alleles[1];
    return_res["comparison_data"] = comparison_data;

    

    
    map <string,int> reads_to_primary_A;
    string tmp_ref = v_db + "/adjustment/cds/" + gene + "/" + alleles[0] + ".fas";
    string log = v_output + "/adjustment/" + gene + "/tmp_" + alleles[0] + ".log";;
    cmd = v_bwa + " mem -B2 -L 1 " + tmp_ref + " " + r1 + " " + " 2> " + log;
    string sam = GetStdoutFromCommand(cmd.c_str());
    vector <string> data;
    boost::split(data,sam,boost::is_any_of("\n"));
    for (auto item : data)
    {
        if (item == "") {continue;}
        if (item.substr(0,1) == "[") {continue;}
        vector <string> fields;
        boost::split(fields,item,boost::is_any_of("\t"));
        string read = fields[0];
        if (fields[1] == "*") {continue;}
        if (fields[2] == "*") {continue;}
        if (fields.size() >= 11)
        {
            string nm = fields[11];
            boost::replace_all(nm, "NM:i:", "");
            if (reads_to_primary_A.find(read) != reads_to_primary_A.end())
            {
                reads_to_primary_A[read] = reads_to_primary_A[read] + stoi(nm);
            }
            if (reads_to_primary_A.find(read) == reads_to_primary_A.end())
            {
                reads_to_primary_A[read] = stoi(nm);
            }
        }
    }
    tmp_ref = v_db + "/adjustment/cds/" + gene + "/" + alleles[0] + ".fas";
    log = v_output + "/adjustment/" + gene + "/tmp_" + alleles[0] + ".log";;
    cmd = v_bwa + " mem -B2 -L 1 " + tmp_ref + " " + r2 + " " + " 2> " + log;
    sam = GetStdoutFromCommand(cmd.c_str());
    data.clear();
    boost::split(data,sam,boost::is_any_of("\n"));
    for (auto item : data)
    {
        if (item == "") {continue;}
        if (item.substr(0,1) == "[") {continue;}
        vector <string> fields;
        boost::split(fields,item,boost::is_any_of("\t"));
        string read = fields[0];
        if (fields[1] == "*") {continue;}
        if (fields[2] == "*") {continue;}
        if (fields.size() >= 11)
        {
            string nm = fields[11];
            boost::replace_all(nm, "NM:i:", "");
            if (reads_to_primary_A.find(read) != reads_to_primary_A.end())
            {
                reads_to_primary_A[read] = reads_to_primary_A[read] + stoi(nm);
            }
            if (reads_to_primary_A.find(read) == reads_to_primary_A.end())
            {
                reads_to_primary_A[read] = stoi(nm);
            }
        }
    }
    
    
    map <string,int> reads_to_primary_B;
    tmp_ref = v_db + "/adjustment/cds/" + gene + "/" + alleles[1] + ".fas";
    log = v_output + "/adjustment/" + gene + "/tmp_" + alleles[1] + ".log";;
    cmd = v_bwa + " mem -B2 -L 1 " + tmp_ref + " " + r1 + " " + " 2> " + log;
    sam = GetStdoutFromCommand(cmd.c_str());
    data.clear();
    boost::split(data,sam,boost::is_any_of("\n"));
    for (auto item : data)
    {
        if (item == "") {continue;}
        if (item.substr(0,1) == "[") {continue;}
        vector <string> fields;
        boost::split(fields,item,boost::is_any_of("\t"));
        string read = fields[0];
        if (fields[1] == "*") {continue;}
        if (fields[2] == "*") {continue;}
        if (fields.size() >= 11)
        {
            string nm = fields[11];
            boost::replace_all(nm, "NM:i:", "");
            
            if (reads_to_primary_B.find(read) != reads_to_primary_B.end())
            {
                reads_to_primary_B[read] = reads_to_primary_B[read] + stoi(nm);
            }
            if (reads_to_primary_B.find(read) == reads_to_primary_B.end())
            {
                reads_to_primary_B[read] = stoi(nm);
            }
        }
    }
    
    tmp_ref = v_db + "/adjustment/cds/" + gene + "/" + alleles[1] + ".fas";
    log = v_output + "/adjustment/" + gene + "/tmp_" + alleles[1] + ".log";;
    cmd = v_bwa + " mem -B2 -L 1 " + tmp_ref + " " + r2 + " " + " 2> " + log;
    sam = GetStdoutFromCommand(cmd.c_str());
    data.clear();
    boost::split(data,sam,boost::is_any_of("\n"));
    for (auto item : data)
    {
        if (item == "") {continue;}
        if (item.substr(0,1) == "[") {continue;}
        vector <string> fields;
        boost::split(fields,item,boost::is_any_of("\t"));
        string read = fields[0];
        if (fields[1] == "*") {continue;}
        if (fields[2] == "*") {continue;}
        if (fields.size() >= 11)
        {
            string nm = fields[11];
            boost::replace_all(nm, "NM:i:", "");
            
            if (reads_to_primary_B.find(read) != reads_to_primary_B.end())
            {
                reads_to_primary_B[read] = reads_to_primary_B[read] + stoi(nm);
            }
            if (reads_to_primary_B.find(read) == reads_to_primary_B.end())
            {
                reads_to_primary_B[read] = stoi(nm);
            }
        }
    }
    
    for (auto item : reads_to_primary_A)
    {
        string read = item.first;
        int nmA = item.second;
        int nmB = 1000;
        if (reads_to_primary_B.find(read) != reads_to_primary_B.end())
        {
            nmB = reads_to_primary_B[read];
        }
        if (nmA < nmB)
        {
            return_res["reads"].append("," + read + ";" + to_string(nmA));
        }
        if (nmA > nmB)
        {
            return_res["reads"].append("," + read + ";" + to_string(nmB));
        }
 //       if (nmA == nmB)
 //       {
 //           return_res["reads"].append("," + read + ";" + to_string(nmB));
 //       }
    }
    for (auto item : reads_to_primary_B)
    {
        string read = item.first;
        int nmB = item.second;
        int nmA = 1000;
        if (reads_to_primary_A.find(read) != reads_to_primary_A.end())
        {
            nmA = reads_to_primary_A[read];
        }
        if (nmA < nmB)
        {
            return_res["reads"].append("," + read + ";" + to_string(nmA));
        }
        if (nmA > nmB)
        {
            return_res["reads"].append("," + read + ";" + to_string(nmB));
        }
        if (nmA == nmB)
        {
            return_res["reads"].append("," + read + ";" + to_string(nmB));
        }
    }
    return return_res;
}



void fasta_exons (string bedfile, string inputvcf, string referencefile, string gene, string outfile, string master_output)
{

    ifstream bed( bedfile );
    vector <string> line_data;
    map <int,string> bed_data;
    for( std::string line; getline( bed, line ); )
    {
        if (line == "") {continue;}
        boost::split(line_data,line,boost::is_any_of("\t"));
        int start = stoi(line_data[1]);
        bed_data[start] = line;
    }
    bed.close();
    
 

    for(auto &&item: bed_data)
    {
        vector <string> line_data;
        boost::split(line_data,item.second,boost::is_any_of("\t"));
        string v_chr = line_data[0];
        string vstart = line_data[1];
        string vend = line_data[2];

        string outseg = master_output + "/" + gene + "." + vstart + ".fas";

        check_fasta(vstart, vend, inputvcf, referencefile, v_chr, outseg);
    }
 
 
    map <string,string> fasta_data;
    for(auto &&item: bed_data)
    {
        
        vector <string> line_data;
        boost::split(line_data,item.second,boost::is_any_of("\t"));
        string v_chr = line_data[0];
        string vstart = line_data[1];
        string vend = line_data[2];

        string file = master_output + "/" + gene + "." + vstart + ".fas";
       
        ifstream fasta( file );
        vector <string> data;
        for( std::string line; getline( fasta, line ); )
        {
            string id = line.substr(1);
            getline( fasta, line );
            string seq = line;
            
            fasta_data[id] = fasta_data[id] + seq;
        }
        fasta.close();
        string res = GetStdoutFromCommand("rm " + file);
    }

   
    ofstream fasta;
    fasta.open (outfile);
    for(auto &&item: fasta_data)
    {
        fasta << ">" << item.first << endl;
        fasta << item.second << endl;
    }
    fasta.close();
    
}




void check_fasta(string vstart, string vend, string inputvcf, string referencefile, string v_chr, string outfile)
{
  
    int $higher = 0;
    int $lower = 0;
    if (((vstart == "0") || (vend == "0")))
    {
        vector<string> line_data;
        ifstream input( inputvcf );
        for( std::string line; getline( input, line ); )
        {
            if (line.substr(0,1) == "#"){continue;}
            if (line == ""){continue;}
            boost::split(line_data,line,boost::is_any_of("\t"));
            if ($higher == 0) {$higher = stoi(line_data[1]);}
            if ($lower == 0) {$lower = stoi(line_data[1]);}
            if (stoi(line_data[1]) >= $higher){$higher = stoi(line_data[1]);}
        }
        input.close();
    }
    if (vstart != "0") {$lower = stoi(vstart);}
    if (vend != "0") {$higher = stoi(vend);}

    

    string sequence = "";
    ifstream input( referencefile );
    for( std::string line; getline( input, line ); )
    {
        if (line.substr(0,1) == ">"){continue;}
        if (line.substr(0,1) == ""){continue;}
        sequence += line;
    }
    input.close();
    

    sequence = sequence.substr($lower-1, ($higher - $lower)+1);
    
    if ((sequence == "") || (sequence.size() < ($higher - $lower)))
    {
        warnings.push_back("It was not possible to retrieve a proper reference sequence");
     //   PrintWarnings();
        return;
    }



    //loading vcf data
    vector<string> head;
    vector<string> subdata;
    vector<string> genotypes;
    vector<string> positions;
    std::map <pair<string,string>, string> snp_data;
    std::map <string,int> samples;
    std::map <pair<string,string>,string> alt_data;
    
    int snps = 0;
    int format_pos = 0;
    int ref_pos = 0;
    int alt_pos = 0;
    
    ifstream vcf( inputvcf );
    for( std::string line; getline( vcf, line ); )
    {
        if (line == "") {continue;}
        if (line.substr(0,2) == "##"){continue;}
        if (line.substr(0,4) == "#CHR")
        {
            boost::split(head,line,boost::is_any_of("\t"));
            int a = 0;
            for(vector<string>::iterator sample = head.begin();sample!=head.end();++sample)
            {
                if (*sample == "REF"){ref_pos = a;}
                if (*sample == "ALT"){alt_pos = a;}
                if (*sample == "FORMAT"){format_pos = a;break;}
                a++;
            }
            continue;
        }
        
        boost::split(subdata,line,boost::is_any_of("\t"));
        if ((subdata[0] != v_chr) && (v_chr != "")){continue;}
        if ((stoi(subdata[1]) < stoi(vstart)) && (vstart != "0")) {continue;}
        if ((stoi(subdata[1]) > stoi(vend) ) && (vend != "0")) {continue;}
        
        positions.push_back (subdata[1]);
        snps++;
        
        vector <string> alternatives;
        boost::split(alternatives,subdata[alt_pos],boost::is_any_of(","));
        
        pair <string,string> key = make_pair(subdata[1],"0");
        alt_data[key] = subdata[ref_pos];
        
        int a = 1;
        for(vector<string>::iterator alt = alternatives.begin();alt!=alternatives.end();++alt)
        {
            pair <string,string> key = make_pair(subdata[1],to_string(a));
            alt_data[key] = *alt;
            a++;
        }
        
        for (a = format_pos + 1; a < subdata.size(); a++)
        {
            boost::split(genotypes,subdata[a],boost::is_any_of(":"));
            pair <string,string> key = make_pair(head[a],subdata[1]);
            snp_data[key] = genotypes[0];
            
        }
    }
    vcf.close();




    if (format_pos == 0) { 
//        warnings.push_back ("It was not possible to detect the FORMAT field. Please check the VCF file."); 
//        PrintWarnings(); return;
    }
    if (snps == 0) {
 //       warnings.push_back ("No variants within this segment and all samples will present the same reference sequence.");
    }



    
    for (int a = format_pos + 1; a < subdata.size(); a++)
    {
        samples[head[a]] = 1;
    }
    
    ofstream myfile;
    myfile.open (outfile);
    
    for (int a = format_pos + 1; a < subdata.size(); a++)
    {
        myfile << ">" << head[a] << "_h1\n";
        string draft = sequence;
        int correct = 0;
        for(vector<string>::iterator pos = positions.begin();pos!=positions.end();++pos)
        {
            pair <string,string> key = make_pair(head[a],*pos);
            string gen = snp_data[key];
            vector<string> alleles;
            boost::split(alleles,gen,boost::is_any_of("|"));
            
            if (alleles.size() == 2)
            {
                if (alleles[0] == "0") {continue;}
                else
                {
                    pair <string,string> key1 = make_pair(*pos,"0");
                    pair <string,string> key2 = make_pair(*pos,alleles[0]);
                    
                    if (alt_data[key2] == "*") {continue;}
                    
                    int b = stoi(*pos)+correct-$lower;
                    
                    draft.replace(b,alt_data[key1].size(),alt_data[key2]);
                    correct = correct + (alt_data[key2].size() - alt_data[key1].size());
                }
            }
        }
        std::transform(draft.begin(), draft.end(),draft.begin(), ::toupper);
        myfile << draft << endl;
        
        
        myfile << ">" << head[a] << "_h2\n";
        draft = sequence;
        correct = 0;
        for(vector<string>::iterator pos = positions.begin();pos!=positions.end();++pos)
        {
            pair <string,string> key = make_pair(head[a],*pos);
            string gen = snp_data[key];
            vector<string> alleles;
            boost::split(alleles,gen,boost::is_any_of("|"));
            
            if (alleles.size() == 2)
            {
                
                if (alleles[1] == "0") {continue;}
                else
                {
                    pair <string,string> key1 = make_pair(*pos,"0");
                    pair <string,string> key2 = make_pair(*pos,alleles[1]);
                    
                    if (alt_data[key2] == "*") {continue;}
                    
                    
                    int b = stoi(*pos)+correct-$lower;
                    
                    draft.replace(b,alt_data[key1].size(),alt_data[key2]);
                    correct = correct + (alt_data[key2].size() - alt_data[key1].size());
                }
            }
        }
        std::transform(draft.begin(), draft.end(),draft.begin(), ::toupper);
        myfile << draft << endl;
    }
    myfile.close();
    
}

double Average(std::vector<double> samples)
{
     double sum;
	 double count;
	 for (auto sub : samples)
	 {
		 sum = sum + sub;
		 count++;
	 }
	 double mean = sum / count;
	 return mean;
}



double Variance(std::vector<double> samples)
{
     int size = samples.size();

     double variance = 0;
     double t = samples[0];
     for (int i = 1; i < size; i++)
     {
          t += samples[i];
          double diff = ((i + 1) * samples[i]) - t;
          variance += (diff * diff) / ((i + 1.0) *i);
     }

     return variance / (size - 1);
}

double StandardDeviation(std::vector<double> samples)
{
     return sqrt(Variance(samples));
}
