//  kir-mapper
//
//  Created by Erick C. Castelli
//  2024 GeMBio.Unesp.
//  erick.castelli@unesp.br
// Contributions from code on tht web

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <map>
#include <boost/algorithm/string.hpp>
#include <string>
#include <thread>
#include <unordered_map>
#include <cstring>
#include <cassert>
#include <future>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <mutex>
#include <filesystem>

#include <numeric>
#include <algorithm>
#include <assert.h>

#include "external.hpp"
#include "functions.hpp"
#include "ThreadPool.hpp"
#include "ncopy.hpp"

namespace fs = filesystem;
using namespace std;

mutex mtx_ncopy;



// https://stackoverflow.com/a/12399290/7128154
template <typename t>
std::vector<size_t> sorted_index(const std::vector<t> &v) {

  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}
// https://stackoverflow.com/a/1267878/7128154
template< typename order_iterator, typename value_iterator >
void reorder( order_iterator order_begin, order_iterator order_end, value_iterator v )  {
    typedef typename std::iterator_traits< value_iterator >::value_type value_t;
    typedef typename std::iterator_traits< order_iterator >::value_type index_t;
    typedef typename std::iterator_traits< order_iterator >::difference_type diff_t;

    diff_t remaining = order_end - 1 - order_begin;
    for ( index_t s = index_t(), d; remaining > 0; ++ s ) {
        for ( d = order_begin[s]; d > s; d = order_begin[d] ) ;
        if ( d == s ) {
            -- remaining;
            value_t temp = v[s];
            while ( d = order_begin[d], d != s ) {
                swap( temp, v[d] );
                -- remaining;
            }
            v[s] = temp;
        }
    }
}

// https://stackoverflow.com/a/1267878/7128154
template< typename order_iterator, typename value_iterator >
void reorder_destructive( order_iterator order_begin, order_iterator order_end, value_iterator v )  {
    typedef typename std::iterator_traits< value_iterator >::value_type value_t;
    typedef typename std::iterator_traits< order_iterator >::value_type index_t;
    typedef typename std::iterator_traits< order_iterator >::difference_type diff_t;

    diff_t remaining = order_end - 1 - order_begin;
    for ( index_t s = index_t(); remaining > 0; ++ s ) {
        index_t d = order_begin[s];
        if ( d == (diff_t) -1 ) continue;
        -- remaining;
        value_t temp = v[s];
        for ( index_t d2; d != s; d = d2 ) {
            std::swap( temp, v[d] );
            std::swap( order_begin[d], d2 = (diff_t) -1 );
            -- remaining;
        }
        v[s] = temp;
    }
}



// https://stackoverflow.com/a/29677616/7128154
// https://stackoverflow.com/a/37708864/7128154
template <typename t>
double quantile(double q, std::vector<t> values, std::vector<double> weights = std::vector<double>())
{
    assert( 0. <= q && q <= 1. && "expecting quantile in range [0; 1]");
    if (weights.empty())
    {
        weights = std::vector<double>(values.size(), 1.);
    }
    else
    {
        assert (values.size() == weights.size()  && "values and weights missfit in quantiles");
        std::vector<size_t> inds = sorted_index(values);
        reorder_destructive(inds.begin(), inds.end(), weights.begin());
    }

    stable_sort(values.begin(), values.end());
    // values and weights are sorted now

    std::vector<double> quantiles (weights.size());
    quantiles[0] = weights[0];
    for (int ii = 1; ii < quantiles.size(); ii++)
    {
        quantiles[ii] = quantiles[ii-1] + weights[ii];
    }
    double norm = std::accumulate(weights.begin(), weights.end(), 0.0);
    int ind = 0;
    double qcurrent = 0;
    for (; ind < quantiles.size(); ind++)
    {
        qcurrent = (quantiles[ind] - weights[ind] / 2. ) / norm;
        quantiles[ind] = qcurrent;
        if (qcurrent > q)
        {
            if (ind == 0) {return values[0];}
            double rat = (q - quantiles[ind-1]) / (quantiles[ind] - quantiles[ind-1]);
            return values[ind-1] + (values[ind] - values[ind-1]) * rat;
        }
    }
    return values[values.size()-1];

}

template <typename t>
double quantile(double q, std::vector<t> values, std::vector<int> weights)
{
    std::vector<double> weights_double (weights.begin(), weights.end());
    return quantile(q, values, weights_double);
}


void main_ncopy ()
{



/*
    std::vector<int> vals {96,35,59,52,14,51,28,8,13,10,75,83,20};
    std::cout << "quantile(0.25, vals)=" << quantile(0.25, vals) << std::endl;
    std::cout << "quantile(.75, vals)=" << quantile(.75, vals) << std::endl;

    std::vector<int> vals2 {43,36,55,34,15,89,68,46,11,12,67,73,15};
    std::cout << "quantile(0.25, vals)=" << quantile(0.25, vals2) << std::endl;
    std::cout << "quantile(.75, vals)=" << quantile(.75, vals2) << std::endl;


    return;

    std::vector<int> vals2 {1, 2, 3};
    std::vector<double> ws2 {1, 2, 3};
    std::cout << "quantile(.13, vals2, ws2)=" << quantile(.13, vals2, ws2) << std::endl;
*/
    
    
    int v_check = 0;
    if (v_output != "") {v_check = 1;}
 
    if (((v_db == "") || (v_check == 0)))
    {
        screen_message (screen_size, 0, "", 1, v_quiet);
        v_message = "Program:   " + Program_name + "::ncopy";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        v_message = "Version:   " + Program_version + ", " + Program_date;
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        v_message = "Usage:     kir-mapper ncopy -output map_output_folder <options>";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Mandatory options:", 1, v_quiet);
        screen_message (screen_size, 2, "-output      output folder (same as map and ncopy)", 1, v_quiet);
 
        screen_message (screen_size, 2, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Other options:", 1, v_quiet);
        screen_message (screen_size, 2, "-db          path to the kir-mapper database", 1, v_quiet);
        screen_message (screen_size, 2, "-threads     number of threads [" + v_threads + "]", 1, v_quiet);
        screen_message (screen_size, 2, "-reference   KIR3DL3,5UPKIR,HLA-E,HLA-G [default: KIR3DL3]", 1, v_quiet);
        screen_message (screen_size, 2, "-samples     text file listing the samples to consider", 1, v_quiet);
        screen_message (screen_size, 2, "--exome     only exons", 1, v_quiet);
        screen_message (screen_size, 2, "--quiet     quiet mode", 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        return;
    }
    
    float mincov = 15;
    float mincovref = 5;
    float mindepthratio = 0.45f;
    int number_of_data_points = 5;
    


    debug_message("Checking database - start");
    // checking database
    string v_db_info = v_db + "/db_dna.info";
    boost::replace_all(v_db_info, "\\ ", " ");

    if (! fileExists(v_db_info))
    {
        v_message = "Error accessing database " + v_db_info;
        cout << v_message << endl;
        warnings.push_back (v_message);
        v_sample = "";
        v_r1 = "";
        main_ncopy();
        return;
    }


    ifstream file_db(v_db_info.c_str());
    int db_version_ok = 0;
    for( std::string line; getline( file_db, line ); )
    {
        if (line == "kir-mapper:1") {
            db_version_ok = 1;
            continue;
        }
        vector<string> db_data;
        boost::split(db_data,line,boost::is_any_of(":"));
        if (db_data[0] == "GENE")
        {
            v_genes = v_genes + db_data[1] + ",";
            chr_hg38[db_data[1]] = db_data[2];
            position_hg38[db_data[1]] = stoi(db_data[3]);
            chr_size_hg38[db_data[1]] = db_data[4];
            gene_opt_start[db_data[1]] = stoi(db_data[5]);
            gene_opt_end[db_data[1]] = stoi(db_data[6]);
            gene_type[db_data[1]] = db_data[7];
        }


        if (db_data[0] == "SIZE")
        {
            v_size = stoi(db_data[1]);
        }
    }
    file_db.close();
    
    
    if (db_version_ok == 0) {
        v_message = "The database is not compatible with this version of kir-mapper.";
        warnings.push_back (v_message);
        v_sample = "";
        v_r1 = "";
        main_ncopy();
        return;
    }
    
    if (v_output == "") {
        v_message = "You must indicate a valid map output folder";
        warnings.push_back (v_message);
        main_ncopy();
        return;
    }
    if ((! fileExists(v_output)) || (! fileExists(v_output + "/map")))
    {
        v_message = "You must indicate a valid folder, with 'map'";
        warnings.push_back (v_message);
        v_output = "";
        main_ncopy();
        return;
    }
	
	
    debug_message("Checking database - end");


    debug_message("Checking thresholds - start");

    string thresholds = v_db + "/ncopy/thresholds.txt";;
    if (! fileExists(thresholds)) {
        v_message = "There is something wrong with this database!";
        warnings.push_back (v_message);
        main_ncopy();
        return;
    }
    debug_message("Checking thresholds - end");
   
   
   
   
    
    debug_message("Creating output structure - start");

    v_command = "mkdir " + v_output;
    v_system_out = GetStdoutFromCommand(v_command);
    v_command = "mkdir " + v_output + "/ncopy/";
    v_system_out = GetStdoutFromCommand(v_command);
    v_command = "mkdir " + v_output + "/ncopy/plots/";
    v_system_out = GetStdoutFromCommand(v_command);
  
    if (! fileExists(v_output))
    {
        v_message = "Error creating the output folder.";
        warnings.push_back (v_message);
        v_output = "";
        main_ncopy();
        return;
    }
    debug_message("Creating output structure - end");

        
    debug_message("Copying thresholds when it doesn't exist - start");
    string test_thresholds = v_output + "/ncopy/thresholds.txt";
    if (fileExists(test_thresholds)) {thresholds = test_thresholds;}
    if (! fileExists(test_thresholds))
    {
        v_command = "cp " + thresholds + " " + test_thresholds;
        v_system_out = GetStdoutFromCommand(v_command);
        thresholds = test_thresholds;
    }
   debug_message("Copying thresholds when it doesn't exist - end");
 

    debug_message("Copying R app");
    string R_app_source = v_db + "/ncopy/kir-mapper_plot_app.R";
    string R_app_dest = v_output + "/ncopy/kir-mapper_plot_app.R";

    if (fileExists(R_app_source)) {
        v_command = "cp " + R_app_source + " " + R_app_dest;
        v_system_out = GetStdoutFromCommand(v_command);
    }
    debug_message("Copying R app - end");


    
   map <string,int> samples_to_include;
    if (v_sample_list != "")
    {
         debug_message("loading sample file - start");
         ifstream samplefile(v_sample_list.c_str());
         for( std::string line; getline( samplefile, line ); )
         {
            samples_to_include[line] = 1;
         }
         samplefile.close();
         debug_message("loading sample file - end");
    }




    
        

 
    
    
    
    
    
    screen_message (screen_size, 0, "", 1, v_quiet);
    v_message = "Program:   " + Program_name + "::ncopy";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "Version:   " + Program_version + ", " + Program_date;
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "Cores:     using " + v_threads + " thread(s)";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    if (v_exome == 1) {
        v_message = "Mode:      exome - only exons";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
    }
    if (v_exome == 0) {
        v_message = "Mode:      full genome - include introns";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
    }

    





    debug_message("Loading bams - start");

    v_message = " > Loading BAM files...";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    map <string,string> filedb;
    string path_map = v_output + "/map/";
    for (const auto & entry : fs::directory_iterator(path_map))
    {
        string dir_path = entry.path();
        string sample = base_name(dir_path);

        if (v_sample_list != ""){
            if (samples_to_include.find(sample) == samples_to_include.end()) {continue;}
        }

        string bamunique = "";
        string bamuniquenodup = "";

        if (v_exome == 1) {
            bamunique = dir_path + "/" + sample + ".unique.bam";
            bamuniquenodup = dir_path + "/" + sample + ".unique.nodup.bam";
        }
        if (v_exome == 0) {
            bamunique = dir_path + "/" + sample + ".unique.bam";
            bamuniquenodup = dir_path + "/" + sample + ".unique.nodup.bam";
        }

 

        if (fileExists(bamuniquenodup))
        {
            string cmd = v_samtools + " view -H " + bamuniquenodup;
            string head = GetStdoutFromCommand(cmd.c_str());
            int validbam = 0;
            if (head.find("kir-mapper") != std::string::npos) {validbam = 1;}
            if (head.find("hla-mapper") != std::string::npos) {validbam = 1;}
            if (validbam == 1) {filedb[sample] = bamuniquenodup;}
        }
        if (! fileExists(bamuniquenodup))
        {
            if (! fileExists(bamunique)) {continue;}
            string cmd = v_samtools + " view -H " + bamunique;
            string head = GetStdoutFromCommand(cmd.c_str());
            int validbam = 0;
            if (head.find("kir-mapper") != std::string::npos) {validbam = 1;}
            if (head.find("hla-mapper") != std::string::npos) {validbam = 1;}
            if (validbam == 1) {filedb[sample] = bamunique;}
        }
    }
    
    string v_output_map = v_output;
    v_output = v_output + "/ncopy/";
 
    

    
    
    
    if (filedb.size() == 0)
    {
        v_message = "The list of BAM files is empty!";
        warnings.push_back (v_message);
        return;
    }
    v_message = " > Processing " + to_string(filedb.size()) + " BAM files.";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    
    
    string targetfile = "";
    if (v_exome == 1) {targetfile = v_db + "/ncopy/targets.cds.txt";}
    if (v_exome == 0) {targetfile = v_db + "/ncopy/targets.full.txt";}

    if (! fileExists(targetfile)) {
        v_message = "There is something wrong with this database!";
        warnings.push_back (v_message);
        main_ncopy();
        return;
    }
    debug_message("Loading bams - end");

   

    debug_message("Loading targets and references - start");

    map <string,string> targetsdb;
    string reference = "";
    ifstream targets(targetfile.c_str());
    for( std::string line; getline( targets, line ); )
    {
        if (line == "") {continue;}
        if (line.substr(0,1) == "#") {continue;}
        vector <string> data;
        boost::split(data,line,boost::is_any_of(":"));
        if (data[0] == "TARGET")
        {
            targetsdb[data[1]].append("," + data[2] + ":" + data[3] + "-" + data[4]);
        }
        if (data[0] == "REF")
        {
            reference = data[1];
        }
    }
    if (v_reference_name != "") {reference = v_reference_name;}
    targets.close();
    v_message = " > Evaluating " + to_string(targetsdb.size()) + " genes.";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = " > Using " + reference + " as reference.";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    debug_message("Loading targets and references - end");

    
    
    


    map <string,float> thresh_1;
    map <string,float> thresh_2;
    map <string,float> thresh_3;
    map <string,float> thresh_4;
    
    ifstream thre(thresholds.c_str());
    for( std::string line; getline( thre, line ); )
    {
        if (line == "") {continue;}
        vector <string> data;
        boost::split(data,line,boost::is_any_of(":"));
        thresh_1[data[0]] = stof(data[1]);
        thresh_2[data[0]] = stof(data[2]);
        thresh_3[data[0]] = stof(data[3]);
        thresh_4[data[0]] = stof(data[4]);
    }
    thre.close();
    
    
    map <pair<string,string>,float> copyresults;
//    map <pair<string,string>,string> hetresults;
    map <pair<string,string>,int> hetcounts;



    debug_message("reindexing bams - start");
    v_message = " > Re-indexing bams if necessary ...";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    for (auto file : filedb)
    {
        string bai = file.second + ".bai";
        if (! fileExists(bai))
        {
            string cmd = v_samtools + " index -@ " + v_threads + " " + file.second;
            GetStdoutFromCommand(cmd);
            continue;
        }
        vector <string> list;
        list.push_back(file.second);
        if (! isNewer(bai, list))
        {
            string cmd = v_samtools + " index -@ " + v_threads + " " + file.second;
            GetStdoutFromCommand(cmd);
        }
    }
    debug_message("reindexing bams - end");
    









    map <pair<string,string>,string> depth_values_list;

  if ((v_exome == 0) || (v_exome == 1))
  {
        debug_message("Evaluating depth - start");
        for (auto item : targetsdb)
        {
            string gene = item.first;
            
            v_message = " > Processing " + gene + " ...";
            screen_message (screen_size, 0, v_message, 1, v_quiet);

            vector <string> regions;
            string subtext = targetsdb[gene].substr(1);
            boost::split(regions,subtext,boost::is_any_of(","));

            
            ThreadPool pool(stoi(v_threads));
            std::vector< std::future<int> > results;
            int loop = 1;
            
            for (auto file : filedb)
            {
                
                results.emplace_back(
                pool.enqueue([loop,file,gene,regions,mincov,mindepthratio,&hetcounts,&copyresults,&depth_values_list]
                {
                
                    if (file.second == "") {return 1;}
                    string sample = file.first;
                    float sum = 0;
                    int hetcount = 0;
                    int hetvalid = 0;

                    for (auto region : regions)
                    {
                        if (region == "") {continue;}
                        string v_cmd = v_samtools + " coverage --ff SECONDARY,DUP -q 1 -r " + region + " " + file.second;
                        string depth = GetStdoutFromCommand(v_cmd);
                        
                        if (depth == "") {
                            depth = GetStdoutFromCommand(v_cmd);
                        }
                        if (depth == "") {
                            depth = GetStdoutFromCommand(v_cmd);
                        }

                        vector <string> data;
                        boost::split(data,depth,boost::is_any_of("\n"));
                        if (data.size() < 1) {return 1;}
                        
                        int index = 0;
                        for (auto item : data)
                        {
                            if (item.find("startpos") == std::string::npos)
                            {
                                break;
                            }
                            index++;
                        }
                        
                        vector <string> fields;
                        boost::split(fields,data[index],boost::is_any_of("\t"));
    
                        if (fields.size() >= 6) {
                            sum = sum + stof(fields[6]);

                            pair <string,string> key = make_pair(gene,sample);
                            mtx_ncopy.lock();
                            depth_values_list[key].append(";" + fields[6]);
                            mtx_ncopy.unlock();

                            if (stof(fields[6]) < mincov) {continue;}
                        }

                        
                        v_cmd = v_samtools + " mpileup -r " + region + " --ff SECONDARY,DUP -q 1 " + file.second;
                        string pile = GetStdoutFromCommand(v_cmd);
                        if (pile == "")
                        {
                            pile = GetStdoutFromCommand(v_cmd);
                        }
                        if (pile == "")
                        {
                            pile = GetStdoutFromCommand(v_cmd);
                        }
                        
                        
                        boost::split(data,pile,boost::is_any_of("\n"));
                        for (auto pileline : data)
                        {
                            if (pileline == "") {continue;}
                            if (pileline.substr(0,1) == "[") {continue;}
                            vector <string> sub;
                            boost::split(sub,pileline,boost::is_any_of("\t"));
                            if (sub.size() > 3)
                            {
                                string cov = sub[3];
                                if ((cov == "") || (cov == "0")) {continue;}
                                
                                if (stof(cov) < mincov) {continue;}
                                hetvalid++;
                                std::transform(sub[4].begin(), sub[4].end(),sub[4].begin(), ::toupper);
                                string bases = sub[4];
                                if (bases.find("+") != std::string::npos) {continue;}
                                if (bases.find("-") != std::string::npos) {continue;}
                                boost::replace_all(bases , "^", "");
                                boost::replace_all(bases , "[", "");
                                boost::replace_all(bases , "]", "");
                                boost::replace_all(bases , "$", "");
                                boost::replace_all(bases , "A", ",A");
                                boost::replace_all(bases , "G", ",G");
                                boost::replace_all(bases , "C", ",C");
                                boost::replace_all(bases , "T", ",T");
                                
                                vector <string> nucs;
                                string subtext = bases.substr(1);
                                boost::split(nucs,subtext,boost::is_any_of(","));
                                map <string,float> freq;
                                for (auto nuc : nucs)
                                {
                                    freq[nuc]++;
                                }
                                
    
                                int valid = 0;
                                if ((freq["A"] / stof(cov)) >= mindepthratio) {valid++;}
                                if ((freq["T"] / stof(cov)) >= mindepthratio) {valid++;}
                                if ((freq["C"] / stof(cov)) >= mindepthratio) {valid++;}
                                if ((freq["G"] / stof(cov)) >= mindepthratio) {valid++;}
    
                                
                                if (valid >= 2)
                                {
                                    hetcount++;
                                }
                            }
                        }
                    }

                    if (hetvalid > 0) {
                        pair <string,string> key = make_pair(gene,sample);
                        mtx_ncopy.lock();
                        hetcounts[key] = hetcount;
                        mtx_ncopy.unlock();
                    }
                    
                    float depth = sum / float(regions.size());
                    pair <string,string> key = make_pair(gene,sample);


                    if (v_exome == 1) {
                        string list = depth_values_list[key];
                        vector <string> loop;
                        vector <float> values;
                        boost::split(loop,list,boost::is_any_of(";"));
                        for (auto item : loop)
                        {
                            if (item == "") {continue;}
                            values.push_back(stof(item));
                        }


                        depth = median(values);

                    }




                    mtx_ncopy.lock();
                    copyresults[key] = depth;
                    mtx_ncopy.unlock();
                    return 1;
                })
                );
                loop++;
            }
            for(auto && result: results){result.get();} // waiting for all threads
        }
        debug_message("Evaluating depth - end");

  }







    v_message = " > Writing tables and VCF ...";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    
    map <string,int> failed_samples;
    for (auto file : filedb)
    {
        if (file.second == "") {continue;}
        string sample = file.first;
        pair <string,string> key = make_pair(reference,sample);
        float ref = copyresults[key];
        if (ref < mincovref) {
            failed_samples[sample] = 1; 
                v_message = " > warning: sample " + sample + " failed due to low depth at the reference (" + to_string(ref) + ")";
                screen_message (screen_size, 0, v_message, 1, v_quiet);
        }
    }
    
    
    string v_out_depth = v_output + "/depth_values.txt";
    string v_out_ratios = v_output + "/ratio_values.txt";
    string v_out_cp = v_output + "/copy_numbers.txt";

    
    ofstream DEPTH;
    DEPTH.open (v_out_depth.c_str());
    ofstream RATIOS;
    RATIOS.open (v_out_ratios.c_str());
    ofstream CNS;
    CNS.open (v_out_cp.c_str());

    
    DEPTH << "SAMPLE";
    RATIOS << "SAMPLE";
    CNS << "SAMPLE\tGENE\tCOPY_NUMBER" << endl;
    
    for (auto item : targetsdb)
    {
        string gene = item.first;
        DEPTH << "\t" << gene;
        RATIOS << "\t" << gene;
    }
    DEPTH << endl;
    RATIOS << endl;
    


    map <string,float> sample_size;
    map <pair<string,int>,float> n_thresh;
    map <pair<string,string>,string> cn_data;

    
    for (auto file : filedb)
    {
        if (file.second == "") {continue;}
        string sample = file.first;
        
        pair <string,string> key = make_pair(reference,sample);
        float ref = 0;
        if (v_exome == 0) {ref = copyresults[key];}
        
        if (v_exome == 1) {
            string list = depth_values_list[key];
            vector <string> loop;
            vector <float> values;
            boost::split(loop,list,boost::is_any_of(";"));
            for (auto item : loop)
            {
                if (item == "") {continue;}
                values.push_back(stof(item));
            }

            ref = median(values);

        }
 
        
        
        DEPTH << sample;
        RATIOS << sample;
 
 
        for (auto item : targetsdb)
        {
            string gene = item.first;
            pair <string,string> key = make_pair(gene,sample);
            float ratio_depth = 0;
            float ratio_motif = 0;
            float ratio = 0;
            if (failed_samples.find(sample) != failed_samples.end())
            {
                DEPTH << "\tNA";
                RATIOS << "\tNA";
            }
            
            if (failed_samples.find(sample) == failed_samples.end())
            {
                
                float depthvalue = 0;
                if (v_exome == 0) {depthvalue = copyresults[key];}
                
                if (v_exome == 1) {
                    string list = depth_values_list[key];
                    vector <string> loop;
                    vector <float> values;
                    boost::split(loop,list,boost::is_any_of(";"));
                    for (auto item : loop)
                    {
                        if (item == "") {continue;}
                        values.push_back(stof(item));
                    }

                    depthvalue = median(values);

                }
                

                DEPTH << "\t" << depthvalue;
                ratio_depth = depthvalue / ref;
                ratio = ratio_depth;
                RATIOS << "\t" << ratio;
                sample_size[gene]++;
            }

            int cn = 0;
            if (ratio >= thresh_1[gene]) {cn = 1;}
            if (ratio >= thresh_2[gene]) {cn = 2;}
            if (ratio >= thresh_3[gene]) {cn = 3;}
            if (ratio >= thresh_4[gene]) {cn = 4;}
            
            if (failed_samples.find(sample) != failed_samples.end())
            {
                CNS << sample << "\t" << gene << "\tNA" << endl;
                pair <string,string> key = make_pair(sample,gene);
                cn_data[key] = "NA";
            }
            else {
                CNS << sample << "\t" << gene << "\t" << cn << endl;
                pair <string,int> key = make_pair(gene,cn);
                n_thresh[key]++;
                pair <string,string> key2 = make_pair(sample,gene);
                cn_data[key2] = to_string(cn);

            }
        }
        DEPTH << endl;
        RATIOS << endl;
    }
    DEPTH.close();
    RATIOS.close();
    CNS.close();
    
    
 

    
    
    string v_out_cp_table = v_output + "/copy_numbers.table.txt";
    ofstream TB;
    TB.open (v_out_cp_table.c_str());
    TB << "SAMPLE";
    for (auto item : targetsdb)
    {
        string gene = item.first;
        TB << "\t" + gene;
    }
    TB << endl;
    
    for (auto file : filedb)
    {
        if (file.second == "") {continue;}
        string sample = file.first;
        
        TB << sample;
        for (auto item : targetsdb)
        {
            string gene = item.first;
            pair <string,string> key = make_pair(sample,gene);
            TB << "\t" << cn_data[key];
        }
        TB << endl;
    }
    TB.close();
    

    
    string v_out_cp_pres = v_output + "/presence.table.txt";
    TB.open (v_out_cp_pres.c_str());
    TB << "SAMPLE";
    for (auto item : targetsdb)
    {
        string gene = item.first;
        TB << "\t" + gene;
    }
    TB << endl;
    
    for (auto file : filedb)
    {
        if (file.second == "") {continue;}
        string sample = file.first;
        
        TB << sample;
        for (auto item : targetsdb)
        {
            string gene = item.first;
            pair <string,string> key = make_pair(sample,gene);
            if (cn_data[key] == "NA") {TB << "\t" << cn_data[key];}
            if (cn_data[key] != "NA") {
                int value = stoi(cn_data[key]);
                if (value == 0) {TB << "\t0";}
                if (value > 0) {TB << "\t1";}
            }
       }
        TB << endl;
    }
    TB.close();
    

    
    vector <string> head;
    vector <string> gene_data;
    string chrom = "";
    
    string model = v_db + "/ncopy/vcf_model.txt";;
    ifstream vcf(model.c_str());
    for( std::string line; getline( vcf, line ); )
    {
        if (line == "") {continue;}
        if (line.substr(0,2) == "##") {head.push_back(line);continue;}
        if (line.substr(0,2) == "#C") {chrom = line;continue;}
        gene_data.push_back(line);
    }
    vcf.close();
    
    string v_out_vcf = v_output + "/presence.table.vcf";
    TB.open (v_out_vcf.c_str());
    for (auto item : head) {TB << item << endl;}
    
    TB << chrom;
    for (auto file : filedb)
    {
        if (file.second == "") {continue;}
        string sample = file.first;
        TB << "\t" << sample;
    }
    TB << endl;
    
    for (auto item : gene_data)
    {
        vector <string> data;
        boost::split(data,item,boost::is_any_of("\t"));
        string gene = data[2];
        
        TB << item;
        for (auto file : filedb)
        {
            if (file.second == "") {continue;}
            string sample = file.first;
            pair <string,string> key = make_pair(sample,gene);
            if (cn_data[key] == "NA") {TB << "\t.";}
            if (cn_data[key] != "NA") {
                int value = stoi(cn_data[key]);
                if (value == 0) {TB << "\t1";}
                if (value > 0) {TB << "\t0";}
            }
        }
        TB << endl;
    }
    TB.close();
    
    

    
 
    debug_message("Building plots - start");

    v_message = " > Building plots ...";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    
    string plotR = v_db + "/ncopy/plot_1gene.R";
    if (! fileExists(plotR)) {
        v_message = "There is something wrong with this kir-mapper database!";
        warnings.push_back (v_message);
        main_ncopy();
        return;
    }
    string script_one;
    ifstream plot(plotR.c_str());
    for( std::string line; getline( plot, line ); )
    {
        script_one.append("\n" + line);
    }
    plot.close();

   string plotR2 = v_db + "/ncopy/plot_2gene.R";
    if (! fileExists(plotR2)) {
        v_message = "There is something wrong with this kir-mapper database!";
        warnings.push_back (v_message);
        main_ncopy();
        return;
    }
    string script_two;
    ifstream plot2(plotR2.c_str());
    for( std::string line; getline( plot2, line ); )
    {
        script_two.append("\n" + line);
    }
    plot2.close();


    for (auto item : targetsdb)
    {
        string gene = item.first;
        
        string secondary = "NA";
        
        if (gene == "KIR2DL2") {secondary = "KIR2DL3";}
        if (gene == "KIR2DL3") {secondary = "KIR2DL2";}

        if (gene == "KIR2DS2") {secondary = "KIR2DL2";}

        if (gene == "KIR3DL1") {secondary = "KIR3DS1";}
        if (gene == "KIR3DS1") {secondary = "KIR3DL1";}
        
        if (gene == "KIR2DS1") {secondary = "KIR2DS4";}
        if (gene == "KIR2DS4") {secondary = "KIR2DS1";}
        
        if (gene == "KIR2DS3") {secondary = "KIR2DS5";}
        if (gene == "KIR2DS5") {secondary = "KIR2DS3";}

        if (gene == "KIR2DL1") {secondary = "KIR2DP1";}
        if (gene == "KIR2DP1") {secondary = "KIR2DL1";}

        

        string v_plotdb = v_output + "/plots/" + gene + "_plot_db.txt";
        ofstream DB;
        DB.open (v_plotdb.c_str());
        DB << "Sample\tGene\tRatio\tHeterozygosis\tOrder" << endl;
        
        map <float,string> ordered_results;
        for (auto file : filedb)
        {
            if (file.second == "") {continue;}
            string sample = file.first;

            if (failed_samples.find(sample) != failed_samples.end()) {continue;}
            
            pair <string,string> key = make_pair(gene,sample);
            pair <string,string> key2 = make_pair(reference,sample);
 
            float depth_ref = copyresults[key2];
            float depth_ratio = copyresults[key] / depth_ref;

 //           float motif_ref = motif_count[key2];
 //           float motif_ratio = motif_count[key] / motif_ref;

            float ratio = depth_ratio;
            
            
            
            string hetdata = "";
            if (hetcounts.find(key) == hetcounts.end()) {hetdata = "heterozygosis: unknown";}
            if (hetcounts.find(key) != hetcounts.end())
            {
                if (hetcounts[key] == 0) {hetdata = "heterozygosis: unknown";}
                if ((hetcounts[key] <= 2) && (hetcounts[key] > 0)) {hetdata = "heterozygosis: yes (1-2)";}
                if ((hetcounts[key] > 2) && (hetcounts[key] > 0)) {hetdata = "heterozygosis: yes (>2)";}
            }
            
            ordered_results[ratio].append(";" + sample + "\t" + gene + "\t" + to_string(ratio) + "\t" + hetdata);
//            ordered_results[ratio].append(";" + sample + "\t" + gene + "\t" + to_string(ratio) + "\t" + hetdata);
        }
        int count = 1;
        for (auto iter : ordered_results)
        {
            string line = iter.second;
            vector <string> lines;
            string subtext = line.substr(1);
            boost::split(lines,subtext,boost::is_any_of(";"));
            for (auto sub : lines)
            {
                DB << sub << "\t" << count << endl;
                
                
                if (secondary != "NA") {
                    vector <string> data;
                    boost::split(data,sub,boost::is_any_of("\t"));
                    pair <string,string> key = make_pair(secondary,data[0]);
                    pair <string,string> key2 = make_pair(reference,data[0]);
                    float ref = copyresults[key2];
                    float ratio = copyresults[key] / ref;
                    DB << data[0] << "\t" + secondary << "\t" + to_string(ratio) << "\theterozygosis: unknown\t" << count << endl;
                }
                count++;
            }
        }
        DB.close();
        
        string scriptgene = "";
        if (secondary == "NA") {scriptgene = script_one;}
        if (secondary != "NA") {scriptgene = script_two;}

        string obs = "";
        if (v_exome == 1) {obs = "Calculated using the exome mode";}
        if (v_exome == 0) {obs = "Calculated using the wgs mode";}
        
        boost::replace_all(scriptgene, "#GENE", gene);
        boost::replace_all(scriptgene, "#SEC", secondary);
        boost::replace_all(scriptgene, "#PATH", v_output + "/plots/");
        boost::replace_all(scriptgene, "#FILE", gene + "_plot_db.txt");
        boost::replace_all(scriptgene, "#THR1", to_string(thresh_1[gene]));
        boost::replace_all(scriptgene, "#THR2", to_string(thresh_2[gene]));
        boost::replace_all(scriptgene, "#THR3", to_string(thresh_3[gene]));
        boost::replace_all(scriptgene, "#THR4", to_string(thresh_4[gene]));
        boost::replace_all(scriptgene, "#REF", reference);
        boost::replace_all(scriptgene, "#CAP",    obs);
        
        string plotfile = v_output + "/plots/" + gene + ".R";
        ofstream PLOTRout;
        PLOTRout.open (plotfile.c_str());
        PLOTRout << scriptgene;
        PLOTRout.close();
        string cmd = "Rscript " + plotfile;
        v_system_out = GetStdoutFromCommand(cmd);
        
    }
    debug_message("Building plots - end");

    
    v_message = " > Computation done.";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "You should explore the plots to evaluate";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "   if you need to change thresholds.";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
     
     
}
    
   