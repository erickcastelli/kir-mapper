//  kir-mapper
//
//  Created by Erick C. Castelli
//  2024 GeMBio.Unesp.
//  erick.castelli@unesp.br


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
#include <filesystem>

#include "setup.hpp"
#include "external.hpp"
#include "functions.hpp"
#include "exome_bias.hpp"
#include "ThreadPool.hpp"


namespace fs = filesystem;
using namespace std;

mutex mtx_exomebias;

void main_exome_bias ()
{

    int v_check = 0;
    if (v_output != "") {v_check = 1;}
 
    if (((v_db == "") || (v_check == 0)))
    {
        screen_message (screen_size, 0, "", 1, v_quiet);
        v_message = "Program:   " + Program_name + "::exomebias";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        v_message = "Version:   " + Program_version + ", " + Program_date;
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        v_message = "Usage:     kir-mapper exomebias -output map_output_folder <options>";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Mandatory options:", 1, v_quiet);
        screen_message (screen_size, 2, "-output      output folder (same as map)", 1, v_quiet);
 
        screen_message (screen_size, 2, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Other options:", 1, v_quiet);
        screen_message (screen_size, 2, "-db          path to the kir-mapper database", 1, v_quiet);
        screen_message (screen_size, 2, "-threads     number of threads [" + v_threads + "]", 1, v_quiet);
        screen_message (screen_size, 2, "--quiet     quiet mode", 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        return;
    }

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
        main_exome_bias();
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
        main_exome_bias();
        return;
    }
    
    if (v_output == "") {
        v_message = "You must indicate a valid map output folder";
        warnings.push_back (v_message);
        main_exome_bias();
        return;
    }
    if (! fileExists(v_output))
    {
        v_message = "You must indicate a valid map output folder";
        warnings.push_back (v_message);
        v_output = "";
        main_exome_bias();
        return;
    }
    debug_message("Checking database - end");





    debug_message("Creating output structure - start");
    v_command = "mkdir " + v_output;
    v_system_out = GetStdoutFromCommand(v_command);
    v_command = "mkdir " + v_output + "/exomebias/";
    v_system_out = GetStdoutFromCommand(v_command);
    
     if (! fileExists(v_output))
    {
        v_message = "Error creating the output folder.";
        warnings.push_back (v_message);
        v_output = "";
        main_exome_bias();
        return;
    }
    debug_message("Creating output structure - end");



    screen_message (screen_size, 0, "", 1, v_quiet);
    v_message = "Program:   " + Program_name + "::exomebias";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "Version:   " + Program_version + ", " + Program_date;
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "Cores:     using " + v_threads + " thread(s)";
    screen_message (screen_size, 0, v_message, 1, v_quiet);






    debug_message("Loading bams - start");

    v_message = " > Loading BAM files...";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    map <string,string> filedb;
    string path_map = v_output + "/map/";
    for (const auto & entry : fs::directory_iterator(path_map))
    {
        string dir_path = entry.path();
        string sample = base_name(dir_path);

        string bamunique = dir_path + "/" + sample + ".unique.bam";
        string bamuniquenodup = dir_path + "/" + sample + ".unique.nodup.bam";
        
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
    v_output = v_output + "/exomebias/";

    if (filedb.size() == 0)
    {
        v_message = "The list of BAM files is empty!";
        warnings.push_back (v_message);
        return;
    }
    v_message = " > Processing " + to_string(filedb.size()) + " BAM files.";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    
    debug_message("Loading bams - end");




    string targetfile = "";
    if (v_exome == 1) {targetfile = v_db + "/ncopy/targets.cds.txt";}
    if (v_exome == 0) {targetfile = v_db + "/ncopy/targets.full.txt";}



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





   



    vector <string> list_of_bam_files;
    map <string,string> sample_bam;


    for (auto file : filedb)
    {
        list_of_bam_files.push_back(file.second);
        sample_bam[file.second] = file.first;
    }

    for (auto item : sample_bam)
    {
        string sample = item.second;
        string file = item.first;
        string presencefile = v_output_map + "/map/" + sample + "/" + sample + "_presence_report.txt";
        if (! fileExists(presencefile)) {cout << "Error loading gene presence file" << endl;}
        ifstream in(presencefile.c_str());
        for( std::string line; getline( in, line ); )
        {
            if (line == ""){continue;}
            vector <string> data;
            boost::split(data,line,boost::is_any_of("\t"));
            string gene = data[0];
            pair <string,string> key = make_pair(file,gene);
            genes_avaliable_in_sample[key] = data[2];

        }
        in.close();
            pair <string,string> key = make_pair(file,"HLA-G");
            genes_avaliable_in_sample[key] = "present";
            key = make_pair(file,"HLA-E");
            genes_avaliable_in_sample[key] = "present";
    }




    //v_message = " > Checking depth ...";
    //screen_message (screen_size, 0, v_message, 1, v_quiet);
    //calculate_depth(list_of_bam_files);



    v_message = " > Checking probe capture bias ...";
    screen_message (screen_size, 0, v_message, 1, v_quiet);


    ThreadPool poolbias(stoi(v_threads));
    std::vector< std::future<int> > results_bias;

    int loop = 0;
    for (auto item : targetsdb)
    {
        string gene = item.first;
        string fileout = v_output + gene + ".capture_bias.txt";
        if (fileExists(fileout)) {continue;}

        results_bias.emplace_back(
        poolbias.enqueue([loop,gene, v_output_map, list_of_bam_files,fileout]
            {
                calculate_capture_bias(list_of_bam_files, gene, v_output, fileout, 50);
                return 1;
            })
            );
            loop++;
    }
    for(auto && result: results_bias){result.get();} // waiting for all threads




    string selected_regions = v_output + "/selected_regions.txt";
    ifstream selreg(selected_regions.c_str());
    for( std::string line; getline( selreg, line ); )
    {
        if (line == "") {continue;}
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));
        targetsdb[data[0]] = "," + data[1];
    }
    selreg.close();

   
}



   
string calculate_capture_bias (std::vector<string> files, string gene, string v_out, string fileout, int number_of_tests)
{

    debug_message("Processing " + gene );
    string bed = v_db + "/genotype/bed/" + gene + ".CDS.bed";
    if (! fileExists(bed)) {return "error";}
    float mincov = 5;

    map <int,string> ratios;
    string chr = "";
    int loop = 0;
    for (auto bam : files)
    {
        
        pair <string,string> key = make_pair(bam,gene);
        if (genes_avaliable_in_sample.find(key) == genes_avaliable_in_sample.end()){continue;}
        if (genes_avaliable_in_sample[key] != "present"){continue;}
        
        //float master_cov = exome_depth[bam];
        //if (master_cov < 10) {continue;}

        debug_message("Processing " + bam + " : start depth distribution" );
        string cmd = v_samtools + " depth -a -Q 1 -b " + bed + " --verbosity 0 " + bam;
        string out = GetStdoutFromCommand(cmd);
        vector <string> data;
        boost::split(data,out,boost::is_any_of("\n"));
        if (data.size() == 0) {return "error"; }
        
        float count = 0;
        map <int, float> cov_values;

        float sum = 0;
        float nucs = 0;
        for (auto item : data)
        {
            if (item == "") {continue;}
            vector <string> sub;
            boost::split(sub,item,boost::is_any_of("\t"));
            float cov = stof(sub[2]);
            int pos = stoi(sub[1]);
            chr = sub[0];
            sum = sum + cov;
            nucs++;
            if (cov > mincov) {count++;}
            cov_values[pos] = cov;
        }
        if ((count / float(data.size())) < 0.8) {continue;}
        float master_cov = sum / nucs;


        for (auto item : cov_values)
        {
            int pos = item.first;
            float cov = item.second;
            float ratio = cov / master_cov;
            ratios[pos].append(";" + to_string(ratio));
        }
        debug_message("Processing " + bam + " : end depth distribution" );

        loop++;
//        if (loop == number_of_tests) {break;}
    }

    if (ratios.size() == 0) {return "error";}
 
    
    ofstream bias;
    bias.open (fileout.c_str());
    bias << "GENE\tCHR\tPOS\tMIN\tMAX\tRATIO\tSTD\tVAR\tVALUES" << endl;
    map <int,double> variances;

    for (auto item : ratios)
    {
        int pos = item.first;
        debug_message("calculating statistics for each position: starting pos " + to_string(pos) );
        string covs = item.second;
        vector <string> data;
        boost::split(data,covs,boost::is_any_of(";"));
        vector <double> values;
        for (auto value : data) {if (value == "") {continue;} values.push_back(stod(value));}
        float min = 10000000;
        float max = 0;
        for (auto value : values)
        {
            if (min > value) {min = value;}
            if (max < value) {max = value;}
        }
        double ave = Average(values);
        double std = StandardDeviation(values);
        double var = Variance(values);

        double count = 0;
        double sum = 0;
        string values_line = "";
        for (auto item : values) {
            values_line.append("," + to_string(item));
            sum = sum + item;
            count++;
        }
        double ratio = sum / count;
   


        bias << gene << "\t" << chr << "\t" << pos << "\t" << min << "\t" << max << "\t" << ratio << "\t" << std << "\t" << var << "\t" << values_line.substr(1) << endl;
        variances[pos] = var;
        debug_message("calculating statistics for each position: ending pos " + to_string(pos) );
 
    }
    bias.close();

    int size = 150;
    int best_start = 0;
    int best_end = 0;
    double best_var = 10000;
    vector <int> positions;
    for (auto item : variances)
    {
        positions.push_back(item.first);
    }

    debug_message("searching for most stable segment...start " );
    for (int a = 0; a < (positions.size() - size); a++)
    {
        vector <double> list_of_values;
        int start = 0;
        int end = 0;
        for (int b = 0; b < size; b++)
        {
            int key = a + b;
            int pos = positions[key];
            if (start == 0) {start = pos;}
            if (end == 0) {end = pos;}
            if ((pos - end) > 1){list_of_values.clear(); break;}
            end = pos;
            list_of_values.push_back(variances[pos]);
        }
        if (list_of_values.size() <= 2) {continue;}
        double var = Variance(list_of_values);        
        if (var < best_var)
        {
            best_var = var;
            best_start = start;
            best_end = end;
        }
    }
    debug_message("searching for most stable segment...end ");

    if ((best_start == 0) || (best_end ==  0)) {return "error";}

    mtx_exomebias.lock();
    capture_bias_region[gene] = chr + ":" + to_string(best_start) + "-" + to_string(best_end);
    string selected_regions = v_output + "/capture_bias/selected_regions.txt";
    ofstream selreg(selected_regions, ios::app);
    selreg << gene << "\t" << chr << ":" + to_string(best_start) << "-" << to_string(best_end) << endl;
    selreg.close();
    mtx_exomebias.unlock();

    return "end";
}



string calculate_depth (std::vector<string> files)
{

    debug_message("Processing coverage");
    string bedcov = v_db + "/ncopy/coverage.bed";
    if (! fileExists(bedcov)) {return "error";}

    ThreadPool poolbias(stoi(v_threads));
    std::vector< std::future<int> > results_bias;

    int loop = 0;
    for (auto bam : files)
    {
        
        results_bias.emplace_back(
            poolbias.enqueue([loop, bam, bedcov]
                {

        debug_message("Processing " + bam + " : start main cov" );
        string cmd = v_samtools + " depth -a -Q 1 -b " + bedcov + " --verbosity 0 " + bam;
        string out = GetStdoutFromCommand(cmd);
        vector <string> data;
        boost::split(data,out,boost::is_any_of("\n"));
        if (data.size() == 0) {return 1; }
        float sum = 0;
        float div = 0;
        for (auto item : data)
        {
            if (item == "") {continue;}
            vector <string> sub;
            boost::split(sub,item,boost::is_any_of("\t"));
            float cov = stof(sub[2]);
            sum = sum + cov;
            div++;
        }
        float master_cov = sum / div;
        mtx_exomebias.lock();
        exome_depth[bam] = master_cov;
        debug_message("Processing " + bam + " : end main cov" );
        mtx_exomebias.unlock();
        return 1;
        })
                );
                loop++;
    }
    for(auto && result: results_bias){result.get();} // waiting for all threads
    return "end";
}