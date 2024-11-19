//  kir-mapper
//
//  Created by Erick C. Castelli
//  2024 GeMBio.Unesp.
//  erick.castelli@unesp.br
// Contributions from code on the web


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

#include "map_dna.hpp"
#include "external.hpp"
#include "functions.hpp"
#include "ThreadPool.hpp"
#include "preselect.hpp"

mutex mtx_hg38_select;
mutex mtx_selection_general;
mutex mtxg;

vector <string> selected_results_r1;
vector <string> selected_results_r2;
vector <string> selected_results_trim_r1;
vector <string> selected_results_trim_r2;
//unordered_map <string,int> motifs_master;

// for Nanopore
int long_minreadsize = 350;
map <pair<string,string>,int> long_score;
map <string,int> long_size;
map <string,int> long_subreads;


void read_sam_dna (int frame, string sam, string gene)
{
    ifstream input( sam.c_str() );
    for( std::string line; getline( input, line ); )
    {
        if (line.substr(0,1) == "@") {continue;}
        vector<string> sam_line;
        boost::split(sam_line,line,boost::is_any_of("\t"));
  
        int mm = 10000;
        if ((sam_line[2] == "*") || (sam_line[5] == "*")) {continue;}
        
        vector<string> nm_value;
        boost::split(nm_value,sam_line[12],boost::is_any_of(":"));
        mm = stoi(nm_value[2]);
        

        string cigar = splitcigar(sam_line[5]);
        vector<string> cigardata;
        boost::split(cigardata,cigar,boost::is_any_of(","));
        string c1 = cigardata[0];
        string c2 = cigardata[cigardata.size()-1];
        
        
        
        if ((c1 != "") && (mm != 10000))
        {
            if ((c1.substr(c1.size()-1,1) == "S") || (c1.substr(c1.size()-1,1) == "H"))
            {mm = mm + stoi(c1.substr(0,c1.size()));}
        }
        if ((c2 != "") && (mm != 10000))
        {
            if ((c2.substr(c2.size()-1,1) == "S") || (c2.substr(c2.size()-1,1) == "H"))
            {mm = mm + stoi(c2.substr(0,c2.size()));}
        }
        
        if ((sam_line[0].substr(sam_line[0].size()-2,2) == "/1") || (sam_line[0].substr(sam_line[0].size()-2,2) == "/2")){sam_line[0] = sam_line[0].substr(0,sam_line[0].size()-2);}

        
        pair <string,string> key = make_pair(gene,sam_line[0]);
        
        if (frame == 1) {
            if (sequence_list_r1.find(key) != sequence_list_r1.end())
            {
                if (sequence_list_r1[key] > mm + 1)
                {
                    sequence_list_r1[key] = mm + 1;
//                    sequence_size_r1[sam_line[0]] = sam_line[9].size();
                }
            }
            else {
                sequence_list_r1[key] = mm + 1;
 //               sequence_size_r1[sam_line[0]] = sam_line[9].size();
            }
        }
        
        if (frame == 2) {
            if (sequence_list_r2.find(key) != sequence_list_r2.end())
            {
                if (sequence_list_r2[key] > mm +1)
                {
                    sequence_list_r2[key] = mm +1;
   //                 sequence_size_r2[sam_line[0]] = sam_line[9].size();
                }
            }
            else {
                sequence_list_r2[key] = mm +1;
  //              sequence_size_r2[sam_line[0]] = sam_line[9].size();
            }
        }
        
    }
    input.close();
    removefile(sam,v_debug);
}








void main_dna_map ()
{

    int v_check = 0;
    if (v_r0 != "") {v_check = 1;}
    if ((v_r1 != "") && (v_r2 != "")) {v_check = 1;}
    if (v_bam != "") {v_check = 1;}
	

	if ((v_r1 != "") && (v_sample == "")) 
	{
		string basefile = v_r1;
		string basename = basefile.substr(basefile.find_last_of("/\\") + 1);
		v_sample = basename; 
		boost::algorithm::replace_all(v_sample, ".fastq.gz", "");
		boost::algorithm::replace_all(v_sample, ".fastq", "");
		boost::algorithm::replace_all(v_sample, ".fq.gz", "");
		boost::algorithm::replace_all(v_sample, ".fq", "");
		v_sample = v_sample + "_";
	}

	if ((v_r0 != "") && (v_sample == "")) 
	{
		string basefile = v_r0;
		string basename = basefile.substr(basefile.find_last_of("/\\") + 1);
		v_sample = basename; 
		boost::algorithm::replace_all(v_sample, ".fastq.gz", "");
		boost::algorithm::replace_all(v_sample, ".fastq", "");
		boost::algorithm::replace_all(v_sample, ".fq.gz", "");
		boost::algorithm::replace_all(v_sample, ".fq", "");
		v_sample = v_sample + "_";
	}


	if ((v_bam != "") && (v_sample == "")) 
	{
		string basefile = v_bam;
		string basename = basefile.substr(basefile.find_last_of("/\\") + 1);
		v_sample = basename;
		boost::algorithm::replace_all(v_sample, ".bam", "");
		boost::algorithm::replace_all(v_sample, ".BAM", "");
		v_sample = v_sample + "_";
	}
	
	
  
    if (((v_db == "") || (v_check == 0)) || (v_sample == ""))
    {
        screen_message (screen_size, 0, "", 1, v_quiet);
        v_message = "Program:   " + Program_name + "::map";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        v_message = "Version:   " + Program_version + ", " + Program_date;
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        v_message = "Usage:     kir-mapper map -r1 R1.gz -r2 R2.gz -sample sample_name <options>";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        v_message = "Usage:     kir-mapper map -bam your.bam -sample sample_name <options>";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Mandatory options:", 1, v_quiet);
        screen_message (screen_size, 2, "-bam         a BAM/CRAM file (ignore r0/r1/r2)", 1, v_quiet);
        screen_message (screen_size, 2, "-r1          a paired-ended forward fastq (fq, fastq, or gz)", 1, v_quiet);
        screen_message (screen_size, 2, "-r2          a paired-ended reverse fastq (fq, fastq, or gz)", 1, v_quiet);
        screen_message (screen_size, 2, "-r0          a single-ended fastq (fq, fastq, or gz)", 1, v_quiet);
        screen_message (screen_size, 2, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Other options:", 1, v_quiet);
        screen_message (screen_size, 2, "-sample      sample name [default: same as input file]", 1, v_quiet);
        screen_message (screen_size, 2, "-db          path to the kir-mapper database", 1, v_quiet);
        screen_message (screen_size, 2, "-output      output folder", 1, v_quiet);
        screen_message (screen_size, 2, "-threads     number of threads [" + v_threads + "]", 1, v_quiet);
        screen_message (screen_size, 2, "-buffer      number of sequences in buffer [" + to_string(v_buffer) + "]", 1, v_quiet);
        screen_message (screen_size, 2, "-error       the threshold for nucleotide quality trimming [" + to_string(v_mtrim_error).substr(0,4) + "]", 1, v_quiet);
        screen_message (screen_size, 2, "-tolerance   fraction of mismatches allowed [0.05 Illumina, 0.1 Nonapore]", 1, v_quiet);
        screen_message (screen_size, 2, "-downsample  downsampling for adjustment [" + to_string(downsampling) + "]", 1, v_quiet);

        screen_message (screen_size, 2, "", 1, v_quiet);
        screen_message (screen_size, 2, "--skip-unmapped   skip retrieving unmapped reads [not recommended]", 1, v_quiet);
        screen_message (screen_size, 2, "--skip-adjust     skip the adjustment procedure [not recommended]", 1, v_quiet);
        screen_message (screen_size, 2, "--low-mem         force low memory mode for sequence selection", 1, v_quiet);
        screen_message (screen_size, 2, "--exome           this is an exome (only exons)", 1, v_quiet);
        //screen_message (screen_size, 2, "--nanopore        this is nanopore data (must indicate -bam)", 1, v_quiet);
        //screen_message (screen_size, 2, "--rna             (beta) this is RNAseq data", 1, v_quiet);
        screen_message (screen_size, 2, "--quiet           quiet mode", 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        return;
    }
    
    if (v_sample == "") {
        v_r1 = "";
        v_r0 = "";
        v_r2 = "";
        v_bam = "";
        warnings.push_back("You must indicate a sample name or id (-sample sample_id)");
        main_dna_map();
        return;
    }
    
	if (v_nanopore == 1){
        if (v_bam == "")
		{
			v_r1 = "";
			v_r0 = "";
			v_r2 = "";
			warnings.push_back("You must indicate -bam when dealing with nanopore data");
			main_dna_map();
			return;
		}
	}

    
	
	
    if ((v_r1 != "") && (v_r1 == v_r2))
    {
        warnings.push_back("r1 and r2 must be different fastq files");
        main_dna_map();
        return;
    }

    if (v_bam != "") {
        
		string file = v_bam;
		boost::to_upper(file);
		if (! ends_with(file,".BAM")) 
		{
			warnings.push_back ("This is not a BAM file."); 
            v_sample = "";
            v_r1 = "";
            v_r0 = "";
            v_r2 = "";
            v_bam = "";
            main_dna_map();
			return;
		}
		
		
		ifstream file_bam(v_bam.c_str());
        if (!file_bam)
        {
            warnings.push_back ("Error accessing the .bam file.");
            v_sample = "";
            v_r1 = "";
            v_r0 = "";
            v_r2 = "";
            v_bam = "";
            main_dna_map();
            return;
        }
        file_bam.close();
    }

    if (v_r1 != "") {
		
		
		string file = v_r1;
		boost::to_upper(file);
		if (((ends_with(file,".FQ"))  || (ends_with(file,".FQ.GZ"))) || ((ends_with(file,".FASTQ"))  || (ends_with(file,".FASTQ.GZ"))))
		{}
		else {
			warnings.push_back ("This is not a FASTQ file."); 
            v_sample = "";
            v_r1 = "";
            v_r0 = "";
            v_r2 = "";
            v_bam = "";
            main_dna_map();
			return;
		}
		
        ifstream file_r1(v_r1.c_str());
        if (!file_r1)
        {
            warnings.push_back ("Error accessing the r1 file.");
            v_sample = "";
            v_r1 = "";
            v_r0 = "";
            main_dna_map();
            return;
        }
        file_r1.close();
    }

    if (v_r0 != "") {
		
		string file = v_r0;
		boost::to_upper(file);
		if (((ends_with(file,".FQ"))  || (ends_with(file,".FQ.GZ"))) || ((ends_with(file,".FASTQ"))  || (ends_with(file,".FASTQ.GZ"))))
		{}
		else {
			warnings.push_back ("This is not a FASTQ file."); 
            v_sample = "";
            v_r1 = "";
            v_r0 = "";
            v_r2 = "";
            v_bam = "";
            main_dna_map();
			return;
		}
		
        ifstream file_r1(v_r0.c_str());
        if (!file_r1)
        {
            warnings.push_back ("Error accessing the r0 file.");
            v_sample = "";
            v_r1 = "";
            v_r0 = "";
            main_dna_map();
            return;
        }
        file_r1.close();
    }
    
    
    if (v_r2 != "") {
        
		string file = v_r2;
		boost::to_upper(file);
		if (((ends_with(file,".FQ"))  || (ends_with(file,".FQ.GZ"))) || ((ends_with(file,".FASTQ"))  || (ends_with(file,".FASTQ.GZ"))))
		{}
		else {
			warnings.push_back ("This is not a FASTQ file."); 
            v_sample = "";
            v_r1 = "";
            v_r0 = "";
            v_r2 = "";
            v_bam = "";
            main_dna_map();
			return;
		}
		
		ifstream file_r2(v_r2.c_str());
        if (!file_r2)
        {
            warnings.push_back ("Error accessing the r2 file.");
            v_sample = "";
            v_r2 = "";
            main_dna_map();
            return;
        }
        file_r2.close();
    }
    
    
    if (((ends_with(v_r1,".fq")) && (ends_with(v_r2,".gz"))) || ((ends_with(v_r2,".fq")) && (ends_with(v_r1,".gz"))))
    {
        
        warnings.push_back ("kir-mapper cannot deal with .fastq and .gz files simultaneously.");
        v_sample = "";
        v_r1 = "";
        main_dna_map();
        return;
    }
    
    if (((ends_with(v_r1,".fastq")) && (ends_with(v_r2,".gz"))) || ((ends_with(v_r2,".fastq")) && (ends_with(v_r1,".gz"))))
    {
        warnings.push_back ("kir-mapper cannot deal with .fastq and .gz files simultaneously.");
        v_sample = "";
        v_r1 = "";
        main_dna_map();
        return;
    }
    
    if (((v_r1 != "") && (v_r2 == "")) || ((v_r1 == "") && (v_r2 != "")))
    {
        warnings.push_back ("You must indicate both r1 and r2 when in paired-end mode.");
        v_sample = "";
        v_r1 = "";
        main_dna_map();
        return;
    }
    
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
        main_dna_map();
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
        v_message = "This database is not compatible with this version of kir-mapper.";
        warnings.push_back (v_message);
        v_sample = "";
        v_r1 = "";
        main_dna_map();
        return;
    }
    
    
    if (v_bam != "")
    {
        if (v_bed != "")
        {
            ifstream file_bed(v_bed.c_str());
            if (!file_bed)
            {
                warnings.push_back ("Error accessing the BED file.");
                v_sample = "";
                v_r1 = "";
                v_bam = "";
                v_r2 = "";
                main_dna_map();
                return;
            }
            file_bed.close();
        }
    }
    
    
    
    
    v_callint = 1;
  
    string v_sample_sub = v_sample.substr(0, v_sample.size()-1);
  
    
    // criando output folder
    if (v_output == "")
    {
        if (v_bam != "") {v_output = findfilepath(v_bam) + "/kir-mapper/";}
        if (v_r0 != "") {v_output = findfilepath(v_r0) + "/kir-mapper/";}
        if (v_r1 != "") {v_output = findfilepath(v_r1) + "/kir-mapper/";}
    }
    
    v_command = "mkdir " + v_output;
    v_system_out = GetStdoutFromCommand(v_command);
    v_command = "mkdir " + v_output + "/map/";
    v_system_out = GetStdoutFromCommand(v_command);
    v_command = "mkdir " + v_output + "/map/" + v_sample_sub;
    v_system_out = GetStdoutFromCommand(v_command);
    v_command = "mkdir " + v_output + "/map/" + v_sample_sub + "/log";
    v_system_out = GetStdoutFromCommand(v_command);
 
    v_output = v_output + "/map/" + v_sample_sub + "/";
    
    
    if (! fileExists(v_output))
    {
        warnings.push_back ("Error creating the output folder.");
        v_sample = "";
        v_r1 = "";
        v_bam = "";
        v_r2 = "";
        main_dna_map();
        return;
    }
    
 
    string v_general_log = v_output + v_sample_sub + ".kir-mapper.log";
    ofstream general_log;
    general_log.open (v_general_log.c_str());
    
    srand(time(NULL));
    
    screen_message (screen_size, 0, "", 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);
    v_message = Program_name + "::map";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    general_log << v_message << endl;
 
    v_message = "Version " + Program_version;
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    general_log << v_message << endl;
   

    v_message = "Number of threads: " + v_threads;
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    v_message = "Sample: " + v_sample_sub;
    screen_message (screen_size, 0, v_message, 1, v_quiet);

       
    if (v_exome == 1) {
        v_message = "Mode: only exons or exome";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
		general_log << v_message << endl;
    }
    
	if (((v_exome == 0) && (v_rna == 0)) && (v_nanopore == 0)) {
        v_message = "Mode: full genome";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
		general_log << v_message << endl;
    }
	
	if (((v_exome == 0) && (v_rna == 0)) && (v_nanopore == 1)) {
        v_message = "Mode: Nanopore";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
		general_log << v_message << endl;
		v_map_type = "single";
    }
	
    if (v_rna == 1) {
        v_message = "Mode: RNAseq";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
		general_log << v_message << endl;
    }
    
    v_message = "Sample: " + v_sample_sub;
    general_log << v_message << endl;
    v_message = "Database: " + v_db;
    general_log << v_message << endl;
    v_message = "Output: " + v_output;
    general_log << v_message << endl;
    
    v_message = "Files (R0): " + v_r0;
    general_log << v_message << endl;
    v_message = "Files (R1): " + v_r1;
    general_log << v_message << endl;
    v_message = "Files (R2): " + v_r2;
    general_log << v_message << endl;
    v_message = "Files (BAM): " + v_bam;
    general_log << v_message << endl;

        
    v_message = "Trimming error: " + to_string(v_mtrim_error);
    general_log << v_message << endl;
    v_message = "Tolerance: " + to_string(v_tolerance);
    general_log << v_message << endl;
    v_message = "Minimum read size: " + to_string(v_size);
    general_log << v_message << endl;

    if (v_bam != "")
    {
        if (v_skip_unmapped == 1) {
            v_message = "Retrieving unmapped: no";
        }
        if (v_skip_unmapped == 0) {
            v_message = "Retrieving unmapped: yes";
        }
        general_log << v_message << endl;
    }

    if (v_skiptyping == 1) {
        v_message = "Performing adjustments: no";
    }
    if (v_skiptyping == 0) {
        v_message = "Performing adjustments: yes";
    }
    general_log << v_message << endl;
    


    
    main_preselect();
	
    if (v_callint == 2) {return;}
//    general_log << v_message << endl;
    
    
	
	
	

    
	screen_message (screen_size, 0, "Sorting sequences ...", 2, v_quiet);

	int loop = stoi(v_threads);
	boost::algorithm::erase_all(v_genes, " ");
	boost::split(v_gene_list,v_genes,boost::is_any_of(","));

	string selectcleanR1 = v_output + v_sample + "selected.trim.R1.fastq";
	string selectcleanR2 = v_output + v_sample + "selected.trim.R2.fastq";
	if (v_map_type == "single") {selectcleanR1 = v_output + v_sample + "selected.trim.R0.fastq";}
	
	if (v_nanopore == 1) {
		v_map_type == "single";
		selectcleanR1 = v_output + v_sample + "selected.R0.fastq";
	}
		
	selected_results_r1.clear();
	selected_results_r2.clear();
	selected_results_trim_r1.clear();
	selected_results_trim_r2.clear();
	
	map <string,int> r1_trim_id_map;
	unordered_map <string,string> r1_trim_seq_map;
	unordered_map <string,string> r2_trim_seq_map;
	unordered_map <string,string> r1_trim_qual_map;
	unordered_map <string,string> r2_trim_qual_map;

	ifstream r1sort (selectcleanR1.c_str());
	ifstream r2sort (selectcleanR2.c_str());
	
	double length_sum;
	double seqcount;
	
	for( std::string line; getline( r1sort, line ); )
	{
		string idr1 = line;
		getline( r1sort, line );
		string seqr1 = line;
		
		getline( r1sort, line );
		string infor1 = line;
		getline( r1sort, line );
		string qualr1 = line;

		getline( r2sort, line );
		string idr2 = line;
		getline( r2sort, line );
		string seqr2 = line;
		getline( r2sort, line );
		string infor2 = line;
		getline( r2sort, line );
		string qualr2 = line;

		r1_trim_id_map[idr1] = 1;
		r1_trim_seq_map[idr1] = seqr1;
		r2_trim_seq_map[idr1] = seqr2;
		r1_trim_qual_map[idr1] = qualr1;
		r2_trim_qual_map[idr1] = qualr2;
		
		length_sum = length_sum + double(seqr1.size());
		seqcount++;
		
	}
	r1sort.close();
	r2sort.close();
	
	
	
	
	if (v_nanopore == 0) {
		int readsize = int(length_sum / seqcount);
		v_message = "Mean read size after quality trimming: " + to_string(readsize);
		screen_message (screen_size, 0, v_message, 1, v_quiet);
		general_log << v_message << endl;
		
		if (readsize < 100) 
		{
			v_message = " > Warning: mean read size after quality trimming < 100.";
			screen_message (screen_size, 0, v_message, 1, 0);
			general_log << v_message << endl;

			if (readsize < 80) {
				v_message = " > Warning: We do not recommend using kir-mapper with this data.";
				screen_message (screen_size, 0, v_message, 1, 0);
				general_log << v_message << endl;

			}
			v_message = " > Warning: Results might be biased due to low read size.";
			screen_message (screen_size, 0, v_message, 1, 0);
			general_log << v_message << endl;
		}
	}
	
	
	
	
	for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
	{
		if (*i ==  "") {continue;}
		screen_message (screen_size, 0, "Sorting sequences for " + *i + " ...", 2, v_quiet);
		
		unordered_map <string,int> motifs;

		string motifdbfile = "";
		if (v_rna == 1) {motifdbfile = v_db + "/motif/rna/loci/" + *i + ".txt";}
		if (v_rna == 0) {motifdbfile = v_db + "/motif/dna/loci/" + *i + ".txt";}
		

		boost::replace_all(motifdbfile, "\\ ", " ");

		ifstream motifdb (motifdbfile.c_str());
		for( std::string line; getline( motifdb, line ); )
		{
			if (line == "") {continue;}
			motifs[line] = 1;
		}
		motifdb.close();
		
		vector <string> select_data_sorted_r1;
		vector <string> select_data_sorted_r2;

		if (v_map_type == "paired") {
 
			ThreadPool pool(loop);
			std::vector< std::future<int> > results;
			  
			for( auto & item : r1_trim_id_map)
			  {
				  string id = item.first;
				  string r1seq = r1_trim_seq_map[id];
				  string r2seq = r2_trim_seq_map[id];
				  string r1qual = r1_trim_qual_map[id];
				  string r2qual = r2_trim_qual_map[id];

				  results.emplace_back(
										 pool.enqueue([id, r1seq, r2seq, r1qual, r2qual, &motifs, &select_data_sorted_r1, &select_data_sorted_r2]
						  {
						  if (r1seq.size() < motif_size) {return 1;}
						  if (r2seq.size() < motif_size) {return 1;}
						  for (int c = 0; c < (r1seq.size() - motif_size); c++)
						  {
							  string sub = r1seq.substr(c,motif_size);
							  if ( motifs.find(sub) != motifs.end() ) {
								  mtx_selection_general.lock();
								  select_data_sorted_r1.push_back(id + "\n" + r1seq + "\n+\n" + r1qual);
								  select_data_sorted_r2.push_back(id + "\n" + r2seq + "\n+\n" + r2qual);
								  mtx_selection_general.unlock();
								  return 1;
							  }
							}
							  return 1;
						  })
						  );
			  }
			  for(auto && result: results){result.get();} // waiting for all threads
			
			ofstream OUTselectR1;
			string selectR1 = v_output + v_sample + "sorted_" + *i + "_R1.fastq";
			OUTselectR1.open (selectR1.c_str());

			ofstream OUTselectR2;
			string selectR2 = v_output + v_sample + "sorted_" + *i + "_R2.fastq";
			OUTselectR2.open (selectR2.c_str());
			
			for (int a = 0; a < select_data_sorted_r1.size(); a++)
			{
				  OUTselectR1 << select_data_sorted_r1[a] << endl;
				  OUTselectR2 << select_data_sorted_r2[a] << endl;
			}
			OUTselectR1.close();
			OUTselectR2.close();
			select_data_sorted_r1.clear();
			select_data_sorted_r2.clear();

			screen_message (screen_size, 0, "Sorting sequences for " + *i + " ... done", 2, v_quiet);
			
		}
		

 
		   if (v_map_type == "single") {
	
			   ThreadPool pool(loop);
			   std::vector< std::future<int> > results;
				 
			   for( auto & item : r1_trim_id_map)
				 {
					 string id = item.first;
					 string r1seq = r1_trim_seq_map[id];
					 string r1qual = r1_trim_qual_map[id];
 
					 results.emplace_back(
											pool.enqueue([id, r1seq, r1qual, &motifs, &select_data_sorted_r1]
							 {
							 if (r1seq.size() < motif_size) {return 1;}
							 for (int c = 0; c < (r1seq.size() - motif_size); c++)
							 {
								 string sub = r1seq.substr(c,motif_size);
								 if ( motifs.find(sub) != motifs.end() ) {
									 mtx_selection_general.lock();
									 select_data_sorted_r1.push_back(id + "\n" + r1seq + "\n+\n" + r1qual);
									 mtx_selection_general.unlock();
									 return 1;
								 }
							   }
								 return 1;
							 })
							 );
				 }
				 for(auto && result: results){result.get();} // waiting for all threads
			   
			   ofstream OUTselectR1;
			   string selectR1 = v_output + v_sample + "sorted_" + *i + "_R0.fastq";
			   OUTselectR1.open (selectR1.c_str());

			   for (int a = 0; a < select_data_sorted_r1.size(); a++)
			   {
					 OUTselectR1 << select_data_sorted_r1[a] << endl;
			   }
			   OUTselectR1.close();
			   select_data_sorted_r1.clear();
			   select_data_sorted_r2.clear();

			   screen_message (screen_size, 0, "Sorting sequences for " + *i + " ... done", 2, v_quiet);
			   
		   }
		
	}

	screen_message (screen_size, 0, "Sorting sequences: done", 1, v_quiet);
    

 
    


    
    
 

    map <string,float> locus_depth;
    std::vector< std::future<int> > presence_result;
    map <string,int> locus_presence;

    if (v_rna == 0) {
        screen_message (screen_size, 0, "Detecting KIR genes ...", 2, v_quiet);
        ThreadPool presence_pool(stoi(v_threads));
        int loopp = 0;

        
        for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
        {
            if (*i ==  "") {continue;}
            string gene = *i;
            loopp++;
            string motif_list = v_db + "/presence/cds/" + gene + ".txt";
            if (! fileExists(motif_list)) {continue;}
            
            presence_result.emplace_back(
                                        presence_pool.enqueue([loopp, gene,&locus_depth]
            {
            

                string motif_list = v_db + "/presence/cds/" + gene + ".txt";
                ifstream file( motif_list.c_str());
                map <string,float> motifs;
                for( std::string line; getline( file, line ); )
                {
                    if (line == "") {continue;}
                    vector <string> data;
                    boost::split(data,line,boost::is_any_of("\t"));
                    motifs[data[1]] = 0;
//                    string rev = reverse_and_complement(data[1]);
//                    motifs[rev] = 0;
                }
                file.close();
                
				string selectR1 = v_output + v_sample + "sorted_" + gene + "_R1.fastq";
				string selectR2 = v_output + v_sample + "sorted_" + gene + "_R2.fastq";
				string selectR0 = v_output + v_sample + "sorted_" + gene + "_R0.fastq";
                
                if (v_map_type == "paired")
                {
                    ifstream file( selectR1.c_str());
                    for( std::string line; getline( file, line ); )
                    {
                        if (line == "") {continue;}
                        string id = line;
                        getline( file, line );
                        string seq = line;
                        getline( file, line );
                        getline( file, line );
                        for (int a = 0; a < (seq.size()-30); a++)
                        {
                            string motif = seq.substr(a,30);
                            if (motifs.find(motif) != motifs.end())
                            {
                                motifs[motif]++;
                            }
                        }
						string revseq = reverse_and_complement(seq);
						for (int a = 0; a < (revseq.size()-30); a++)
                        {
                            string motif = revseq.substr(a,30);
                            if (motifs.find(motif) != motifs.end())
                            {
                                motifs[motif]++;
                            }
                        }
                        
                    }
                    file.close();
                }
                
				if (v_map_type == "paired")
                {
                    ifstream file( selectR2.c_str());
                    for( std::string line; getline( file, line ); )
                    {
                        if (line == "") {continue;}
                        string id = line;
                        getline( file, line );
                        string seq = line;
                        getline( file, line );
                        getline( file, line );
                        for (int a = 0; a < (seq.size()-30); a++)
                        {
                            string motif = seq.substr(a,30);
                            if (motifs.find(motif) != motifs.end())
                            {
                                motifs[motif]++;
                            }
                        }
						string revseq = reverse_and_complement(seq);
						for (int a = 0; a < (revseq.size()-30); a++)
                        {
                            string motif = revseq.substr(a,30);
                            if (motifs.find(motif) != motifs.end())
                            {
                                motifs[motif]++;
                            }
                        }
                        
                    }
                    file.close();
                }
                
				if (v_map_type == "single")
                {
                    ifstream file( selectR0.c_str());
                    for( std::string line; getline( file, line ); )
                    {
                        if (line == "") {continue;}
                        string id = line;
                        getline( file, line );
                        string seq = line;
                        getline( file, line );
                        getline( file, line );
                        for (int a = 0; a < (seq.size()-30); a++)
                        {
                            string motif = seq.substr(a,30);
                            if (motifs.find(motif) != motifs.end())
                            {
                                motifs[motif]++;
                            }
                        }
						string revseq = reverse_and_complement(seq);
						for (int a = 0; a < (revseq.size()-30); a++)
                        {
                            string motif = revseq.substr(a,30);
                            if (motifs.find(motif) != motifs.end())
                            {
                                motifs[motif]++;
                            }
                        }
                        
                    }
                    file.close();
                }
                
                float sum = 0;
                float count = 0;
                for (auto item : motifs)
                {
                    count++;
                    sum = sum + item.second;
                }
                float ratio = sum / count;
                
                mtxg.lock();
                locus_depth[gene] = ratio;
                mtxg.unlock();
                return 1;
            })
            );
            
        }
        for(auto && result: presence_result){result.get();} // waiting for all threads
        screen_message (screen_size, 0, "Detecting KIR genes: done", 1, v_quiet);
        

        v_message = "Depth estimation for KIR3DL3: about " + to_string(int(locus_depth["KIR3DL3"]));
        screen_message (screen_size, 0, v_message, 1, 0);
		general_log << v_message << endl;


        if (locus_depth["KIR3DL3"] < 5)
        {
            v_message = " Very low depth for KIR3DL3 (< 5). Quitting!";
            warnings.push_back (v_message);
            return;
        }


        if ((locus_depth["KIR3DL3"] < 20) && (v_exome == 0))
        {
            v_message = " > Warning: low depth, about " + to_string(int(locus_depth["KIR3DL3"])) + " for KIR3DL3";
            screen_message (screen_size, 0, v_message, 1, 0);
			general_log << v_message << endl;
		
            v_message = " > Warning: minimum recommended is 20 for WGS";
            screen_message (screen_size, 0, v_message, 1, 0);
			
			v_message = " > Warning: Results might be biased due to low depth";
			screen_message (screen_size, 0, v_message, 1, 0);
			general_log << v_message << endl;
        }
		
		if ((locus_depth["KIR3DL3"] < 50) && (v_exome == 1))
        {
            v_message = " > Warning: low depth, about " + to_string(int(locus_depth["KIR3DL3"])) + " for KIR3DL3";
            screen_message (screen_size, 0, v_message, 1, 0);
			general_log << v_message << endl;
		
            v_message = " > Warning: minimum recommended is 50 for Exomes";
            screen_message (screen_size, 0, v_message, 1, 0);
			
			v_message = " > Warning: Results might be biased due to low depth";
			screen_message (screen_size, 0, v_message, 1, 0);
			general_log << v_message << endl;
        }


        string presence_report = v_output + v_sample + "presence_report.txt";
        ofstream OUTLOCUS;
        OUTLOCUS.open (presence_report.c_str());
        OUTLOCUS << "GENE\tRATIO\tSTATUS" << endl;

        float ref = locus_depth["KIR3DL3"];
        for (auto item : locus_depth)
        {
            float ratio = item.second / ref;
            if (ratio >= 0.15) {
                OUTLOCUS << item.first << "\t" << ratio << "\tpresent" << endl;
                locus_presence[item.first] = 1;
            }
            else {
                locus_presence[item.first] = 0;
                OUTLOCUS << item.first << "\t" << ratio << "\tabsent" << endl;
            }
            if (v_debug == 1)
            {
                cout << item.first << "\t" << ratio << endl;
            }
        }
        OUTLOCUS.close();
    }
    









	screen_message (screen_size, 0, "Scoring sequences ...", 2, v_quiet);

	boost::algorithm::erase_all(v_genes, " ");
	boost::split(v_gene_list,v_genes,boost::is_any_of(","));
    
	
	if (v_nanopore == 1) 
	{
		for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
		{
			if (*i ==  "") {continue;}
			
			if (v_rna == 0) {
				if (locus_presence.find(*i) != locus_presence.end())
				{
					if (locus_presence[*i] == 0) {continue;}
				}
			}
			screen_message (screen_size, 0, "Scoring sequences for " + *i + " ...", 2, v_quiet);
			
			string v_ref = "'" + v_db + "/mapper/dna-nanopore/select/" + *i + "/" + *i + ".fas' "; 
			boost::replace_all(v_ref, "\\ ", " ");
			string v_r0_sort = v_output + v_sample + "sorted_" + *i + "_R0.fastq";

			string v_log = v_output + "/log/" + *i + "_minimap.log";
			ofstream minilog;
			minilog.open (v_log.c_str());
				
			string outsamtmp =  v_output + v_sample + *i + ".tmp.sam";
			v_command = v_minimap + " -a -t " + v_threads + " -x map-ont -o " + outsamtmp + " " + v_ref + " " + v_r0_sort;
			debug_message(v_command);
			string data = GetStdoutFromCommand(v_command.c_str());
			minilog << data;
			minilog.close();
			
			
			string outfqtmp =  v_output + v_sample + *i + ".tmp.fq";
			debug_message("openning " + outfqtmp);
			ofstream fqtmp;
			fqtmp.open (outfqtmp.c_str());
			
			ifstream SAM( outsamtmp.c_str());
			int linenumber = 1;
			for( std::string line; getline( SAM, line ); )
			{
				string item = line;
				debug_message(*i + ", sam line number: " + to_string(linenumber));
				linenumber++;
				if (item == "") {continue;}
				if (item.substr(0,1) == "@") {continue;}
				if (item.substr(0,1) == "[") {continue;}
				vector <string> fields;
				boost::split(fields,item,boost::is_any_of("\t"));
				if (fields.size() < 11) {continue;}
				
				if ((fields[1] != "0") && (fields[1] != "16")) {continue;}
				
				
				string read = fields[0];
				string seq = fields[9];
				string qual = fields[10];
				string nm = fields[11];
				string cigar = fields[5];
				
				std::size_t found = nm.find("NM:");
				if (found==std::string::npos) {continue;}
				boost::replace_all(nm, "NM:i:", "");
				debug_message("Testing NM start");
				int nmint = stoi(nm);
				debug_message("Testing NM end");

				
				boost::replace_all(cigar, "S", ",S;");
				boost::replace_all(cigar, "H", ",H;");
				boost::replace_all(cigar, "M", ",M;");
				boost::replace_all(cigar, "I", ",I;");
				boost::replace_all(cigar, "D", ",D;");
				boost::replace_all(cigar, "X", ",X;");
				

				vector <string> cigars;
				boost::split(cigars,cigar,boost::is_any_of(";"));
				
				int upperH = 0;
				int upperS = 0;
				int lowerH = 0;
				int lowerS = 0;
				
				debug_message("Finding coordinates start");
				for (auto subcigar : cigars)
				{
					if (subcigar == "") {continue;}
					vector <string> sub;
					boost::split(sub,subcigar,boost::is_any_of(","));
					if (sub.size() != 2) {continue;}
					if (! stoi(sub[0])) {continue;}
					if ((sub[1] == "H") && (upperH == 0))	{upperH = stoi(sub[0]);continue;}
					if ((sub[1] == "S") && (upperS == 0))	{upperS = stoi(sub[0]);continue;}
					if ((sub[1] == "S") && (lowerS == 0))	{lowerS = stoi(sub[0]);continue;}
					if ((sub[1] == "H") && (lowerH == 0))	{lowerH = stoi(sub[0]);continue;}
				}
				debug_message("Finding coordinates end");


				debug_message("Ratio calculation start ");
				int cutupper = upperH + upperS;
				int cutlower = lowerH + lowerS;
				int size = seq.length() - (cutupper + cutlower);
				if (size <= 0) {continue;}
				debug_message("Sub sequence size: " + to_string(size));
				string subseq = seq.substr(cutupper,size);
				string subqual = qual.substr(cutupper,size);
				if (subseq.length() < long_minreadsize) {continue;}
				float ratio = float(nmint) / float(subseq.length());
				debug_message("Ratio: " + to_string(ratio));
				if (ratio > v_tolerance) {continue;}
				debug_message("Ratio calculation end ");
				
				if (fields[1] != "16")
				{
					string rev = reverse_and_complement(subseq);
					subseq = rev;
					
					string revqual = subqual; 
					reverse(revqual.begin(), revqual.end()); 
					subqual = revqual;
				}

	
				debug_message("writing fastq tmp start ");
				string newreadname = read + ";" + *i;
				long_subreads[newreadname] = 1;
				fqtmp << "@" << newreadname << endl;
				fqtmp << subseq << endl;
				fqtmp << "+" << endl;
				fqtmp << subqual << endl;
				
				pair <string,string> key = make_pair (newreadname,*i);
				
				long_score[key] = nmint;
				long_size[newreadname] = size;
				debug_message("writing fastq tmp end ");
			}
			fqtmp.close();
			SAM.close();
			removefile(outsamtmp,v_debug);
			
			
			for( std::vector<string>::const_iterator j = v_gene_list.begin(); j != v_gene_list.end(); ++j)
			{
				if (*j ==  "") {continue;}
				if (*j ==  *i) {continue;}
				
				if (v_rna == 0) {
					if (locus_presence.find(*j) != locus_presence.end())
					{
						if (locus_presence[*j] == 0) {continue;}
					}
				}
				
				string v_ref = "'" + v_db + "/mapper/dna-nanopore/score/" + *j + "/" + *j + ".fas' "; 
				boost::replace_all(v_ref, "\\ ", " ");
			
				string outsamtmp =  v_output + v_sample + *i + "." + *j + ".tmp.sam";
				v_command = v_minimap + " -a -t " + v_threads + " -x map-ont -o " + outsamtmp + " " + v_ref + " " + outfqtmp;
				debug_message(v_command);
				string data = GetStdoutFromCommand(v_command.c_str());
				
				
				ifstream SAMSEC( outsamtmp.c_str());
				int linenumber = 1;
				for( std::string line; getline( SAMSEC, line ); )
				{
					string item = line;
					debug_message(*j + ", sam line number: " + to_string(linenumber));
					linenumber++;
					if (item == "") {continue;}
					if (item.substr(0,1) == "@") {continue;}
					if (item.substr(0,1) == "[") {continue;}
					vector <string> fields;
					boost::split(fields,item,boost::is_any_of("\t"));
					if (fields.size() < 11) {continue;}
					
					if ((fields[1] != "0") && (fields[1] != "16")) {continue;}
					
					
					string read = fields[0];
					string nm = fields[11];
					string cigar = fields[5];
				
					boost::replace_all(cigar, "S", ",S;");
					boost::replace_all(cigar, "H", ",H;");
					boost::replace_all(cigar, "M", ",M;");
					boost::replace_all(cigar, "I", ",I;");
					boost::replace_all(cigar, "D", ",D;");
					boost::replace_all(cigar, "X", ",X;");
				

					vector <string> cigars;
					boost::split(cigars,cigar,boost::is_any_of(";"));
				
					int upperH = 0;
					int upperS = 0;
					int lowerH = 0;
					int lowerS = 0;
				
					debug_message("Finding coordinates start");
					for (auto subcigar : cigars)
					{
						if (subcigar == "") {continue;}
						vector <string> sub;
						boost::split(sub,subcigar,boost::is_any_of(","));
						if (sub.size() != 2) {continue;}
						if (! stoi(sub[0])) {continue;}
						if ((sub[1] == "H") && (upperH == 0))	{upperH = stoi(sub[0]);continue;}
						if ((sub[1] == "S") && (upperS == 0))	{upperS = stoi(sub[0]);continue;}
						if ((sub[1] == "S") && (lowerS == 0))	{lowerS = stoi(sub[0]);continue;}
						if ((sub[1] == "H") && (lowerH == 0))	{lowerH = stoi(sub[0]);continue;}
					}
					debug_message("Finding coordinates end");
					
					
					std::size_t found = nm.find("NM:");
					if (found==std::string::npos) {continue;}
					boost::replace_all(nm, "NM:i:", "");
					debug_message("Testing NM start");
					int nmint = stoi(nm);
					nmint = nmint + upperH + upperS + lowerS + lowerH;
			
					pair <string,string> key1 = make_pair (read,*j);
					long_score[key1] = nmint;
				}
				SAMSEC.close();
				removefile(outsamtmp,v_debug);
			}
			
				
			v_ref = "'" + v_db + "/others/dna/kir_others.fasta' "; 
			boost::replace_all(v_ref, "\\ ", " ");
		
			outsamtmp =  v_output + v_sample + *i + ".others.tmp.sam";
			v_command = v_minimap + " -a -t " + v_threads + " -x map-ont -o " + outsamtmp + " " + v_ref + " " + outfqtmp;
			debug_message(v_command);
			data = GetStdoutFromCommand(v_command.c_str());
			
			
			ifstream SAMSEC( outsamtmp.c_str());
			linenumber = 1;
			for( std::string line; getline( SAMSEC, line ); )
			{
				string item = line;
				debug_message("others, sam line number: " + to_string(linenumber));
				linenumber++;
				if (item == "") {continue;}
				if (item.substr(0,1) == "@") {continue;}
				if (item.substr(0,1) == "[") {continue;}
				vector <string> fields;
				boost::split(fields,item,boost::is_any_of("\t"));
				if (fields.size() < 11) {continue;}
				
				if ((fields[1] != "0") && (fields[1] != "16")) {continue;}
				
				
				string read = fields[0];
				string nm = fields[11];
				string cigar = fields[5];
			
				boost::replace_all(cigar, "S", ",S;");
				boost::replace_all(cigar, "H", ",H;");
				boost::replace_all(cigar, "M", ",M;");
				boost::replace_all(cigar, "I", ",I;");
				boost::replace_all(cigar, "D", ",D;");
				boost::replace_all(cigar, "X", ",X;");
			

				vector <string> cigars;
				boost::split(cigars,cigar,boost::is_any_of(";"));
			
				int upperH = 0;
				int upperS = 0;
				int lowerH = 0;
				int lowerS = 0;
			
				debug_message("Finding coordinates start");
				for (auto subcigar : cigars)
				{
					if (subcigar == "") {continue;}
					vector <string> sub;
					boost::split(sub,subcigar,boost::is_any_of(","));
					if (sub.size() != 2) {continue;}
					if (! stoi(sub[0])) {continue;}
					if ((sub[1] == "H") && (upperH == 0))	{upperH = stoi(sub[0]);continue;}
					if ((sub[1] == "S") && (upperS == 0))	{upperS = stoi(sub[0]);continue;}
					if ((sub[1] == "S") && (lowerS == 0))	{lowerS = stoi(sub[0]);continue;}
					if ((sub[1] == "H") && (lowerH == 0))	{lowerH = stoi(sub[0]);continue;}
				}
				debug_message("Finding coordinates end");
				
				
				std::size_t found = nm.find("NM:");
				if (found==std::string::npos) {continue;}
				boost::replace_all(nm, "NM:i:", "");
				debug_message("Testing NM start");
				int nmint = stoi(nm);
				nmint = nmint + upperH + upperS + lowerS + lowerH;
		
				pair <string,string> key1 = make_pair (read,"others");
				long_score[key1] = nmint;
			}
			SAMSEC.close();
			removefile(outsamtmp,v_debug);
			
			
			
			
			
			
			
			
		}
		/*
		string testfile = v_output + v_sample + ".testelist.txt";
		ofstream test;
		test.open (testfile.c_str());
		for (auto item : long_score)
		{
			test << item.first.first << "\t" << item.first.second << "\t" << item.second << endl;
		}
		test.close();
		*/
    
	}
	

	
	
	if (v_nanopore == 0) 
	{
	
		for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
		{
			if (*i ==  "") {continue;}
			
			if (v_rna == 0) {
				if (locus_presence.find(*i) != locus_presence.end())
				{
					if (locus_presence[*i] == 0) {continue;}
				}
			}
			screen_message (screen_size, 0, "Scoring sequences for " + *i + " ...", 2, v_quiet);
			
		   string v_ref = "";
		   if (v_rna == 1) {
				v_ref = "'" + v_db + "/mapper/rna/" + *i + "/" + *i + ".fas' ";
		   }
		   if (v_rna == 0) {
				v_ref = "'" + v_db + "/mapper/dna/" + *i + "/" + *i + ".fas' "; 
		   }
			
		   boost::replace_all(v_ref, "\\ ", " ");
		   string v_sam1 = v_output + *i + ".mapper.1.sam";
		   string v_sam2 = v_output + *i + ".mapper.2.sam";
			   
		   string v_r0_sort = v_output + v_sample + "sorted_" + *i + "_R0.fastq";
		   string v_r1_sort = v_output + v_sample + "sorted_" + *i + "_R1.fastq";
		   string v_r2_sort = v_output + v_sample + "sorted_" + *i + "_R2.fastq";
			
		   string v_sai1 = v_output + "tmp1.sai ";
		   string v_sai2 = v_output + "tmp2.sai ";


		   string v_log = v_output + "/log/" + *i + "_mapper_aln.1.log";
		   v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_r1_sort + " > " + v_sai1 + " 2>" + v_log;
		   
			if (v_map_type == "single") {v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_r0_sort + " > " + v_sai1 + " 2>" + v_log;}
		   system (v_command.c_str());
		   
		   v_log = v_output + "/log/" + *i + "_mapper_sam.1.log";
		   v_command = v_bwa + " samse " + v_ref + v_sai1 + v_r1_sort + " > " + v_sam1 + " 2>" + v_log;
		   
		   if (v_map_type == "single") {v_command = v_bwa + " samse " + v_ref + v_sai1 + v_r0_sort + " > " + v_sam1 + " 2>" + v_log;}
		   system (v_command.c_str());
		   v_command = " rm " + v_sai1;
		   system (v_command.c_str());

			if (v_map_type == "paired") {
				v_log = v_output + "/log/" + *i + "_mapper_aln.2.log";
				v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_r2_sort + " > " + v_sai2 + " 2>" + v_log;
				system (v_command.c_str());
				v_log = v_output + "/log/" + *i + "_mapper_sam.2.log";
				v_command = v_bwa + " samse " + v_ref + v_sai2 + v_r2_sort + " > " + v_sam2 + " 2>" + v_log;
				system (v_command.c_str());
				v_command = " rm " + v_sai2;
				system (v_command.c_str());
			}
			
			thread r1 (read_sam_dna, 1, v_sam1, *i);
			thread r2 (read_sam_dna, 2, v_sam2, *i);
			r1.join();
			r2.join();
			
		}
		
	  
		v_message = "Scoring reads considering other genes ...";
		screen_message (screen_size, 0, v_message, 2, v_quiet);

		 
		string v_ref = "'" + v_db + "/others/dna/kir_others.fasta' ";
		boost::replace_all(v_ref, "\\ ", " ");

		string v_sam1 = v_output + "others.1.sam";
		string v_sam2 = v_output + "others.2.sam";
		string v_r0_sort = v_output + v_sample + "selected.trim.R0.fastq";
		string v_r1_sort = v_output + v_sample + "selected.trim.R1.fastq";
		string v_r2_sort = v_output + v_sample + "selected.trim.R2.fastq";
		string v_sai1 = v_output + "tmp1.sai ";
		string v_sai2 = v_output + "tmp2.sai ";

		 
		string v_log = v_output + "/log/others_aln.1.log";
		v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_r1_sort + " > " + v_sai1 + " 2>" + v_log;
		if (v_map_type == "single") {v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_r0_sort + " > " + v_sai1 + " 2>" + v_log;}
		system (v_command.c_str());
		v_log = v_output + "/log/others_sam.1.log";
		v_command = v_bwa + " samse " + v_ref + v_sai1 + v_r1_sort + " > " + v_sam1 + " 2>" + v_log;
		if (v_map_type == "single") {v_command = v_bwa + " samse " + v_ref + v_sai1 + v_r0_sort + " > " + v_sam1 + " 2>" + v_log;}
		system (v_command.c_str());
		v_command = " rm " + v_sai1;
		system (v_command.c_str());

		if (v_map_type == "paired"){
			v_log = v_output + "/log/others_aln.2.log";
			v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_r2_sort + " > " + v_sai2 + " 2>" + v_log;
			system (v_command.c_str());
			v_log = v_output + "/log/others_sam.2.log";
			v_command = v_bwa + " samse " + v_ref + v_sai2 + v_r2_sort + " > " + v_sam2 + " 2>" + v_log;
			system (v_command.c_str());
			v_command = " rm " + v_sai2;
			system (v_command.c_str());
		}
		thread r1 (read_sam_dna, 1, v_sam1, "others");
		thread r2 (read_sam_dna, 2, v_sam2, "others");
		r1.join();
		r2.join();
		 

    }
	
	v_message = "Scoring reads: done";
	screen_message (screen_size, 0, v_message, 1, v_quiet);
    
    
    



   

    
    
    
    
    
    
    
    
    
    
    //Addressing sequences
	
	
     v_message = "Comparing scores ...";
     screen_message (screen_size, 0, v_message, 2, v_quiet);
	 map <string,string> address_list_long;
	 std::map <string, string> address_list_r1;
     
	 if (v_nanopore == 1)
	 {
		vector<string> v_target_list_sub;
		boost::algorithm::erase_all(v_genes, " ");
		boost::split(v_target_list_sub,v_genes,boost::is_any_of(","));
		
		string v_add = v_output + v_sample + "addressing_table.txt";
		ofstream ADD;
		ADD.open (v_add.c_str());
			
		ADD << "subRead\tSize";
		for (auto gene : v_target_list_sub)
		{
			if (gene == "") {continue;}
			ADD << "\t" << gene;
		}
		ADD << "\tOthers";
		ADD << "\tTarget" << endl;
		
		
		vector <string> v_target_list_sub_tmp = v_target_list_sub;
		v_target_list_sub_tmp.push_back("others");
		
		for (auto item : long_subreads)
		{
			string read = item.first;
			
			ADD << read;
			vector <string> read_data;
			boost::split(read_data,read,boost::is_any_of(";"));
			string valid = read_data[1];
			ADD << "\t" << long_size[read];
			
			int minscore = 10000;
			string mingene = "";
			
			
			for (auto gene : v_target_list_sub_tmp)
			{
				if (gene == "") {continue;}
				pair <string,string> key = make_pair (read,gene);
				string nm = "-";
				
				
				if (long_score.find(key) != long_score.end()) {nm = to_string(long_score[key]);}
				ADD << "\t" << nm;
				
				if (nm != "-") {
					int currentscore = stoi(nm);
					if (minscore > currentscore){minscore = currentscore; mingene = gene; continue;}
					if (minscore == currentscore){mingene.append("," + gene); continue;}
				}
			}
			
			float ratio = float(minscore) / float(long_size[read]);
			
			if ((mingene == valid) && (ratio <= v_tolerance))
			{
				ADD << "\t" << valid << endl; 
				address_list_long[read] = valid;
			}
			else {ADD << "\t-" << endl;}
			
		}
 		 v_message = "Comparing scores: done";
		 screen_message (screen_size, 0, v_message, 1, v_quiet);
		 
	 }
	 

	 

	 
	 if (v_nanopore == 0) {
		 vector<string> v_target_list_sub;
		 boost::algorithm::erase_all(v_genes, " ");
		 string tmp = v_genes + ",others";
		 boost::split(v_target_list_sub,tmp,boost::is_any_of(","));
		 
			 
		float progress = 0;
		float totalreads = sequence_size_r1.size();
	 
		ThreadPool poolscore(stoi(v_threads));
		std::vector< std::future<int> > resultscore;
		int readcount = 0;
		vector <string> table;
		
		if (v_map_type == "paired") {
			for(map<string,int>::iterator it = sequence_size_r1.begin(); it != sequence_size_r1.end(); ++it)
			{
				string id = it->first;
				pair <string,string> key = make_pair ("others",id);
				int sub1 = 0;
				int sub2 = 0;
				if (sequence_list_r1.find(key) != sequence_list_r1.end()){sub1 = sequence_list_r1[key];}
				if (sequence_list_r2.find(key) != sequence_list_r2.end()){sub2 = sequence_list_r2[key];}
				if ((sub1 != 0) && (sub2 == 0)){sequence_list_r2[key] = sub1;}
				if ((sub2 != 0) && (sub1 == 0)){sequence_list_r1[key] = sub2;}
			}
		}
		
		 for(map<string,int>::iterator it = sequence_size_r1.begin(); it != sequence_size_r1.end(); ++it)
		 {
			 string id = it->first;
			 
			 resultscore.emplace_back(
			 poolscore.enqueue([readcount, id, v_target_list_sub, &address_list_r1, &table, &progress, totalreads]
			   {
			 
				 map <string,int> nm_sum;
				 string table_str = id;
			 
				 int min_nm_sum = 100000;
				 string address_to = "";
				 int count = 0;
			 
				 for( std::vector<string>::const_iterator g = v_target_list_sub.begin(); g != v_target_list_sub.end(); ++g)
				 {
					 if (*g == "") {continue;}
				 
					 nm_sum[*g] = 10000;
					 pair <string,string> key = make_pair (*g,id);
				 
					 string scores = "-";
				 
					 if (v_map_type == "paired"){
						 
						 if ((sequence_list_r1.find(key) != sequence_list_r1.end()) && (sequence_list_r2.find(key) != sequence_list_r2.end())){
							 if (((sequence_list_r1[key] != 10000) && (sequence_list_r2[key] != 10000)) && ((sequence_list_r1[key] != 0) && (sequence_list_r2[key] != 0)))
							 {
								 int v_max_r1 = 0;
								 int v_max_r2 = 0;
								 v_max_r1 = int((sequence_size_r1[id] * v_tolerance));
								 v_max_r2 = int((sequence_size_r2[id] * v_tolerance));
								 
	//                           if (*g == "others"){
	//                           if ((sequence_list_r1[key] == 1) && (sequence_list_r2[key] > v_max_r2))
	//                             {mtxg.lock();sequence_list_r2[key] = 1;mtxg.unlock();}
	//                             if ((sequence_list_r2[key] == 1) && (sequence_list_r1[key] > v_max_r1))
	//                             {mtxg.lock();sequence_list_r1[key] = 1;mtxg.unlock();}
	//                            }
								 
								 if (((sequence_list_r1[key]-1) <= v_max_r1) && ((sequence_list_r2[key]-1) <= v_max_r2))
								 {
									 scores = to_string(sequence_list_r1[key]) + "," + to_string(sequence_list_r2[key]);
									 if (min_nm_sum > (sequence_list_r1[key] + sequence_list_r2[key])) {min_nm_sum = (sequence_list_r1[key] + sequence_list_r2[key]);count++;}
									 nm_sum[*g] = sequence_list_r1[key] + sequence_list_r2[key];
								 }
							 }
						 }
					 }
					 
					 if (v_map_type == "single"){
						 if (sequence_list_r1.find(key) != sequence_list_r1.end()){
							 if ((sequence_list_r1[key] != 10000) && (sequence_list_r1[key] != 0))
							 {
								 int v_max_r1 = 0;
								 v_max_r1 = int((sequence_size_r1[id] * v_tolerance));
								 if ((sequence_list_r1[key]-1) <= v_max_r1)
								 {
									 scores = to_string(sequence_list_r1[key]);
									 if (min_nm_sum > sequence_list_r1[key] ) {min_nm_sum = (sequence_list_r1[key]);count++;}
									 nm_sum[*g] = sequence_list_r1[key];
								 }
							 }
						}
					 }
					
					 table_str.append ("\t" + scores);
			 }
			 
				 if (count >= 1)
				 {
					 int count_min_nm_sum = 0;
					 for( std::vector<string>::const_iterator g = v_target_list_sub.begin(); g != v_target_list_sub.end(); ++g)
					 {
						 if (nm_sum[*g] == min_nm_sum)
						 {
							 count_min_nm_sum++;
							 address_to = address_to + "," + *g;
							 pair <string,string> key = make_pair (*g,id);
							 mtxg.lock();
							 sequence_address[key] = 1;
							 mtxg.unlock();
						 }
					 }
					 
					 if (address_to != "") {
						 mtxg.lock();
						 address_list_r1[id] = address_to;
						 mtxg.unlock();
					 }
					 
				 }
				
				 
				 string subaddress = address_to;
				 if (subaddress == "") {subaddress = ",-";}
				 subaddress = subaddress.substr(1);
				 
				 table_str.append("\t" + subaddress);
				 table_str.append("\n");
				 mtxg.lock();
				 table.push_back(table_str);
				 progress++;
				 float ratio = ((progress / totalreads) * 100);
				 v_message = "Comparing scores ... " + to_string(ratio) + " % ";
				 screen_message (screen_size, 0, v_message, 2, v_quiet);
				 mtxg.unlock();
				 return 1;
			   })
			 );
		 }
		for(auto && result: resultscore){result.get();} // waiting for all threads

		
		
		string v_add = v_output + v_sample + "addressing_table.txt";
		ofstream ADD;
		ADD.open (v_add.c_str());
			
		ADD << "Read";
		for( std::vector<string>::const_iterator g = v_target_list_sub.begin(); g != v_target_list_sub.end(); ++g)
		{
		   if (*g == "") {continue;}
		   ADD << "\t" << *g;
		}
		ADD << "\tTarget" << endl;
		for (auto & item : table)
		{
			ADD << item;
		}
		ADD.close();
		 
		 v_message = "Comparing scores: done";
		 screen_message (screen_size, 0, v_message, 1, v_quiet);
		 
	  
		sequence_list_r1.clear();
		sequence_list_r2.clear();
		sequence_size_r1.clear();
		sequence_size_r2.clear();
    
	}


    
	
	
	
	
	
	// Creating fastq after filtering
	
	if (v_nanopore == 1) {
		for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
		{
		   if (*i == "") {continue;}

		   string gene = *i;
		   
		   v_message = "Making fastq files for " + *i + " ...";
		   screen_message (screen_size, 0, v_message, 2, v_quiet);
		   
		   string outfqtmp =  v_output + v_sample + *i + ".tmp.fq";
		   string outfqfinal = v_output + v_sample + *i + "_R0.fastq";
		   
			ofstream OUT_FINAL;
			OUT_FINAL.open (outfqfinal.c_str());
			ifstream FQ( outfqtmp.c_str());
			for( std::string line; getline( FQ, line ); )
			{
				string read = line;
				string subread = read.substr(1);
				getline( FQ, line );
				string seq = line;
				getline( FQ, line );
				string info = line;
				getline( FQ, line );
				string qual = line;
				
				if (address_list_long.find(subread) != address_list_long.end())
				{
					if (address_list_long[subread] == gene){
						OUT_FINAL << read << endl;
						OUT_FINAL << seq << endl;
						OUT_FINAL << info << endl;
						OUT_FINAL << qual << endl;
					}
				}
				
			}
			FQ.close();
			OUT_FINAL.close();
		   
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	

    
    
    
    
	
	
	
    
    
    if (v_nanopore == 0) {
    
		// Creating fastq after filtering
		unordered_map <string,string> original_data_r1;
		unordered_map <string,string> original_data_r2;
		unordered_map <string,string> original_data_r1_trim;
		unordered_map <string,string> original_data_r2_trim;
		vector <int> list_of_read_sizes;
		
		v_message = "Making fastq files ...";
		screen_message (screen_size, 0, v_message, 2, v_quiet);

		string r1_original = v_output + v_sample + "selected.R1.fastq";
		string r2_original = v_output + v_sample + "selected.R2.fastq";
		string r1_original_trim = v_output + v_sample + "selected.trim.R1.fastq";
		string r2_original_trim = v_output + v_sample + "selected.trim.R2.fastq";
		if (v_map_type == "single") {r1_original = v_output + v_sample + "selected.R0.fastq";}
		if (v_map_type == "single") {r1_original_trim = v_output + v_sample + "selected.trim.R0.fastq";}

	  
	  

		
		ifstream R1( r1_original.c_str());
		ifstream R2( r2_original.c_str());
		for( std::string line; getline( R1, line ); )
		{
			string id = line;
			getline( R1, line );
			string seq = line;
			getline( R1, line );
			string info = line;
			getline( R1, line );
			string qual = line;
			list_of_read_sizes.push_back(seq.size());
			vector<string> read_id;
			boost::split(read_id,id,boost::is_any_of(" "));
	//        if ((read_id[0].substr(read_id[0].size()-2,2) == "/1") || (read_id[0].substr(read_id[0].size()-2,2) == "/2")){read_id[0] = read_id[0].substr(0,read_id[0].size()-2);}

			
			original_data_r1[read_id[0]] = id + "\n" + seq + "\n" + info + "\n" + qual;
	 
			getline( R2, line );
			id = line;
			getline( R2, line );
			seq = line;
			getline( R2, line );
			info = line;
			getline( R2, line );
			qual = line;
			original_data_r2[read_id[0]] = id + "\n" + seq + "\n" + info + "\n" + qual;
		}
		R1.close();
		R2.close();
		
		int total_size = 0;
		int count_reads = 0;
		for (auto &&item: list_of_read_sizes)
		{
			total_size = total_size + item;
			count_reads++;
		}
		if (count_reads > 0) {read_mean_size = (total_size / count_reads);}
		if (count_reads == 0) {read_mean_size = 0;}
		list_of_read_sizes.clear();

		
		ifstream R1t( r1_original_trim.c_str());
		ifstream R2t( r2_original_trim.c_str());
		for( std::string line; getline( R1t, line ); )
		{
			string id = line;
			getline( R1t, line );
			string seq = line;
			getline( R1t, line );
			string info = line;
			getline( R1t, line );
			string qual = line;
			vector<string> read_id;
			boost::split(read_id,id,boost::is_any_of(" "));
	//        if ((read_id[0].substr(read_id[0].size()-2,2) == "/1") || (read_id[0].substr(read_id[0].size()-2,2) == "/2")){read_id[0] = read_id[0].substr(0,read_id[0].size()-2);}

			original_data_r1_trim[read_id[0]] = seq + "\n" + info + "\n" + qual;
			getline( R2t, line );
			id = line;
			getline( R2t, line );
			seq = line;
			getline( R2t, line );
			info = line;
			getline( R2t, line );
			qual = line;
			original_data_r2_trim[read_id[0]] = seq + "\n" + info + "\n" + qual;
		 }
		 R1t.close();
		 R2t.close();
		


		if (v_map_type == "paired") {
		   for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
		   {
			   if (*i == "") {continue;}

			   string gene = *i;
			   
			   v_message = "Making fastq files for " + *i + " ...";
			   screen_message (screen_size, 0, v_message, 2, v_quiet);
			   
			   string sort = v_output + v_sample + "sorted_" + *i + "_R1.fastq";
			   
			   map <string,int> ids;
			   ifstream idsrec( sort.c_str());
			   for( std::string line; getline( idsrec, line ); )
			   {
				   string id = line;
				   getline( R1t, line );
				   string seq = line;
				   getline( R1t, line );
				   string info = line;
				   getline( R1t, line );
				   string qual = line;
				   
				   vector <string> read_id;
				   boost::split(read_id,id,boost::is_any_of(" "));
	//               if ((read_id[0].substr(read_id[0].size()-2,2) == "/1") || (read_id[0].substr(read_id[0].size()-2,2) == "/2")){read_id[0] = read_id[0].substr(0,read_id[0].size()-2);}
				   ids[read_id[0]] = 1;
			   }
			   idsrec.close();
	  
			   int size = gene_opt_end[gene] - gene_opt_start[gene];
			   int downS = (downsampling * size) / read_mean_size;
			   int downS_counter = 0;
			   
			   string v_out_R1 = v_output + v_sample + *i + "_R1.fastq";
			   ofstream OUT_R1;
			   OUT_R1.open (v_out_R1.c_str());
			   
			   string v_out_R2 = v_output + v_sample + *i + "_R2.fastq";
			   ofstream OUT_R2;
			   OUT_R2.open (v_out_R2.c_str());

			   string v_out_clean_R1 = v_output + v_sample + *i + ".trim.R1.fastq";
			   ofstream OUT_clean_R1;
			   OUT_clean_R1.open (v_out_clean_R1.c_str());

			   string v_out_clean_R2 = v_output + v_sample + *i + ".trim.R2.fastq";
			   ofstream OUT_clean_R2;
			   OUT_clean_R2.open (v_out_clean_R2.c_str());

			   
			   for (auto & id : ids)
			   {
				   
				   vector <string> data_split_r1;
				   boost::split(data_split_r1,original_data_r1[id.first],boost::is_any_of("\n"));

				   pair <string,string> key = make_pair (gene,id.first.substr(1));
				   if (sequence_address.find(key) == sequence_address.end()) {continue;}
	   
				   vector <string> data_split_r1_trim;
				   boost::split(data_split_r1_trim,original_data_r1_trim[id.first],boost::is_any_of("\n"));
	   
					vector <string> data_split_r2;
					boost::split(data_split_r2,original_data_r2[id.first],boost::is_any_of("\n"));
					
					vector <string> data_split_r2_trim;
					boost::split(data_split_r2_trim,original_data_r2_trim[id.first],boost::is_any_of("\n"));
				   
				   

				   OUT_R1 << data_split_r1[0] << endl;
				   OUT_R1 << data_split_r1[1] << endl;
				   OUT_R1 << data_split_r1[2] << endl;
				   OUT_R1 << data_split_r1[3] << endl;
				   OUT_R2 << data_split_r2[0] << endl;
				   OUT_R2 << data_split_r2[1] << endl;
				   OUT_R2 << data_split_r2[2] << endl;
				   OUT_R2 << data_split_r2[3] << endl;

				   if (downS_counter < downS) {
							   OUT_clean_R1 << id.first + " 1:trim" << endl;
							   OUT_clean_R1 << data_split_r1_trim[0] << endl;
							   OUT_clean_R1 << data_split_r1_trim[1] << endl;
							   OUT_clean_R1 << data_split_r1_trim[2] << endl;
							   OUT_clean_R2 << id.first + " 2:trim" << endl;
							   OUT_clean_R2 << data_split_r2_trim[0] << endl;
							   OUT_clean_R2 << data_split_r2_trim[1] << endl;
							   OUT_clean_R2 << data_split_r2_trim[2] << endl;
				   }
					downS_counter++;
					continue;
			   }
			   OUT_R1.close();
			   OUT_R2.close();
			   OUT_clean_R1.close();
			   OUT_clean_R2.close();
		   }
		}

		
		if (v_map_type == "single") {
			for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
			{
				if (*i == "") {continue;}

				string gene = *i;
				
				v_message = "Making fastq files for " + *i + " ...";
				screen_message (screen_size, 0, v_message, 2, v_quiet);

				string sort = v_output + v_sample + "sorted_" + *i + "_R0.fastq";
				
				map <string,int> ids;
				ifstream idsrec( sort.c_str());
				for( std::string line; getline( idsrec, line ); )
				{
					string id = line;
					getline( R1t, line );
					string seq = line;
					getline( R1t, line );
					string info = line;
					getline( R1t, line );
					string qual = line;
					vector <string> read_id;
					boost::split(read_id,id,boost::is_any_of(" "));
					ids[read_id[0]] = 1;
				}
				idsrec.close();
				
				
				
				

				int size = gene_opt_end[gene] - gene_opt_start[gene];
				int downS = (downsampling * size) / read_mean_size;
				int downS_counter = 0;

				
				string v_out_R1 = v_output + v_sample + *i + "_R0.fastq";
				ofstream OUT_R1;
				OUT_R1.open (v_out_R1.c_str());
				
				string v_out_singlet = v_output + v_sample + *i + ".trim.R0.fastq";
				ofstream OUT_singlet;
				OUT_singlet.open (v_out_singlet.c_str());
				
				for (auto & id : ids)
				{
					vector <string> data_split_r1;
					boost::split(data_split_r1,original_data_r1[id.first],boost::is_any_of("\n"));
					
					pair <string,string> key = make_pair (gene,id.first.substr(1));
					if (sequence_address.find(key) == sequence_address.end()) {continue;}
		
					vector <string> data_split_r1_trim;
					boost::split(data_split_r1_trim,original_data_r1_trim[id.first],boost::is_any_of("\n"));
	  
	 
					OUT_R1 << data_split_r1[0] << endl;
					OUT_R1 << data_split_r1[1] << endl;
					OUT_R1 << data_split_r1[2] << endl;
					OUT_R1 << data_split_r1[3] << endl;

					if (downS_counter < downS) {
						OUT_singlet << id.first << " 1:trim" << endl;
						OUT_singlet << data_split_r1_trim[0] << endl;
						OUT_singlet << data_split_r1_trim[1] << endl;
						OUT_singlet << data_split_r1_trim[2] << endl;
					}
					downS_counter++;
					continue;
				}
				OUT_R1.close();
				OUT_singlet.close();
			}
		}
		
	/*    
		original_data_r2.clear();
		original_data_r1.clear();
		original_data_r1_trim.clear();
		original_data_r2_trim.clear();
	*/
	   
		v_message = "Making fastq files: done";
		screen_message (screen_size, 0, v_message, 1, v_quiet);
	}

    

 
    
    

   
     
    
    
    
    
    
    //Adjustments
    
    
    int mincov = 0;
    for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
    {
        if (*i == "KIR2DL4") {mincov = 5;}
    }

    map <string,string> typing_result_alleleA;
    map <string,string> typing_result_alleleB;

    
    if (v_rna == 1) {v_skiptyping = 1;}

    
    if (v_skiptyping == 0) {
        v_message = "Performing adjustments ...";
        screen_message (screen_size, 0, v_message, 2, v_quiet);
        string v_tag_sub = v_sample.substr(0, v_sample.size()-1);
        
        string cmd = "mkdir " + v_output + "/adjustment/";
        GetStdoutFromCommand(cmd.c_str());
        
        ThreadPool adj(stoi(v_threads));
        std::vector< std::future<int> > results_adj;
        
        map <string,string> typing_result_alleleA;
        map <string,string> typing_result_alleleB;
        vector <string> gen_log;
        int count = 0;
        
        for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
        {
            if (*i == "") {continue;}
            
            if (gene_type.find(*i) != gene_type.end()) {
                if (gene_type[*i]  != "TYPE") {continue;}
            }
            
            if (locus_presence.find(*i) != locus_presence.end())
            {
                if (locus_presence[*i] == 0) {continue;}
            }
            string gene = *i;
            
            v_message = "Adjusting " + *i + " ...";
            screen_message (screen_size, 0, v_message, 2, v_quiet);

            
            if (v_typeall == 0)
            {
                if (gene_type[gene] != "TYPE")
                {
                    gen_log.push_back(gene + ": alleles used for adjustments - disabled");
                    continue;
                }
            }

            
                string cmd = "mkdir " + v_output + "/adjustment/" + gene;
                GetStdoutFromCommand(cmd.c_str());
                
                string outadj = v_output + "/adjustment/" + gene;
                
                string fq1 = "";
                string fq2 = "";
                if (v_map_type == "paired")
                {
                    fq1 = v_output + v_sample + gene + "_R1.fastq";
                    fq2 = v_output + v_sample + gene + "_R2.fastq";
                }
                if (v_map_type == "single")
                {
                    fq1 = v_output + v_sample + gene + "_R0.fastq";
                    fq2 = "";
                }
                

            string adjtype = "cds";
            if (v_exome == 0){adjtype = "full";}
            string alleles_adjustment = adjust_freeb(gene,fq1,fq2,outadj, v_threads,adjtype,v_db);
            
            if (alleles_adjustment.find("ailed") != string::npos)
            {
                gen_log.push_back( gene + ": alleles used for adjustments - failed");
                typing_result_alleleA[gene] = "";
                typing_result_alleleB[gene] = "";
                continue;
            }
                
                vector <string> alleles;
                boost::split(alleles,alleles_adjustment,boost::is_any_of("\t"));
                
                typing_result_alleleA[gene] = alleles[0];
                typing_result_alleleB[gene] = alleles[1];
                gen_log.push_back( gene + ": alleles used for adjustments: " + alleles[0] + " and " + alleles[1]);
                
                count++;
                continue;
            
        }
        
        for (auto item : gen_log)
        {
            general_log << item << endl;
        }
        
        v_message = "Performing adjustments: done";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        
    }

    
    for (auto & item : reads_possibly_mapped)
    {
        string gene = item.first.first;
        string read = item.first.second;
        int nm = item.second;
        
        if (reads_mapped_nm.find(read) == reads_mapped_nm.end())
        {
            reads_mapped_nm[read] = nm;
            reads_mapped_gene[read] = gene;
        }
        else
        {
            if (reads_mapped_nm[read] == nm)
            {
                reads_mapped_gene[read].append(";" + gene);
            }
            if (reads_mapped_nm[read] > nm)
            {
                reads_mapped_nm[read] = nm;
                reads_mapped_gene[read] = gene;
            }
        }
    }
    
    

    








    
    
    // Final mapping
	 if (v_nanopore == 1) {
	
		for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
		{
			
			if (*i == "") {continue;}
	
			v_message = "Mapping " + *i + " ...";
			screen_message (screen_size, 0, v_message, 2, v_quiet);
			
			string v_tag_sub = v_sample.substr(0, v_sample.size()-1);
			string v_reference = "";
			v_reference = "'" + v_db + "/reference/dna/loci/" + *i + ".fas'";
			boost::replace_all(v_reference, "\\ ", " ");

			string in = v_output + v_sample + *i + "_R0.fastq";
			string v_sam1 = v_output + v_sample + *i + ".tmp.sam";
			string v_sam2 = v_output + v_sample + *i + ".unique.sam";
			string v_log = v_output + "log/" + v_sample + *i + ".map.log";
			
			string readgroup = "'@RG\\tID:" + v_sample_sub + "\\tSM:" + v_sample_sub + "'";
			
			v_command = v_minimap + " -a -t " + v_threads + " -x map-ont -A 1 -B 2 -O 2,12 -R " + readgroup + " -o " + v_sam1 + " " + v_reference + " " + in;
			string out = GetStdoutFromCommand (v_command.c_str());
			
			ofstream logmap;
			logmap.open (v_log.c_str());
			logmap << out << endl;
			logmap.close();
			
			
			ifstream reference (v_reference.c_str());
			string refseq = "";
			for( std::string line; getline( reference, line ); )
			{
				if (line.substr(0,1) == ">") {continue;}
				if (line == "") {continue;}
				refseq = refseq + line;
			}
			reference.close();
			ref_size[*i] = refseq.length();
			refseq = "";
			
				 
			fstream sam (v_sam1.c_str());
			ofstream OUT1;
			OUT1.open (v_sam2.c_str());
		
				 
			int count_reads = 0;
			int v_pos_correct = position_hg38[*i] - 1;
			
			for( std::string line; getline( sam, line ); )
			{
				if (line == "") {continue;}
				vector<string> sam_line;
				boost::split(sam_line,line,boost::is_any_of("\t"));
				 
				if (sam_line[0] == "@SQ")
				{
					if (usechr == 0) {
						OUT1 << "@SQ\tSN:" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;
					}

					if (usechr == 1) {
						size_t found = chr_hg38[*i].find("chr");
						if (found==std::string::npos){
							OUT1 << "@SQ\tSN:" << "chr" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;
						}
						else {
							OUT1 << "@SQ\tSN:" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;
						}
					}

					OUT1 << "@CO\tkir-mapper " << Program_version << ", human genome version hg38 " << endl;
					continue;
				}
				 
				if ((sam_line[0].substr(0,1) == "@") || (sam_line[2] == "*"))
				{
					OUT1 << line << endl;
					continue;
				}


				string seq = sam_line[9];
				if (seq.length() < v_size) {continue;}

				sam_line[2] = chr_hg38[*i];
				if (usechr == 1)
				{
					size_t found = chr_hg38[*i].find("chr");
					if (found==std::string::npos){sam_line[2] = "chr" + chr_hg38[*i];}
				}
			 
				string tmp = address_list_r1[sam_line[0]] + ",";
				vector<string> address_multi_hits;
				boost::split(address_multi_hits,tmp,boost::is_any_of(","));
				 
				OUT1 << sam_line[0] << "\t" << sam_line[1] << "\t" << sam_line[2] << "\t";
				OUT1 << (stoi(sam_line[3]) + v_pos_correct) << "\t";
				OUT1 << sam_line[4] << "\t" << sam_line[5] << "\t" << sam_line[6] << "\t";
				OUT1 << (stoi(sam_line[7]) + v_pos_correct);
				for (int helper = 8; helper < sam_line.size(); helper++)
				{OUT1 << "\t" << sam_line[helper];}
				OUT1 << endl;
				count_reads++;
				continue;
			}
			sam.close();
			OUT1.close();
		}
		v_message = "Mapping: done";
		screen_message (screen_size, 0, v_message, 1, v_quiet);
	 }
	 
	 
	
    if (v_nanopore == 0) {
		for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
		{
			
			if (*i == "") {continue;}
	/*
			if (locus_presence.find(*i) != locus_presence.end())
			{
				if (locus_presence[*i] == 0) {continue;}
			}
	*/
			v_message = "Mapping " + *i + " ...";
			screen_message (screen_size, 0, v_message, 2, v_quiet);
			
			string v_tag_sub = v_sample.substr(0, v_sample.size()-1);
			string v_reference = "";
			v_reference = "'" + v_db + "/reference/dna/loci/" + *i + ".fas'";
			boost::replace_all(v_reference, "\\ ", " ");

			string in0 = v_output + v_sample + *i + "_R0.fastq";
			string in1 = v_output + v_sample + *i + "_R1.fastq";
			string in2 = v_output + v_sample + *i + "_R2.fastq";
			string v_sam1 = v_output + v_sample + *i + ".tmp.sam";
			string v_sam2 = v_output + v_sample + *i + ".unique.sam";
			string v_sam3 = v_output + v_sample + *i + ".adjusted.sam";
			string v_log = v_output + "log/" + v_sample + *i + ".map.log";
				 
			if ((v_rna == 0) && (v_nanopore == 0)) {     
				if (v_map_type == "paired") {
					v_command = v_bwa + " mem -t " + v_threads + " -B 2 -O 3,3 -L 3,3 -R '@RG\\tID:" + v_tag_sub + "\\tLB:" + v_tag_sub + "\\tSM:" + v_tag_sub + "\\tPL:illumina\\tPU:" + v_tag_sub + "' " + v_reference + " " + in1 + " " + in2 + " > " + v_sam1 + " 2>" + v_log;
					if (*i == "HLA-DRB1")
					{
						v_command = v_bwa + " mem -t " + v_threads + " -k 15 -B 2 -O 3,3 -L 3,3 -R '@RG\\tID:" + v_tag_sub + "\\tLB:" + v_tag_sub + "\\tSM:" + v_tag_sub + "\\tPL:illumina\\tPU:" + v_tag_sub + "' " + v_reference + " " + in1 + " " + in2 + " > " + v_sam1 + " 2>" + v_log;
					}
					
				}
				if (v_map_type == "single") {
					v_command = v_bwa + " mem -t " + v_threads + " -B 2 -O 3,3 -L 3,3 -R '@RG\\tID:" + v_tag_sub + "\\tLB:" + v_tag_sub + "\\tSM:" + v_tag_sub + "\\tPL:illumina\\tPU:" + v_tag_sub + "' " + v_reference + " " + in0 + " > " + v_sam1 + " 2>" + v_log;
				}
				system (v_command.c_str());
			}
			

			 if (v_rna == 1) {     
			 
				string refstar = v_db + "/reference/rna/" + *i;
				v_sam1 = v_output + v_sample + *i + ".";
				if (v_map_type == "single") {
					v_command = v_star + " --runThreadN " + v_threads + " --genomeDir " + refstar + " --readFilesIn " + in0 + " --outFileNamePrefix " + v_sam1 + " --outFilterMismatchNoverReadLmax 0.08 --alignIntronMax 3000  2>" + v_log;
				}
				if (v_map_type == "paired") {
					v_command = v_star + " --runThreadN " + v_threads + " --genomeDir " + refstar + " --readFilesIn " + in1 + " " + in2 + " --outFileNamePrefix " + v_sam1 + " --outFilterMismatchNoverReadLmax 0.08 --alignIntronMax 3000 2>" + v_log;
				}
				string rg_string = "ID:" + v_sample_sub + " LB:" + v_sample_sub + " SM:" + v_sample_sub;
				v_command.append(" --outSAMattrRGline " + rg_string + " --outSAMunmapped Within");
				debug_message(v_command);
				GetStdoutFromCommand(v_command);
				 v_sam1 = v_output + v_sample + *i + ".Aligned.out.sam";
				 string tmpdir = v_output + v_sample + *i + "._STARtmp";
				 GetStdoutFromCommand("rm -rf " + tmpdir);
			 }


			ifstream reference (v_reference.c_str());
			string refseq = "";
			for( std::string line; getline( reference, line ); )
			{
				if (line.substr(0,1) == ">") {continue;}
				if (line == "") {continue;}
				refseq = refseq + line;
			}
			reference.close();
			ref_size[*i] = refseq.length();
			refseq = "";
				 
				 
			fstream sam (v_sam1.c_str());
		
			ofstream OUT1;
			OUT1.open (v_sam2.c_str());
			ofstream OUT2;
			OUT2.open (v_sam3.c_str());
		
				 
			int count_reads = 0;
			int v_pos_correct = position_hg38[*i] - 1;
			unordered_map <string,int> adjusted_for_secondary;
			unordered_map <string,int> adjusted_for_primary;
				 
			for( std::string line; getline( sam, line ); )
			{
				if (line == "") {continue;}
				vector<string> sam_line;
				boost::split(sam_line,line,boost::is_any_of("\t"));
				 
				if (sam_line[0] == "@SQ")
				{
					if (usechr == 0) {
						OUT1 << "@SQ\tSN:" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;
					}

					if (usechr == 1) {
						size_t found = chr_hg38[*i].find("chr");
						if (found==std::string::npos){
							OUT1 << "@SQ\tSN:" << "chr" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;
						}
						else {
							OUT1 << "@SQ\tSN:" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;
						}
					}

					OUT1 << "@CO\tkir-mapper " << Program_version << ", human genome version hg38 " << endl;
					
					if (usechr == 0) {
						OUT2 << "@SQ\tSN:" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;
					}
					if (usechr == 1) {
						size_t found = chr_hg38[*i].find("chr");
						if (found==std::string::npos){OUT2 << "@SQ\tSN:" << "chr" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;}
						else {OUT2 << "@SQ\tSN:" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;}
					}
					OUT2 << "@CO\tkir-mapper " << Program_version << ", human genome version hg38 " << endl;
					continue;
				}
				 
				if ((sam_line[0].substr(0,1) == "@") || (sam_line[2] == "*"))
				{
					OUT1 << line << endl;
					OUT2 << line << endl;
					continue;
				}


				string seq = sam_line[9];
				if (seq.length() < v_size) {continue;}

				sam_line[2] = chr_hg38[*i];
				if (usechr == 1)
				{
					size_t found = chr_hg38[*i].find("chr");
					if (found==std::string::npos){sam_line[2] = "chr" + chr_hg38[*i];}
				}
			 
				string tmp = address_list_r1[sam_line[0]] + ",";
				vector<string> address_multi_hits;
				boost::split(address_multi_hits,tmp,boost::is_any_of(","));
				 

				if (address_multi_hits.size() < 4)
				{
					OUT1 << sam_line[0] << "\t" << sam_line[1] << "\t" << sam_line[2] << "\t";
					OUT1 << (stoi(sam_line[3]) + v_pos_correct) << "\t";
					OUT1 << sam_line[4] << "\t" << sam_line[5] << "\t" << sam_line[6] << "\t";
					OUT1 << (stoi(sam_line[7]) + v_pos_correct);
					for (int helper = 8; helper < sam_line.size(); helper++)
					{OUT1 << "\t" << sam_line[helper];}
					OUT1 << endl;
					OUT2 << sam_line[0] << "\t" << sam_line[1] << "\t" << sam_line[2] << "\t";
					OUT2 << (stoi(sam_line[3]) + v_pos_correct) << "\t";
					OUT2 << sam_line[4] << "\t" << sam_line[5] << "\t" << sam_line[6] << "\t";
					OUT2 << (stoi(sam_line[7]) + v_pos_correct);
					for (int helper = 8; helper < sam_line.size(); helper++)
					{OUT2 << "\t" << sam_line[helper];}
					OUT2 << endl;
					count_reads++;
					continue;
				}
		
				
				if (address_multi_hits.size() >= 4)
				{
					int adjust = 0;

					pair <string,string> key = make_pair(*i,sam_line[0]);
					if (reads_mapped_gene.find(sam_line[0]) != reads_mapped_gene.end())
					{
						
						if (reads_mapped_gene[sam_line[0]].find(*i) != std::string::npos)
						{
							vector <string> sub;
							string subtext = reads_mapped_gene[sam_line[0]].substr(1);
							boost::split(sub,subtext,boost::is_any_of(";"));
							int v1 = rand() % sub.size();
							if (sub[v1] == *i) {
								adjust = 1;
							}
						}
						if (reads_mapped_gene[sam_line[0]] == *i)
						{
							adjust = 1;
						}
					}

					if (adjust == 1) {
							OUT2 << sam_line[0] << "\t" << sam_line[1] << "\t" << sam_line[2] << "\t";
							OUT2 << (stoi(sam_line[3]) + v_pos_correct) << "\t";
							OUT2 << sam_line[4] << "\t" << sam_line[5] << "\t" << sam_line[6] << "\t";
							OUT2 << (stoi(sam_line[7]) + v_pos_correct);
							for (int helper = 8; helper < sam_line.size(); helper++)
							{OUT2 << "\t" << sam_line[helper];}
							OUT2 << endl;
							count_reads++;
							adjusted_for_primary[sam_line[0]] = 1;
					}
						 
						// for paired
						if (sam_line[1] == "163") {sam_line[1] = "419";}
						if (sam_line[1] == "83") {sam_line[1] = "339";}
						if (sam_line[1] == "99") {sam_line[1] = "355";}
						if (sam_line[1] == "147") {sam_line[1] = "403";}
						if (sam_line[1] == "81") {sam_line[1] = "337";}
						if (sam_line[1] == "167") {sam_line[1] = "423";}
						if (sam_line[1] == "145") {sam_line[1] = "401";}
						if (sam_line[1] == "97") {sam_line[1] = "353";}
						if (sam_line[1] == "185") {sam_line[1] = "441";}
						if (sam_line[1] == "73") {sam_line[1] = "329";}
						if (sam_line[1] == "113") {sam_line[1] = "369";}
						if (sam_line[1] == "117") {sam_line[1] = "373";}
						if (sam_line[1] == "121") {sam_line[1] = "377";}
						if (sam_line[1] == "133") {sam_line[1] = "389";}
						if (sam_line[1] == "137") {sam_line[1] = "393";}
						if (sam_line[1] == "177") {sam_line[1] = "433";}
						if (sam_line[1] == "181") {sam_line[1] = "437";}
						if (sam_line[1] == "161") {sam_line[1] = "417";}

						// for single
						if (sam_line[1] == "16") {sam_line[1] = "272";}
						if (sam_line[1] == "0") {sam_line[1] = "256";}

						 
						adjusted_for_secondary[sam_line[0]] = 1;
						 
						OUT1 << sam_line[0] << "\t" << sam_line[1] << "\t" << sam_line[2] << "\t";
						OUT1 << (stoi(sam_line[3]) + v_pos_correct) << "\t";
						OUT1 << sam_line[4] << "\t" << sam_line[5] << "\t" << sam_line[6] << "\t";
						OUT1 << (stoi(sam_line[7]) + v_pos_correct);
						for (int helper = 8; helper < sam_line.size(); helper++)
						{OUT1 << "\t" << sam_line[helper];}
						OUT1 << endl;
						 
						if (adjust == 0) {
							   OUT2 << sam_line[0] << "\t" << sam_line[1] << "\t" << sam_line[2] << "\t";
							   OUT2 << (stoi(sam_line[3]) + v_pos_correct) << "\t";
							   OUT2 << sam_line[4] << "\t" << sam_line[5] << "\t" << sam_line[6] << "\t";
							   OUT2 << (stoi(sam_line[7]) + v_pos_correct);
							   for (int helper = 8; helper < sam_line.size(); helper++)
							   {OUT2 << "\t" << sam_line[helper];}
							   OUT2 << endl;
						   }
					 
				}
			}
				sam.close();
				OUT1.close();
				OUT2.close();

				string log = v_output + v_sample + *i + ".log";
				ofstream log_mapping;
				log_mapping.open (log.c_str());

				log_mapping << endl << "List of reads marked as secondary" << endl;
				for (auto & item : adjusted_for_secondary)
				{
					 log_mapping << item.first << endl;
				}
				log_mapping << endl << endl;
				
				if (typing_result_alleleA[*i] != "") {
					log_mapping << *i << ": alleles considered for adjustments: " << typing_result_alleleA[*i] << ", " << typing_result_alleleB[*i] << endl;
				}
			
				log_mapping << endl;
				log_mapping << "List of reads adjusted for primary:" << endl;
				for (auto & item : adjusted_for_primary)
				{
					 log_mapping << item.first << endl;
				}
				log_mapping << endl;
				 
				 
				log_mapping << "Number of mapped reads (only primary mappings): " << to_string(count_reads) << endl;
				 

	//            double cov = ((double(count_reads) * read_mean_size) / ref_size[*i]);
			
	//           log_mapping << "Depth: " << to_string(cov) << endl;
	//            general_log << *i << " depth: " << to_string(cov) << endl;
				read_count_gene[*i] = double(count_reads);
				 
				log_mapping.close();
			
		}
		v_message = "Mapping: done";
		screen_message (screen_size, 0, v_message, 1, v_quiet);
	}
    
    
    
    
    reads_mapped_nm.clear();
    reads_mapped_gene.clear();
    read_count_gene.clear();
    reads_possibly_mapped.clear();
    sequence_list_r1.clear();
    sequence_list_r2.clear();
    sequence_size_r1.clear();
    sequence_size_r2.clear();
    sequence_address.clear();
    read_count_gene.clear();

 

    
    
    v_message = "Merging files ...";
    screen_message (screen_size, 0, v_message, 2, v_quiet);
    
    map <string,int> sq_data;
    map <string,int> rg_data;
    unordered_map <string,int> optimized_reads;
    vector <string> adjusted_data;
    vector <string> unique_data;
    vector <string> hg38_data;
    
    
    v_message = "Merging files ... reading new alignments ...";
    screen_message (screen_size, 0, v_message, 2, v_quiet);

    thread t1 ([&sq_data,&optimized_reads,&adjusted_data,&rg_data]{
		if (v_nanopore == 1) {return;}
        for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
        {
            if (*i ==  "") {continue;}
            string file = v_output + v_sample + *i + ".adjusted.sam";
            ifstream input( file.c_str() );
            for( std::string line; getline( input, line ); )
            {
                if (line == "") {continue;}
                if (line.substr(0,3) == "@PG") {continue;}
                if (line.substr(0,3) == "@SQ") {
                    mtxg.lock();sq_data[line] = 1;mtxg.unlock();
                    continue;
                }
                if (line.substr(0,3) == "@RG") {mtxg.lock();rg_data[line] = 1;mtxg.unlock();continue;}
                if (line.substr(0,1) == "@") {continue;}
                vector <string> data;
                boost::split(data,line,boost::is_any_of("\t"));
                if (v_map_type != "single") {
                    if ((data[8] == "0") || (data[6] == "*")) {continue;}
                }
//                if (usechr == 1)
//                {
                    if (data[2] == "6") {data[2] = "chr6";}
                    if (data[2] == "19") {data[2] = "chr19";}
                    line = data[0];
                    for (int a = 1; a < data.size(); a++)
                    {
                        line.append("\t" + data[a]);
                    }
 //               }
                
                mtxg.lock();
                optimized_reads[data[0]] = 1;
                mtxg.unlock();
                adjusted_data.push_back(line);
            }
            input.close();
        }
    });
    
    
    thread t2 ([&sq_data,&optimized_reads,&unique_data,&rg_data]{

        for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
        {
            if (*i ==  "") {continue;}
            string file = v_output + v_sample + *i + ".unique.sam";
            ifstream input( file.c_str() );
            for( std::string line; getline( input, line ); )
            {
                if (line == "") {continue;}
                if (line.substr(0,3) == "@PG") {continue;}
                if (line.substr(0,3) == "@SQ")
                {
                    mtxg.lock();sq_data[line] = 1;mtxg.unlock();
                    continue;
                }
                if (line.substr(0,3) == "@RG") {mtxg.lock();rg_data[line] = 1;mtxg.unlock();continue;}
                if (line.substr(0,1) == "@") {continue;}
                vector <string> data;
                boost::split(data,line,boost::is_any_of("\t"));
                if (v_map_type != "single") {
                    if ((data[8] == "0") || (data[6] == "*")) {continue;}
                }
//                if (usechr == 1)
//                {
                    if (data[2] == "6") {data[2] = "chr6";}
                    if (data[2] == "19") {data[2] = "chr19";}
                    line = data[0];
                    for (int a = 1; a < data.size(); a++)
                    {
                        line.append("\t" + data[a]);
                    }
//                }
                
                mtxg.lock();
                optimized_reads[data[0]] = 1;
                mtxg.unlock();
                unique_data.push_back(line);
            }
            input.close();
        }
     });
    
    t1.join();
    t2.join();
    
    
 
    
    
    if ((v_bam != "") && (v_rna == 0))
    {
        general_log << endl << "Reads forced to MQ=0 after optimization:" << endl;
        v_message = "Merging files ... reading original alignments ...";
        screen_message (screen_size, 0, v_message, 2, v_quiet);
        string hg38_select = v_output + v_sample + "hg38.sam";

        string v_bed = "'" + v_db + "/bed/target_dna.bed'";
        boost::replace_all(v_bed, "\\ ", " ");

        v_command = v_samtools + " view -h -@ " + v_threads + " -ML " + v_bed + " " + v_bam;
        v_system_out = GetStdoutFromCommand(v_command);
        
        vector <string> hg38_full;
        boost::split(hg38_full,v_system_out,boost::is_any_of("\n"));
        
        for( auto &line : hg38_full)
        {
            if (line == "") {continue;}
            if (line.substr(0,1) == "[") {continue;}
            if (line.substr(0,3) == "@PG") {continue;}
            if (line.substr(0,3) == "@SQ") {
                continue;
            }
            if (line.substr(0,3) == "@RG") {continue;}
            if (line.substr(0,1) == "@") {continue;}
            vector <string> data;
            boost::split(data,line,boost::is_any_of("\t"));

            
 
            
            if (data[1] == "*") {continue;}
            if (data[2] == "*") {continue;}
            if (data[3] == "0") {continue;}
            if (data[3] == "*") {continue;}
            if (data[4] == "0") {continue;}
            if (data[6] == "*") {continue;}
            if (data[8] == "0") {continue;}
            if (data[8] == "*") {continue;}
            if (data[6] != "=") {continue;}
            
            
            if (optimized_reads.find(data[0]) == optimized_reads.end())
            {
                
                for (int a = 11; a < data.size(); a++)
                {
                    if (data[a].substr(0,5) == "RG:Z:") {data[a] = "RG:Z:" + v_sample_sub;break;}
                }
                
                int pos = stoi(data[3]);
                int helper = 0;
                for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
                {
                    if (*i ==  "") {continue;}
                    if ((pos >= gene_opt_start[*i]) && (pos <= gene_opt_end[*i]))
                    {
                        helper = 1;
                        break;
                    }
                }
                if (helper == 1){
                    data[4] = "0";
                    general_log << data[0] << endl;
                }
                
                string news = data[0];
                for (int a = 1; a < data.size(); a++)
                {news.append("\t" + data[a]);}
                hg38_data.push_back(news);
                continue;
            }
        }
     }
    
    
    
    
    
    v_message = "Merging files ... writing merged files ...";
    screen_message (screen_size, 0, v_message, 2, v_quiet);

	if (v_nanopore == 0) {
		string merged_sam_adjusted = v_output + v_sample_sub + ".adjusted.sam";
		ofstream sam_adjusted;
		sam_adjusted.open (merged_sam_adjusted .c_str());
		sam_adjusted <<  "@CO\tkir-mapper " << Program_version << ", human genome version hg38" << endl;
		
		for (auto &item : sq_data) {sam_adjusted << item.first << endl;}
		for (auto &item : rg_data) {sam_adjusted << item.first << endl;}
		for (auto &item : adjusted_data) {sam_adjusted << item << endl;}
		for (auto &item : hg38_data) {sam_adjusted << item << endl;}
		sam_adjusted.close();

		string v_bam_out = v_output + v_sample_sub + ".adjusted.bam";
		string v_log = v_output + "/log/adjusted_sort.log";
		v_command = v_samtools + " sort -@ " + v_threads + " -m 1g " + merged_sam_adjusted  + " > " +  v_bam_out + " 2>" + v_log;
		system (v_command.c_str());
		v_log = v_output + "/log/adjusted_index.log";
		v_command = v_samtools + " index " + v_bam_out + " 2>" + v_log;
		system (v_command.c_str());
		
		string v_bam_nodup = v_output + v_sample_sub + ".adjusted.nodup.bam";;
		string v_bam_nodup_metrics = v_output + v_sample_sub + ".adjusted.nodup.metrics";
		
		
		
		v_log = v_output + "/log/adjusted_nodup.log";

		if (v_nanopore == 0) {
			std::size_t found = v_picard.find(".jar");
			if (found!=std::string::npos)
			{
				v_command = "java -jar " + v_picard + " MarkDuplicates I=" +  v_bam_out + " O=" + v_bam_nodup + " M=" + v_bam_nodup_metrics + " VALIDATION_STRINGENCY=SILENT 2>" + v_log;
			}
			else {
				v_command = v_picard + " MarkDuplicates I=" +  v_bam_out + " O=" + v_bam_nodup + " M=" + v_bam_nodup_metrics + " VALIDATION_STRINGENCY=SILENT 2>" + v_log;
			}

			if (v_rna == 0) {
				v_log = v_output + "/log/adjusted_nodup_index.log";
				system (v_command.c_str());
				v_command = v_samtools + " index " + v_bam_nodup + " 2>" + v_log;
				system (v_command.c_str());
			}
		}
		
    }
	
    string merged_sam_unique = v_output + v_sample_sub + ".unique.sam";
    ofstream sam_unique;
    sam_unique.open (merged_sam_unique .c_str());
    sam_unique <<  "@CO\tkir-mapper " << Program_version << ", human genome version hg38" << endl;
    for (auto &item : sq_data) {sam_unique << item.first << endl;}
    for (auto &item : rg_data) {sam_unique << item.first << endl;}
    for (auto &item : unique_data) {sam_unique << item << endl;}
    for (auto &item : hg38_data) {sam_unique << item << endl;}
    sam_unique.close();
    
    string v_log = v_output + "/log/unique_sort.log";
    string v_bam_out = v_output + v_sample_sub + ".unique.bam";
    v_log = v_output + "/log/sort.log";
    v_command = v_samtools + " sort -@ " + v_threads + " -m 1g " + merged_sam_unique  + " > " +  v_bam_out + " 2>" + v_log;
    system (v_command.c_str());
    v_log = v_output + "/log/unique_index.log";
    v_command = v_samtools + " index " + v_bam_out + " 2>" + v_log;
    system (v_command.c_str());
    
    
    string v_bam_nodup = v_output + v_sample_sub + ".unique.nodup.bam";;
    string v_bam_nodup_metrics = v_output + v_sample_sub + ".unique.nodup.metrics";
    
    v_log = v_output + "/log/unique_nodup.log";

	if (v_nanopore == 0) {
		std::size_t found = v_picard.find(".jar");
		if (found!=std::string::npos)
		{
			v_command = "java -jar " + v_picard + " MarkDuplicates I=" +  v_bam_out + " O=" + v_bam_nodup + " M=" + v_bam_nodup_metrics + " VALIDATION_STRINGENCY=SILENT 2>" + v_log;
		}
		else {
			v_command = v_picard + " MarkDuplicates I=" +  v_bam_out + " O=" + v_bam_nodup + " M=" + v_bam_nodup_metrics + " VALIDATION_STRINGENCY=SILENT 2>" + v_log;
		}
		
		v_log = v_output + "/log/unique_nodup_index.log";

		if (v_rna == 0) {
			system (v_command.c_str());
			v_command = v_samtools + " index " + v_bam_nodup + " 2>" + v_log;
			system (v_command.c_str());
		}
	}
    v_message = "Merging files: done";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    
    
    v_message = "Cleaning temporary files ... ";
    screen_message (screen_size, 0, v_message, 2, v_quiet);
    
    removefile(v_output + v_sample + "hg38.selected.sam",v_debug);
    removefile(v_output + v_sample + "hg38.bam",v_debug);
    removefile(v_output + v_sample + "hg38.unmapped.sam",v_debug);
    removefile(v_output + v_sample + "hg38.unmapped.bam",v_debug);

    removefile(v_output + v_sample_sub + ".adjusted.sam",v_debug);
    removefile(v_output + v_sample_sub + ".unique.sam",v_debug);
    removefile(v_output + v_sample + "R0.fastq",v_debug);
    removefile(v_output + v_sample + "R1.fastq",v_debug);
    removefile(v_output + v_sample + "R2.fastq",v_debug);
    removefile(v_output + v_sample + "selected.trim.R0.fastq",v_debug);
    removefile(v_output + v_sample + "selected.trim.R1.fastq",v_debug);
    removefile(v_output + v_sample + "selected.trim.R2.fastq",v_debug);
    removefile(v_output + v_sample + "mapped_r0.tmp.fq",v_debug);
    removefile(v_output + v_sample + "mapped_r1.tmp.fq",v_debug);
    removefile(v_output + v_sample + "mapped_r2.tmp.fq",v_debug);
    removefile(v_output + v_sample + "unmapped_r0.tmp.fq",v_debug);
    removefile(v_output + v_sample + "unmapped_r1.tmp.fq",v_debug);
    removefile(v_output + v_sample + "unmapped_r2.tmp.fq",v_debug);
    removefile(v_output + v_sample_sub + "_hg38.sam",v_debug);
    


    
    if (v_map_type == "paired") {
        v_command = "gzip -f " + v_output + v_sample_sub + "_selected.R1.fastq";
        v_system_out = GetStdoutFromCommand(v_command);
        v_command = "gzip -f " + v_output + v_sample_sub + "_selected.R2.fastq";
        v_system_out = GetStdoutFromCommand(v_command);
    }
    if (v_map_type == "single") {
        v_command = "gzip -f " + v_output + v_sample_sub + "_selected.R0.fastq";
        v_system_out = GetStdoutFromCommand(v_command);
    }
    
    v_command = "gzip -f " + v_output + v_sample_sub + "_addressing_table.txt";
    v_system_out = GetStdoutFromCommand(v_command);
    
    
    for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
    {
        
        if (*i == "") {continue;}
        removefile(v_output + v_sample + "sorted_" + *i + "_R1.fastq",v_debug);
        removefile(v_output + v_sample + "sorted_" + *i + "_R2.fastq",v_debug);
        removefile(v_output + v_sample + "sorted_" + *i + "_R0.fastq",v_debug);
        removefile(v_output + v_sample + *i + ".unique.sam",v_debug);
        removefile(v_output + v_sample + *i + ".adjusted.sam",v_debug);
        removefile(v_output + v_sample + *i + ".tmp.sam",v_debug);
        removefile(v_output + v_sample + *i + ".trim.R1.fastq",v_debug);
        removefile(v_output + v_sample + *i + ".trim.R2.fastq",v_debug);
        removefile(v_output + v_sample + *i + ".trim.R0.fastq",v_debug);
		
		if (v_nanopore == 1) {
			removefile(v_output + v_sample + *i + ".tmp.fq",v_debug);
		}
		

        if (v_rna == 1)
        {
            removefile(v_output + v_sample + *i + ".Log.out",v_debug);
            removefile(v_output + v_sample + *i + ".Log.final.out",v_debug);
            removefile(v_output + v_sample + *i + ".Log.progress.out",v_debug);
            removefile(v_output + v_sample + *i + ".SJ.out.tab",v_debug);
            removefile(v_output + v_sample + *i + ".Aligned.out.sam",v_debug);
        }


        v_command = "gzip -f " + v_output + v_sample + *i + "_R1.fastq";
        v_system_out = GetStdoutFromCommand(v_command);
        v_command = "gzip -f " + v_output + v_sample + *i + "_R2.fastq";
        v_system_out = GetStdoutFromCommand(v_command);
        v_command = "gzip -f " + v_output + v_sample + *i + "_R0.fastq";
        v_system_out = GetStdoutFromCommand(v_command);

    }
    
    
    
    general_log.close();
    int quiet_mem = v_quiet;
    
}
