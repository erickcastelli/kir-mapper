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
#include <mutex>



#include "map_dna.hpp"
#include "external.hpp"
#include "functions.hpp"
#include "ThreadPool.hpp"

mutex mtxp;

int minhit = 3;

void main_preselect ()
{

    int v_check = 0;
    if (v_r0 != "") {v_check = 1;}
    if ((v_r1 != "") && (v_r2 != "")) {v_check = 1;}
    if (v_bam != "") {v_check = 1;}
    if ((v_r1 != "") && (v_r1 == v_r2)) {warnings.push_back("r1 and r2 must be different fastq files"); v_check = 0;}
    if (v_sample == "") {warnings.push_back("You must indicate a sample name or id (-sample sample_id)"); v_check = 0;}

    if (((v_db == "") || (v_sample == "")) || (v_check == 0))
    {
        screen_message (screen_size, 0, "", 1, v_quiet);
        v_message = "Program:   " + Program_name + "::select";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        v_message = "Version:   " + Program_version + ", " + Program_date;
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        v_message = "Usage:     kir-mapper select -r1 R1.gz -r2 R2.gz -sample name <options>";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Mandatory options:", 1, v_quiet);
        screen_message (screen_size, 2, "-sample      sample name or identification", 1, v_quiet);
        screen_message (screen_size, 2, "-db          path to the kir-mapper database", 1, v_quiet);
        screen_message (screen_size, 2, "-r0          a single-ended fastq (fq, fastq, or gz)", 1, v_quiet);
        screen_message (screen_size, 2, "-r1          a paired-ended forward fastq (fq, fastq, or gz)", 1, v_quiet);
        screen_message (screen_size, 2, "-r2          a paired-ended reverse fastq (fq, fastq, or gz)", 1, v_quiet);
        screen_message (screen_size, 2, "-bam         a BAM file (ignore r0/r1/r2)", 1, v_quiet);
        screen_message (screen_size, 2, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Other options:", 1, v_quiet);
        screen_message (screen_size, 2, "-output      output folder", 1, v_quiet);
        screen_message (screen_size, 2, "-threads     number of threads [" + v_threads + "]", 1, v_quiet);
//        screen_message (screen_size, 2, "bed         BED file to override the one in the database", 1, v_quiet);
        screen_message (screen_size, 2, "-buffer      number of sequences in buffer [" + to_string(v_buffer) + "]", 1, v_quiet);
        screen_message (screen_size, 2, "-config      path to a kir-mapper configuration file", 1, v_quiet);
                       
        if (v_bwa != "") {v_message = "[found at " + v_bwa + "]";}
        if (v_bwa == "") {v_message = "(!!! bwa not detected !!!)";}
        v_message = "-bwa         path to BWA " + v_message;
        screen_message (screen_size, 2, v_message, 1, v_quiet);

          if (v_samtools != "") {v_message = "[found at " + v_samtools + "]";}
          if (v_samtools == "") {v_message = "[!!!Samtools not detected!!!]";}
          v_message = "-samtools    path to Samtools " + v_message;
          screen_message (screen_size, 2, v_message, 1, v_quiet);
        
        screen_message (screen_size, 2, "", 1, v_quiet);
        screen_message (screen_size, 2, "--rna                 this is RNA-Seq data", 1, v_quiet);
        screen_message (screen_size, 2, "--low-mem             force low memory mode", 1, v_quiet);
        screen_message (screen_size, 2, "--skip-unmapped       skip retrieve unmapped reads from BAM/SAM files", 1, v_quiet);
        screen_message (screen_size, 2, "--quiet               quiet mode", 1, v_quiet);
        screen_message (screen_size, 2, "--noconfig            ignore the pre-configuration file", 1, v_quiet);

        screen_message (screen_size, 0, "", 1, v_quiet);
        return;
    }
    
    
    if (v_callint == 0) {
        
        if (v_bam != "") {
            
			string file = v_bam;
			boost::to_upper(file);
			if (! ends_with(file,".BAM")) 
			{
				warnings.push_back ("This is not a BAM file."); 
				v_sample = "";
                v_r0 = "";
				v_r1 = "";
                v_r2 = "";
                v_bam = "";
                main_preselect();
                return;
			}
			
			ifstream file_bam(v_bam.c_str());
            if (!file_bam)
            {
                warnings.push_back ("Could not access the .bam file or file does not exist.");
                v_sample = "";
                v_r0 = "";
				v_r1 = "";
                v_r2 = "";
                v_bam = "";
                main_preselect();
                return;
            }
            file_bam.close();
        }


        if (v_r1 != "") {
            
			string file = v_r1;
			boost::to_upper(file);
			if (((ends_with(file,".FQ")) || (ends_with(file,".FQ.GZ"))) || ((ends_with(file,".FASTQ")) || (ends_with(file,".FASTQ.GZ"))))
			{}
			else {
				warnings.push_back ("This is not a FASTQ file."); 
				v_sample = "";
                v_r0 = "";
				v_r1 = "";
                v_r2 = "";
                v_bam = "";
                main_preselect();

                return;
			}
			
			ifstream file_r1(v_r1.c_str());
            if (!file_r1)
            {
                warnings.push_back ("Could not access the r1 file or file does not exist.");
                v_sample = "";
                v_r1 = "";
                main_preselect();
                return;
            }
            file_r1.close();
        }

        if (v_r0 != "") {
            
			string file = v_r0;
			boost::to_upper(file);
			if (((ends_with(file,".FQ")) || (ends_with(file,".FQ.GZ"))) || ((ends_with(file,".FASTQ")) || (ends_with(file,".FASTQ.GZ"))))
			{}
			else {
				warnings.push_back ("This is not a FASTQ file."); 
				v_sample = "";
                v_r1 = "";
                v_r2 = "";
                v_bam = "";
				v_r0 = "";
                main_preselect();
                return;
			}
			if (! fileExists(v_r0)) {
                warnings.push_back ("Could not access the r0 file or file does not exist.");
                v_sample = "";
                v_r1 = "";
				v_bam = "";
                main_preselect();
                return;
            }
        }
        
        
        if (v_r2 != "") {
            
			string file = v_r2;
			boost::to_upper(file);
			if (((ends_with(file,".FQ")) || (ends_with(file,".FQ.GZ"))) || ((ends_with(file,".FASTQ")) || (ends_with(file,".FASTQ.GZ"))))
			{}
			else {
				warnings.push_back ("This is not a FASTQ file."); 
				v_sample = "";
                v_r1 = "";
                v_r2 = "";
                v_bam = "";
				v_r0 = "";
                main_preselect();
                return;
			}
			
			ifstream file_r2(v_r2.c_str());
            if (!file_r2)
            {
                warnings.push_back ("Could not access the r2 file or file does not exist.");
                v_sample = "";
                v_r2 = "";
				v_r0 = "";
				v_bam = "";
                main_preselect();
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
		
		if (((ends_with(v_r1,".FQ")) && (ends_with(v_r2,".GZ"))) || ((ends_with(v_r2,".FQ")) && (ends_with(v_r1,".GZ"))))
        {
            
            warnings.push_back ("kir-mapper cannot deal with .fastq and .gz files simultaneously.");
            v_sample = "";
            v_r1 = "";
            main_dna_map();
            return;
        }
        
        if (((ends_with(v_r1,".FASTQ")) && (ends_with(v_r2,".GZ"))) || ((ends_with(v_r2,".FASTQ")) && (ends_with(v_r1,".GZ"))))
        {
            
            warnings.push_back ("kir-mapper cannot deal with .fastq and .gz files simultaneously.");
            v_sample = "";
            v_r1 = "";
            main_dna_map();
            return;
        }
					
        
        if (((v_r1 != "") && (v_r2 == "")) || ((v_r1 == "") && (v_r2 != "")))
        {
            warnings.push_back ("You must indicate both r1 and r2 when using paired-end sequencing data.");
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
            v_message = "Could not access database " + v_db;
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
                    warnings.push_back ("Could not access the BED file or file does not exist.");
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
        
        
        
        // criando output folder
        if (v_output == "")
        {
            if (v_bam != "") {v_output = findfilepath(v_bam) + "kir-mapper/";}
            if (v_r0 != "") {v_output = findfilepath(v_r0) + "kir-mapper/";}
            if (v_r1 != "") {v_output = findfilepath(v_r1) + "kir-mapper/";}
        }
        
        v_command = "mkdir " + v_output;
        v_system_out = GetStdoutFromCommand(v_command);
    //    v_command = "rm -rf " + v_output + "/*";
    //    v_system_out = GetStdoutFromCommand(v_command);
    //    v_command = "rm -rf " + v_output + "/log/";
     //   v_system_out = GetStdoutFromCommand(v_command);
        v_command = "mkdir " + v_output + "/log/";
        v_system_out = GetStdoutFromCommand(v_command);
    
    
        string v_sample_sub = v_sample.substr(0, v_sample.size()-1);
        
        
        string v_general_log = v_output + v_sample_sub + ".log";
        ofstream general_log;
        general_log.open (v_general_log.c_str());
        

        screen_message (screen_size, 0, "", 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        v_message = Program_name + "::select";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        general_log << v_message << endl;
     
        v_message = "Version " + Program_version;
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        general_log << v_message << endl;
       

        v_message = "Number of threads: " + v_threads;
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        general_log << v_message << endl;
            
        
        v_message = "Sample: " + v_sample_sub;
        general_log << v_message << endl;
        v_message = "Database: " + v_db;
        general_log << v_message << endl;
        v_message = "Output: " + v_output;
        general_log << v_message << endl;

   }
    
    
    screen_message (screen_size, 0, "Initializing read selection ... please wait ...", 2, v_quiet);

    unordered_map <string,int> motifs_master_pre;
    string motifdbfile = "";
    if (v_rna == 1) {motifdbfile = v_db + "/motif/rna/all/motif_all.txt";}
    if (v_rna == 0) {motifdbfile = v_db + "/motif/dna/all/motif_all.txt";}


    boost::replace_all(motifdbfile, "\\ ", " ");

    ifstream motifdb (motifdbfile.c_str());
    for( std::string line; getline( motifdb, line ); )
    {
        if (line == "") {continue;}
        motifs_master_pre[line] = 1;
    }
    motifdb.close();
    
    if ((v_r1 != "") && (v_r2 != "")) {v_map_type = "paired";}
    if (v_r0 != "") {v_map_type = "single";}
    if (v_bam != "") {v_map_type = "unknown";}
 
    if (v_bam != "")
    {
        
        screen_message (screen_size, 0, "Extracting mapped reads from BAM file ... please wait ...", 2, v_quiet);
        
        string v_command = v_samtools + " view -H " + v_bam;
        string head = GetStdoutFromCommand(v_command.c_str());
        vector <string> heads;
        boost::split(heads,head,boost::is_any_of("\n"));
        for (auto & item : heads)
        {
            vector <string> data;
            boost::split(data,item,boost::is_any_of("\t"));
            if (data[0] == "@SQ") {if (data[1].substr(3) == "6") {usechr = 0;}}
            if (data[0] == "@SQ") {if (data[1].substr(3) == "chr6") {usechr = 1;}}
        }
        
        v_bed = "'" + v_db + "/bed/select_dna.bed'";
        boost::replace_all(v_bed, "\\ ", " ");

 
        screen_message (screen_size, 0, "Extracting mapped reads from BAM file ... please wait ...", 2, v_quiet);
        
        string v_hg38_sam = v_output + v_sample + "hg38.sam";
        string v_hg38_bam = v_output + v_sample + "hg38.bam";
        string v_hg38_unmapped = v_output + v_sample + "hg38.unmapped.sam";
        string v_hg38_unmapped_bam = v_output + v_sample + "hg38.unmapped.bam";
        string log = v_output + "/log/hg38_extract.log";
        string v_r0_mapped_tmp = v_output + v_sample + "mapped_r0.tmp.fq";;
        string v_r1_mapped_tmp = v_output + v_sample + "mapped_r1.tmp.fq";;
        string v_r2_mapped_tmp = v_output + v_sample + "mapped_r2.tmp.fq";;
        string v_r0_unmapped_tmp = v_output + v_sample + "unmapped_r0.tmp.fq";;
        string v_r1_unmapped_tmp = v_output + v_sample + "unmapped_r1.tmp.fq";;
        string v_r2_unmapped_tmp = v_output + v_sample + "unmapped_r2.tmp.fq";;
        
        v_command = v_samtools + " view -@ " + v_threads + " -h -ML " + v_bed + " " + v_bam + " > " + v_hg38_sam + " 2>" + log;
        system (v_command.c_str());
        if (v_verbose == 1) {screen_message (screen_size, 0, "Converting mapped reads to FASTQ.. please wait ...", 1, v_quiet);}
        
        v_command = v_samtools + " fastq -0 " + v_r0_mapped_tmp + " -1 " + v_r1_mapped_tmp + " -2 " + v_r2_mapped_tmp  + " " + v_hg38_sam + " 2>" + log;
        system (v_command.c_str());


        v_map_type = "paired";
        double sizemR0 = filesize(v_r0_mapped_tmp.c_str());
        double sizemR1 = filesize(v_r1_mapped_tmp.c_str());
        double sizemR2 = filesize(v_r2_mapped_tmp.c_str());
        if (sizemR0 > sizemR1) {v_map_type = "single";}
        if (sizemR0 <= sizemR1) {v_map_type = "paired";}
        double sizeuR0 = 0;
        double sizeuR1 = 0;
        double sizeuR2 = 0;
        
 
        if ((sizemR0 == 0) && (sizemR1 == 0))
        {
            warnings.push_back ("There is not enough data to be processed.");
            warnings.push_back ("No mapped reads retrieved from BAM.");
            v_sample = "";
            v_r1 = "";
            v_bam = "";
            v_r2 = "";
            v_callint = 2;
            return;
        }

        string file = "";
        if (v_map_type == "paired") {file = v_r1_mapped_tmp;}
        if (v_map_type == "single") {file = v_r0_mapped_tmp;}

        screen_message (screen_size, 0, "Selecting mapped reads ... please wait ...", 2, v_quiet);
        ifstream r1 (file.c_str());
        double bit_processed = 0;
        double counter = 0;
        double filelength = filesize(file.c_str());
        
        ThreadPool pool(stoi(v_threads));
        std::vector< std::future<int> > results;
        
        map <string,string> selected_results;

        for( std::string line; getline( r1, line ); )
        {
            string idr1 = line;
            getline( r1, line );
            string seqr1 = line;
            getline( r1, line );
            string infor1 = line;
            getline( r1, line );
            string qualr1 = line;
            
            std::string s = "/1";
            std::string::size_type i = idr1.find(s);
            if (i != std::string::npos)
                idr1.erase(i, s.length());
            
            
            results.emplace_back(
            pool.enqueue([idr1, seqr1, infor1, qualr1, &motifs_master_pre, &selected_results,  &bit_processed, &counter, filelength]
            {
                string idr1i = idr1;
                string seqr1i = seqr1;
                string infor1i = infor1;
                string qualr1i = qualr1;


                if (seqr1i.size() < motif_size) {return 1;}
                if (seqr1i.size() < v_size) {return 1;}
 
                int hitcount = 0;
                for (int c = 0; c < (seqr1i.size() - motif_size); c++)
                {
                    string sub = seqr1i.substr(c,motif_size);

                    if ( motifs_master_pre.find(sub) != motifs_master_pre.end() ) {
                        hitcount++;
                        c = c + 1;
                    }
                    
                    if (hitcount == minhit) {
                        vector <string> idsub;
                        boost::split(idsub,idr1i,boost::is_any_of(" "));
                        string idgen = idsub[0].substr(1);
                        mtxp.lock();
                        selected_results[idgen] = idr1 + "\n" + seqr1 + "\n" + infor1 + "\n" + qualr1;
                        mtxp.unlock();
                        break;
                    }
                }
                mtxp.lock();
                counter++;
                bit_processed = bit_processed + idr1i.length() + seqr1i.length() + qualr1i.length() + infor1i.length() + 4;
                if (counter == v_buffer)
                {
                   v_message = "Selecting mapped reads ... " + to_string((((bit_processed) / filelength) * 100)) + " %";
                   screen_message (screen_size, 0, v_message, 2, v_quiet);
                   counter = 0;
                }
                mtxp.unlock();
                return 1;
                })
                );
        }
        r1.close();
        for(auto && result: results){result.get();} // waiting for all threads

        
        
        
        
        if (v_skip_unmapped == 0) {

            screen_message (screen_size, 0, "Extracting unmapped reads from BAM file ... please wait ...", 2, v_quiet);
            
            v_command = v_samtools + " view -@ " + v_threads + " -f4 -h " + v_bam + " > " + v_hg38_unmapped +  " 2>" + log;
            system (v_command.c_str());
            
            
            v_command = v_samtools + " fastq -0 " + v_r0_unmapped_tmp + " -1 " + v_r1_unmapped_tmp + " -2 " +  v_r2_unmapped_tmp  + " " + v_hg38_unmapped + " 2>" + log;
            system (v_command.c_str());

            
            sizeuR0 = filesize(v_r0_unmapped_tmp.c_str());
            sizeuR1 = filesize(v_r1_unmapped_tmp.c_str());
            sizeuR2 = filesize(v_r2_unmapped_tmp.c_str());
            if (v_map_type == "paired") {file = v_r1_unmapped_tmp;}
            if (v_map_type == "single") {file = v_r0_unmapped_tmp;}

            filelength = filesize(file.c_str());
            bit_processed = 0;
            counter = 0;
 
            screen_message (screen_size, 0, "Selecting unmapped reads ... please wait ...", 2, v_quiet);

            ThreadPool pool(stoi(v_threads));
            std::vector< std::future<int> > results;
            
            ifstream r1 (file.c_str());
            for( std::string line; getline( r1, line ); )
            {
                string idr1 = line;
                getline( r1, line );
                string seqr1 = line;
                getline( r1, line );
                string infor1 = line;
                getline( r1, line );
                string qualr1 = line;
                
                std::string s = "/1";
                std::string::size_type i = idr1.find(s);
                if (i != std::string::npos)
                    idr1.erase(i, s.length());
                
                results.emplace_back(
                pool.enqueue([idr1, seqr1, infor1, qualr1, &motifs_master_pre, &selected_results, &bit_processed, &counter, filelength]
                {
                    string idr1i = idr1;
                    string seqr1i = seqr1;
                    string infor1i = infor1;
                    string qualr1i = qualr1;

                    if (seqr1i.size() < motif_size) {return 1;}
                    if (seqr1i.size() < v_size) {return 1;}
     
                    int hitcount = 0;
                    for (int c = 0; c < (seqr1i.size() - motif_size); c++)
                    {
                        string sub = seqr1i.substr(c,motif_size);

                        if ( motifs_master_pre.find(sub) != motifs_master_pre.end() ) {
                            hitcount++;
                            c = c + 1;
                        }
                        
                        if (hitcount == minhit) {
                            vector <string> idsub;
                            boost::split(idsub,idr1i,boost::is_any_of(" "));
                            string idgen = idsub[0].substr(1);
                            mtxp.lock();
                            selected_results[idgen] = idr1 + "\n" + seqr1 + "\n" + infor1 + "\n" + qualr1;
                            mtxp.unlock();
                            break;
                        }
                    }
                    
                    mtxp.lock();
                    counter++;
                    bit_processed = bit_processed + idr1i.length() + seqr1i.length() + qualr1i.length() + infor1i.length() + 4;
                    if (counter == v_buffer)
                    {
                       v_message = "Selecting unmapped reads ... " + to_string((((bit_processed) / filelength) * 100)) + " %";
                       screen_message (screen_size, 0, v_message, 2, v_quiet);
                       counter = 0;
                    }
                    mtxp.unlock();
                    return 1;
                    })
                    );
            }
            for(auto && result: results){result.get();} // waiting for all threads
            r1.close();
 
        }

        
        
        if (v_map_type == "paired")
        {
            screen_message (screen_size, 0, "Reading mapped file R2 ... please wait ...", 2, v_quiet);

            map <string,string> results_selected_pair;
            
            ThreadPool pool(stoi(v_threads));
            std::vector< std::future<int> > results;
  
            
            file = v_r2_mapped_tmp;
            ifstream r1 (file.c_str());
            for( std::string line; getline( r1, line ); )
            {
                string idr1 = line;
                getline( r1, line );
                string seqr1 = line;
                getline( r1, line );
                string infor1 = line;
                getline( r1, line );
                string qualr1 = line;
                
                std::string s = "/2";
                std::string::size_type i = idr1.find(s);
                if (i != std::string::npos)
                    idr1.erase(i, s.length());

                
                results.emplace_back(
                pool.enqueue([idr1, seqr1, infor1, qualr1, &results_selected_pair, &selected_results]
                {
                    string idr1i = idr1;
                    string seqr1i = seqr1;
                    string infor1i = infor1;
                    string qualr1i = qualr1;
                    
                    vector <string> idsub;
                    boost::split(idsub,idr1i,boost::is_any_of(" "));
                    string idgen = idsub[0].substr(1);

                    if (selected_results.find(idgen) != selected_results.end())
                    {
                        mtxp.lock();
                        results_selected_pair[idgen] = idr1i + "\n" + seqr1i + "\n" + infor1i + "\n" + qualr1i;
                        mtxp.unlock();
                    }
                    return 0;
                })
                );
            }
            for(auto && result: results){result.get();} // waiting for all threads
            r1.close();
            
            
            if (v_skip_unmapped == 0)
            {
                screen_message (screen_size, 0, "Reading unmapped file R2 ... please wait ...", 2, v_quiet);
                file = v_r2_unmapped_tmp;
                ifstream r1 (file.c_str());
                ThreadPool pool(stoi(v_threads));
                std::vector< std::future<int> > results;
     
                for( std::string line; getline( r1, line ); )
                {
                    string idr1 = line;
                    getline( r1, line );
                    string seqr1 = line;
                    getline( r1, line );
                    string infor1 = line;
                    getline( r1, line );
                    string qualr1 = line;
                    
                    std::string s = "/2";
                    std::string::size_type i = idr1.find(s);
                    if (i != std::string::npos)
                        idr1.erase(i, s.length());
                    
                    results.emplace_back(
                    pool.enqueue([idr1, seqr1, infor1, qualr1, &results_selected_pair, &selected_results]
                    {
                        string idr1i = idr1;
                        string seqr1i = seqr1;
                        string infor1i = infor1;
                        string qualr1i = qualr1;
                    
                        vector <string> idsub;
                        boost::split(idsub,idr1i,boost::is_any_of(" "));
                        string idgen = idsub[0].substr(1);

                        if (selected_results.find(idgen) != selected_results.end())
                        {
                            mtxp.lock();
                            results_selected_pair[idgen] = idr1i + "\n" + seqr1i + "\n" + infor1i + "\n" + qualr1i;
                            mtxp.unlock();
                        }
                        return 0;
                    })
                    );
                }
                for(auto && result: results){result.get();} // waiting for all threads
                r1.close();
            }
            
            screen_message (screen_size, 0, "Writing selected reads ... please wait ...", 2, v_quiet);

            ofstream OUTpairedR1;
            ofstream OUTpairedR2;
            string pairedR1 = v_output + v_sample + "selected.R1.fastq";;
            string pairedR2 = v_output + v_sample + "selected.R2.fastq";
            
  
            OUTpairedR1.open (pairedR1.c_str());
            OUTpairedR2.open (pairedR2.c_str());
            
            for (auto & item : selected_results)
            {
                string id = item.first;
                if (results_selected_pair.find(id) != results_selected_pair.end())
                {
                    OUTpairedR1 << item.second <<  endl;
                    OUTpairedR2 << results_selected_pair[id] << endl;
                }
            }
            
            OUTpairedR1.close();
            OUTpairedR2.close();
  
        }

        
        
        if (v_map_type == "single")
        {
 
            screen_message (screen_size, 0, "Writing selected reads ... please wait ...", 2, v_quiet);

            ofstream OUTpairedR1;
            string pairedR1 = v_output + v_sample + "selected.R0.fastq";
 
            OUTpairedR1.open (pairedR1.c_str());
            
            for (auto & item : selected_results)
            {
                string id = item.first;
                OUTpairedR1 << item.second <<  endl;
            }
            
            OUTpairedR1.close();
        }
  
    }
    
    

    
    
    
    
    
    if (v_bam == "") {
    
        screen_message (screen_size, 0, "Selecting reads ...  ... please wait ...", 2, v_quiet);
        
        if (v_map_type == "paired"){

            vector <string> selected_results_r1_pre;
            vector <string> selected_results_r2_pre;
        
            if ( (ends_with(v_r1,".fq")) || (ends_with(v_r1,".fastq")) )
            {
                
                ThreadPool pool(stoi(v_threads));
                std::vector< std::future<int> > results;
                
                double filelength = filesize(v_r1.c_str());
  
                ifstream r1 (v_r1.c_str());
                ifstream r2 (v_r2.c_str());
                double bit_processed = 0;
                
                double counter = 0;
 
                for( std::string line; getline( r1, line ); )
                {
                    
                    if (counter == v_buffer) {for(auto && result: results){result.get();} results.clear();} // waiting for all threads}
                    
                    string idr1 = line;
                    getline( r1, line );
                    string seqr1 = line;
                    getline( r1, line );
                    string infor1 = line;
                    getline( r1, line );
                    string qualr1 = line;
                    
                    boost::replace_all(idr1, "/1", "");
                    boost::replace_all(idr1, "/2", "");

 
                    
                    getline( r2, line );
                    string idr2 = line;
                    getline( r2, line );
                    string seqr2 = line;
                    getline( r2, line );
                    string infor2 = line;
                    getline( r2, line );
                    string qualr2 = line;
                    
                    boost::replace_all(idr2, "/1", "");
                    boost::replace_all(idr2, "/2", "");

                    
                    results.emplace_back(
                    pool.enqueue([idr1, seqr1, infor1, qualr1, idr2, seqr2, infor2, qualr2, &motifs_master_pre, &selected_results_r1_pre, &selected_results_r2_pre, &bit_processed,&counter,filelength]
                    {
                        string idr1i = idr1;
                        string seqr1i = seqr1;
                        string infor1i = infor1;
                        string qualr1i = qualr1;
                        string idr2i = idr2;
                        string seqr2i = seqr2;
                        string infor2i = infor2;
                        string qualr2i = qualr2;

  
                        if (seqr1i.size() < motif_size) {return 1;}
                        if (seqr2i.size() < motif_size) {return 1;}
                        if (seqr1i.size() < v_size) {return 1;}
                        if (seqr2i.size() < v_size) {return 1;}

                        int hitcount = 0;
                        for (int c = 0; c < (seqr1i.size() - motif_size); c++)
                        {
                            string sub = seqr1i.substr(c,motif_size);
        
                            if ( motifs_master_pre.find(sub) != motifs_master_pre.end() ) {hitcount++;c++;}
                            if (hitcount == minhit)
                            {
                                vector <string> idsub;
                                boost::split(idsub,idr1i,boost::is_any_of(" "));
                                string idgen = idsub[0].substr(1);
                                
                                mtxp.lock();
                                selected_results_r1_pre.push_back(idr1i + "\n" + seqr1i + "\n+\n" + qualr1i);
                                selected_results_r2_pre.push_back(idr2i + "\n" + seqr2i + "\n+\n" + qualr2i);
                                mtxp.unlock();
                                break;
                            }
                            c++;
                        }
                        mtxp.lock();
                        counter++;
                        bit_processed = bit_processed + idr1i.length() + seqr1i.length() + qualr1i.length() + infor1i.length() + 4;
                        if (counter == v_buffer)
                        {
                           v_message = "Selecting reads ... " + to_string((((bit_processed) / filelength) * 100)) + " %";
                           screen_message (screen_size, 0, v_message, 2, v_quiet);
                           counter = 0;
                        }
                        mtxp.unlock();
                        return 1;
                        })
                        );
                }
                r1.close();
                r2.close();
                for(auto && result: results){result.get();} // waiting for all threads
            }

            
            if ( ends_with(v_r1,".gz"))
                {
                    std::ifstream r1(v_r1, std::ios_base::in | std::ios_base::binary);
                    std::ifstream r2(v_r2, std::ios_base::in | std::ios_base::binary);

                    ThreadPool pool(stoi(v_threads));
                    std::vector< std::future<int> > results;
     
                    try {
                    boost::iostreams::filtering_istream in1;
                    boost::iostreams::filtering_istream in2;
                    in1.push(boost::iostreams::gzip_decompressor());
                    in1.push(r1);
                    in2.push(boost::iostreams::gzip_decompressor());
                    in2.push(r2);
                        
                    int counter = 0;
                    int master_counter = 0;
                    for(std::string str; std::getline(in1, str); )
                    {
                        string idr1 = str;
                        getline( in1, str );
                        string seqr1 = str;
                        getline( in1, str );
                        string infor1 = str;
                        getline( in1, str );
                        string qualr1 = str;
                        
                        getline( in2, str );
                        string idr2 = str;
                        getline( in2, str );
                        string seqr2 = str;
                        getline( in2, str );
                        string infor2 = str;
                        getline( in2, str );
                        string qualr2 = str;
                        
                        
                        boost::replace_all(idr1, "/1", "");
                        boost::replace_all(idr1, "/2", "");
                        boost::replace_all(idr2, "/1", "");
                        boost::replace_all(idr2, "/2", "");
                        
                        results.emplace_back(
                        pool.enqueue([idr1, seqr1, infor1, qualr1, idr2, seqr2, infor2, qualr2, &motifs_master_pre, &selected_results_r1_pre, &selected_results_r2_pre, &counter, &master_counter]
                        {
                            string idr1i = idr1;
                            string seqr1i = seqr1;
                            string infor1i = infor1;
                            string qualr1i = qualr1;
                            string idr2i = idr2;
                            string seqr2i = seqr2;
                            string infor2i = infor2;
                            string qualr2i = qualr2;

                            if (seqr1i.size() < motif_size) {return 1;}
                            if (seqr2i.size() < motif_size) {return 1;}
                            if (seqr1i.size() < v_size) {return 1;}
                            if (seqr2i.size() < v_size) {return 1;}

                            int hitcount = 0;
                            for (int c = 0; c < (seqr1i.size() - motif_size); c++)
                            {
                                string sub = seqr1i.substr(c,motif_size);
            
                                if ( motifs_master_pre.find(sub) != motifs_master_pre.end() ) {hitcount++;c++;}
                                    
                                if (hitcount == minhit) {
                                    vector <string> idsub;
                                    boost::split(idsub,idr1i,boost::is_any_of(" "));
                                    string idgen = idsub[0].substr(1);

                                     mtxp.lock();
                                     selected_results_r1_pre.push_back(idr1i + "\n" + seqr1i + "\n+\n" + qualr1i);
                                     selected_results_r2_pre.push_back(idr2i + "\n" + seqr2i + "\n+\n" + qualr2i);
                                     mtxp.unlock();
                                     break;
                                }
                                c++;
                            }
                            
                            mtxp.lock();
                            counter++;
                            master_counter++;
                            if (counter == v_buffer)
                            {
                               v_message = "Selecting reads ... " + to_string(master_counter) + " have been processed ...";
                               screen_message (screen_size, 0, v_message, 2, v_quiet);
                               counter = 0;
                            }
                            mtxp.unlock();
                            
                            return 1;
                            })
                            );
                    }
                    r1.close();
                    r2.close();
                    for(auto && result: results){result.get();} // waiting for all threads
                }
                catch(const boost::iostreams::gzip_error& e) {
                    std::cout << e.what() << '\n';
                }
                    
            }
            
            
            
            
            ofstream OUTselectR1;
            ofstream OUTselectR2;
            string selectR1 = v_output + v_sample + "selected.R1.fastq";
            string selectR2 = v_output + v_sample + "selected.R2.fastq";
            OUTselectR1.open (selectR1.c_str());
            OUTselectR2.open (selectR2.c_str());
 
            for (int a = 0; a < selected_results_r1_pre.size(); a++)
            {
                OUTselectR1 << selected_results_r1_pre[a] << endl;
                OUTselectR2 << selected_results_r2_pre[a] << endl;
            }
            
            OUTselectR1.close();
            OUTselectR2.close();
            selected_results_r1_pre.clear();
            selected_results_r2_pre.clear();

        }
 
        
        
        
        
        if (v_map_type == "single"){
            
            vector <string> selected_results_r1_pre;
            vector <string> selected_results_r2_pre;
            
            ThreadPool pool(stoi(v_threads));
            std::vector< std::future<int> > results;
            
            if ((ends_with(v_r0,".fq")) || (ends_with(v_r0,".fastq")))
            {
                ifstream r1 (v_r0.c_str());
                double bit_processed = 0;
                double filelength = filesize(v_r0.c_str());
                double counter = 0;

                for( std::string line; getline( r1, line ); )
                {
                    string idr1 = line;
                    getline( r1, line );
                    string seqr1 = line;
                    getline( r1, line );
                    string infor1 = line;
                    getline( r1, line );
                    string qualr1 = line;
                    
                    results.emplace_back(
                    pool.enqueue([idr1, seqr1, infor1, qualr1,  &motifs_master_pre, &selected_results_r1_pre, &selected_results_r2_pre, &bit_processed,&counter,filelength]
                    {
                        string idr1i = idr1;
                        string seqr1i = seqr1;
                        string infor1i = infor1;
                        string qualr1i = qualr1;

                        if (seqr1i.size() < motif_size) {return 1;}
                        if (seqr1i.size() < v_size) {return 1;}

                        int hitcount = 0;
                        for (int c = 0; c < (seqr1i.size() - motif_size); c++) {
                        
                            string sub = seqr1i.substr(c,motif_size);
        
                            if ( motifs_master_pre.find(sub) != motifs_master_pre.end() ) {hitcount++;c++;}
                            if (hitcount == minhit) {
                                vector <string> idsub;
                                boost::split(idsub,idr1i,boost::is_any_of(" "));
                                string idgen = idsub[0].substr(1);
                                mtxp.lock();
                                selected_results_r1_pre.push_back(idr1i + "\n" + seqr1i + "\n+\n" + qualr1i);
                                mtxp.unlock();
                                break;
                            }
                            c++;
                        }
                        mtxp.lock();
                        counter++;
                        bit_processed = bit_processed + idr1i.length() + seqr1i.length() + qualr1i.length() + infor1i.length() + 4;

                        if (counter == v_buffer)
                        {
                           v_message = "Selecting reads ... " + to_string((((bit_processed) / filelength) * 100)) + " %";
                           screen_message (screen_size, 0, v_message, 2, v_quiet);
                           counter = 0;
                        }
                        mtxp.unlock();
                        return 1;
                        })
                        );
                }
                r1.close();
                for(auto && result: results){result.get();} // waiting for all threads
            }

            
            if ( ends_with(v_r0,".gz"))
                {
                    std::ifstream r1(v_r0, std::ios_base::in | std::ios_base::binary);

                    try {
                    boost::iostreams::filtering_istream in1;
                    in1.push(boost::iostreams::gzip_decompressor());
                    in1.push(r1);
                        
                    int counter = 0;
                    int master_counter = 0;
                    for(std::string str; std::getline(in1, str); )
                    {
                        string idr1 = str;
                        getline( in1, str );
                        string seqr1 = str;
                        getline( in1, str );
                        string infor1 = str;
                        getline( in1, str );
                        string qualr1 = str;
                        
                        results.emplace_back(
                        pool.enqueue([idr1, seqr1, infor1, qualr1, &motifs_master_pre, &selected_results_r1_pre, &selected_results_r2_pre, &counter,&master_counter]
                        {
                            string idr1i = idr1;
                            string seqr1i = seqr1;
                            string infor1i = infor1;
                            string qualr1i = qualr1;

                            if (seqr1i.size() < motif_size) {return 1;}
                            if (seqr1i.size() < v_size) {return 1;}
 
                            int hitcount = 0;
                            for (int c = 0; c < (seqr1i.size() - motif_size); c++)
                            {
                                string sub = seqr1i.substr(c,motif_size);
            
                                if ( motifs_master_pre.find(sub) != motifs_master_pre.end() ) {hitcount++;c++;}
                                if (hitcount == minhit) {
                                    vector <string> idsub;
                                    boost::split(idsub,idr1i,boost::is_any_of(" "));
                                    string idgen = idsub[0].substr(1);
                                    mtxp.lock();
                                    selected_results_r1_pre.push_back(idr1i + "\n" + seqr1i + "\n+\n" + qualr1i);
                                    mtxp.unlock();
                                break;
                                }
                                c++;
                            }
                            
                            mtxp.lock();
                            counter++;
                            master_counter++;
                            if (counter == v_buffer)
                            {
                               v_message = "Selecting reads ... " + to_string(master_counter) + " sequences done ...";
                               screen_message (screen_size, 0, v_message, 2, v_quiet);
                               counter = 0;
                            }
                            mtxp.unlock();
                            
                            return 1;
                            })
                            );
                    }
                    r1.close();
                    for(auto && result: results){result.get();} // waiting for all threads
                }
                catch(const boost::iostreams::gzip_error& e) {
                    std::cout << e.what() << '\n';
                }
            }
            ofstream OUTselectR1;
            ofstream OUTselectcleanR1;
            string selectR1 = v_output + v_sample + "selected.R0.fastq";
            OUTselectR1.open (selectR1.c_str());
            
            for (int a = 0; a < selected_results_r1_pre.size(); a++)
            {
                OUTselectR1 << selected_results_r1_pre[a] << endl;
            }
            OUTselectR1.close();
            selected_results_r1_pre.clear();
        }
         screen_message (screen_size, 0, "Selecting reads: done", 1, v_quiet);
    }
    
    
    
    

    
    
    screen_message (screen_size, 0, "Trimming selected reads ... please wait ...", 2, v_quiet);

    if (v_map_type == "paired")
    {
        string selectR1 = v_output + v_sample + "selected.R1.fastq";
        string selectR2 = v_output + v_sample + "selected.R2.fastq";
        map <string,string> clean_data_r1;
        map <string,string> clean_data_r2;
        
        thread t1 ([selectR1, &clean_data_r1]{
            ifstream input( selectR1.c_str() );
            for( std::string line; getline( input, line ); )
            {
                string idr1 = line;
                getline( input, line );
                string seqr1 = line;
                getline( input, line );
                string infor1 = line;
                getline( input, line );
                string qualr1 = line;
                
                
                string v_limits = mtrim(seqr1, qualr1);
                vector <string> idsub;
                boost::split(idsub,idr1,boost::is_any_of(" "));
                string idgen = idsub[0].substr(1);
                boost::replace_all(idgen, "/1", "");
                boost::replace_all(idgen, "/2", "");
                if (v_limits != "") {
                    mtxp.lock();
                    clean_data_r1[idgen] = idr1 + "\n" + v_limits;
                    mtxp.unlock();
                }
            }
            input.close();
         });

        
        thread t2 ([selectR2, &clean_data_r2]{
            ifstream input( selectR2.c_str() );
            for( std::string line; getline( input, line ); )
            {
                string idr1 = line;
                getline( input, line );
                string seqr1 = line;
                getline( input, line );
                string infor1 = line;
                getline( input, line );
                string qualr1 = line;
                
                
                string v_limits = mtrim(seqr1, qualr1);
                vector <string> idsub;
                boost::split(idsub,idr1,boost::is_any_of(" "));
                string idgen = idsub[0].substr(1);
                boost::replace_all(idgen, "/1", "");
                boost::replace_all(idgen, "/2", "");

                if (v_limits != "") {
                    mtxp.lock();
                    clean_data_r2[idgen] = idr1 + "\n" + v_limits;
                    mtxp.unlock();
                }
            }
            input.close();
         });
 
        t1.join();
        t2.join();
    
        string selectcleanR1 = v_output + v_sample + "selected.trim.R1.fastq";
        string selectcleanR2 = v_output + v_sample + "selected.trim.R2.fastq";

        ofstream OUTR1;
        ofstream OUTR2;
 
        OUTR1.open (selectcleanR1.c_str());
        OUTR2.open (selectcleanR2.c_str());

        for (auto & item : clean_data_r1)
        {
            if ((clean_data_r2[item.first] != "") && (clean_data_r1[item.first] != ""))
            {
                OUTR1 << item.second << endl;
                OUTR2 << clean_data_r2[item.first] << endl;
                vector <string> seqs;
                boost::split(seqs,item.second,boost::is_any_of("\n"));
                sequence_size_r1[item.first] = seqs[1].size();
                
                
                seqs.clear();
                boost::split(seqs,clean_data_r2[item.first],boost::is_any_of("\n"));
                sequence_size_r2[item.first] = seqs[1].size();
            }
        }
        OUTR1.close();
        OUTR2.close();
        
    }

    
    
    if (v_map_type == "single")
    {
        string selectR1 = v_output + v_sample + "selected.R0.fastq";
        map <string,string> clean_data_r1;
        
            ifstream input( selectR1.c_str() );
            for( std::string line; getline( input, line ); )
            {
                string idr1 = line;
                getline( input, line );
                string seqr1 = line;
                getline( input, line );
                string infor1 = line;
                getline( input, line );
                string qualr1 = line;
                
                
                string v_limits = mtrim(seqr1, qualr1);
                vector <string> idsub;
                boost::split(idsub,idr1,boost::is_any_of(" "));
                string idgen = idsub[0].substr(1);
                boost::replace_all(idgen, "/1", "");
                boost::replace_all(idgen, "/2", "");
                if (v_limits != "") {
                    mtxp.lock();
                    clean_data_r1[idgen] = idr1 + "\n" + v_limits;
                    mtxp.unlock();
                }
            }
            input.close();
 
        
        string selectcleanR1 = v_output + v_sample + "selected.trim.R0.fastq";
 
        ofstream OUTR1;
 
        OUTR1.open (selectcleanR1.c_str());
 
        for (auto & item : clean_data_r1)
        {
            if ((clean_data_r1[item.first] != ""))
            {
                OUTR1 << item.second << endl;
                
                vector <string> seqs;
                boost::split(seqs,item.second,boost::is_any_of("\n"));
                sequence_size_r1[item.first] = seqs[0].size();
            }
        }
        OUTR1.close();
    }
    
    
}
    
      
