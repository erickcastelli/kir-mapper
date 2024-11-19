//  kir-mapper
//
//  Created by Erick C. Castelli
//  2024 GeMBio.Unesp.
//  erick.castelli@unesp.br
// Contributions from code on the web


#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <mutex>
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <sys/stat.h>

#include "functions.hpp"
#include "external.hpp"
#include "ThreadPool.hpp"

namespace fs = filesystem;
using namespace std;


string v_regions;
string v_ref;
string v_resources;
mutex mtx_genotype;


void main_genotype() {
	
    
    int v_check = 0;
    if (v_output != "") {v_check = 1;}
 
    if (((v_db == "") || (v_check == 0)))
    {
        screen_message (screen_size, 0, "", 1, v_quiet);
        v_message = "Program:   " + Program_name + "::genotype";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        v_message = "Version:   " + Program_version + ", " + Program_date;
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        v_message = "Usage:     kir-mapper genotype -output output_folder <options>";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Mandatory options:", 1, v_quiet);
        screen_message (screen_size, 2, "-output         output folder (same as map and ncopy)", 1, v_quiet);
        screen_message (screen_size, 2, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Other options:", 1, v_quiet);
        screen_message (screen_size, 2, "-db              path to the kir-mapper database", 1, v_quiet);
        screen_message (screen_size, 2, "-threads         number of threads [" + v_threads + "]", 1, v_quiet);
        screen_message (screen_size, 2, "-target          which gene should be genotyped (e.g., -target KIR2DL1)", 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        screen_message (screen_size, 2, "--full          full genotype, include introns", 1, v_quiet);
        screen_message (screen_size, 2, "--quiet         quiet mode", 1, v_quiet);
        screen_message (screen_size, 2, "--nopolyphase   skip phasing tri/tetraploids", 1, v_quiet);
        screen_message (screen_size, 2, "--two_copies    force two copies when KIR gene is present", 1, v_quiet);
        screen_message (screen_size, 2, "--update_calls  update the allele calls keeping the current VCF", 1, v_quiet);
		
		

        screen_message (screen_size, 0, "", 1, v_quiet);
        return;
    }
    
    if (v_full == 1)
    {
        v_exome = 0;
    }
    if (v_full == 0)
    {
        v_exome = 1;
    }

    int precheck = 1;
    
    debug_message("Checking database - start");

    // checking database
    string v_db_info = v_db + "/db_dna.info";
    boost::replace_all(v_db_info, "\\ ", " ");
    if (! fileExists(v_db_info))
    {
        v_message = "Error accessing database " + v_db_info;
        cout << v_message << endl;
        warnings.push_back (v_message);
        precheck = 0;
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
        precheck = 0;
    }
    
   debug_message("Checking database - done");




    debug_message("Checking output - start");
    
    if (v_output == "") {
        v_message = "You must indicate a valid output folder, with 'map' and 'ncopy'";
        warnings.push_back (v_message);
        precheck = 0;
    }

    if (((! fileExists(v_output)) || (! fileExists(v_output + "/map"))) || (! fileExists(v_output + "/ncopy")))
    {
        v_message = "You must indicate a valid output folder, with 'map' and 'ncopy'";
        warnings.push_back (v_message);
        precheck = 0;
    }
    
    debug_message("Checking output - done");


    
    

    debug_message("Checking hg38 - start");

    v_ref = v_db + "/reference/hg38/reference.fasta";
    if (! fileExists(v_ref))
    {
        v_message = "There is something wrong with this database!";
        warnings.push_back (v_message);
        precheck = 0;
    }
      debug_message("Checking hg38 - done");


    if (precheck == 0){ v_db = ""; main_genotype();return; }
	


    debug_message("Creating output structure - start");

    // Creating output structure
	string cmd = "mkdir " + v_output;
	GetStdoutFromCommand(cmd.c_str());
	cmd = "mkdir " + v_output + "/genotype/";
	GetStdoutFromCommand(cmd.c_str());
    string v_output_adjusted = "";
	
    if (v_exome == 1)
    {
        cmd = "mkdir " + v_output + "/genotype/cds/";
        GetStdoutFromCommand(cmd.c_str());
        v_output_adjusted = v_output + "/genotype/cds/";
    }
    
    if (v_exome == 0)
    {
        cmd = "mkdir " + v_output + "/genotype/full/";
        GetStdoutFromCommand(cmd.c_str());
        v_output_adjusted = v_output + "/genotype/full/";
    }
    v_resources = v_db + "/genotype/";
    int include_sec = 0;

    debug_message("Creating output structure - done");


    screen_message (screen_size, 0, "", 1, v_quiet);
    v_message = "Program:   " + Program_name + "::genotype";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "Version:   " + Program_version + ", " + Program_date;
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "Cores:     using " + v_threads + " thread(s)";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    if (v_exome == 1) {
        v_message = "Mode:      only CDS";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
    }
    if (v_exome == 0) {
        v_message = "Mode:      full genotype";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
    }

    
    

    
    debug_message("Loading regions and targets - start");

	// Loading regions list
	
    v_message = " > Loading targets and regions";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    
    map <pair<string, string>,int> regions;
    map <string,int> valid_genes;
    vector <string> filedb;
    string path_map = v_db + "/genotype/bed/";
    string typegen = "";
    if (v_exome == 1) {typegen = ".CDS.";}
    if (v_exome == 0) {typegen = ".full.";}
 
    for (const auto & entry : fs::directory_iterator(path_map))
    {
        if (base_name(entry.path()).substr(0,1) == ".") {continue;}

        if (base_name(entry.path()).find(typegen) != std::string::npos)
        {
            string gene = base_name(entry.path());
            boost::replace_all(gene, ".CDS.bed", "");
            boost::replace_all(gene, ".full.bed", "");
            ifstream bedfile (entry.path().c_str());
            for( std::string line; getline(bedfile, line ); )
            {
                if (line == "") {continue;}
                vector <string> data;
                boost::split(data,line,boost::is_any_of("\t"));
                if (data.size() == 0) {continue;}
                pair <string,string> k;
                k = make_pair(gene,data[0] + ":" + data[1] + "-" + data[2]);
                valid_genes[gene] = 1;
                regions[k] = 1;
            }
            bedfile.close();
        }
    }
    debug_message("Loading regions and targets - end");
   
    
    
    
    debug_message("Loading copy numbers - start");

    
    
    
    // Loading Copy Numbers
    v_message = " > Loading copy numbers";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    map <pair<string,string>,string> copynumber;
    string copyfilesource = v_output + "/ncopy/copy_numbers.txt";
    if (! fileExists(copyfilesource))
    {
        v_message = "Error detecting the copy number file.";
        warnings.push_back (v_message);
        v_message = "Is this a kir-mapper ncopy output?";
        warnings.push_back (v_message);
        return;
    }
      
    
    map <string,int> failed_samples;
    map <string,int> highcp;
    ifstream copyfile (copyfilesource .c_str());
    for( std::string line; getline( copyfile, line ); )
    {
        if (line == "") {continue;}
        if (line.find("COPY_NUMBER") != std::string::npos) {continue;}
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));
        pair <string,string> k;
        k = make_pair(data[0],data[1]);
        boost::replace_all(data[2], "+", "");
        if (data[2] == "NA")
        {
            failed_samples[data[0]] = 1;
            continue;
        }
        
        string cpcheck = data[2];
        if ((v_force_two_copies == 1) && (cpcheck != "0")) {cpcheck = "2";}

        copynumber[k] = cpcheck;
        if (highcp[data[1]] < stoi(cpcheck)) {highcp[data[1]] = stoi(cpcheck);}
    }
    copyfile.close();


    debug_message("Loading copy numbers - done");

    
    debug_message("Checking targets - start");
    map <string,int> targets;
    if (v_target != "")
    {
        vector <string> sub;
        boost::replace_all(v_target, " ", "");
        boost::replace_all(v_target, ",,", ",");
        boost::split(sub,v_target,boost::is_any_of(","));
        for (auto item : sub) {if (item == ""){continue;} targets[item] = 1;}
    }
    debug_message("Checking targets - done");
    
    
    
    
    
	
    debug_message("Loading bams - start");
	// Loading BAM list
    v_message = " > Loading BAM files";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    filedb.clear();
    path_map = v_output + "/map/";
    for (const auto & entry : fs::directory_iterator(path_map))
    {
        string dir_path = entry.path();
        string sample = base_name(dir_path);
        string bamadj = dir_path + "/" + sample + ".adjusted.bam";
        string bamadjnodup = dir_path + "/" + sample + ".adjusted.nodup.bam";
 
        string bamunique = dir_path + "/" + sample + ".unique.bam";
        string bamuniquenodup = dir_path + "/" + sample + ".unique.nodup.bam";
 
 
        if (fileExists(bamadjnodup))
        {
            string cmd = v_samtools + " view -H " + bamadjnodup;
            string head = GetStdoutFromCommand(cmd.c_str());
            int validbam = 0;
            if (head.find("kir-mapper") != std::string::npos) {validbam = 1;}
            if (head.find("hla-mapper") != std::string::npos) {validbam = 1;}
            if (validbam != 1) {continue;}
            
            cmd = v_samtools + " samples " + bamadjnodup;
            head = GetStdoutFromCommand(cmd.c_str());
            vector <string> fields;
            string sample_id = "";
            if (head.size() != 0)
            {
                boost::split(fields,head,boost::is_any_of("\t"));
                sample_id = fields[0];
            }
            if (sample_id == "") {continue;}
			
			pair <string,string> k;
			k = make_pair(sample_id,"KIR3DL3");
			if (copynumber.find(k) == copynumber.end()) {continue;}
			
            if (failed_samples.find(sample_id) == failed_samples.end()) {filedb.push_back(bamadjnodup);}
            continue;
        }
 

 
        if ((! fileExists(bamadjnodup)) && (fileExists(bamadj)))
        {
            if (! fileExists(bamadj)) {continue;}

            string cmd = v_samtools + " view -H " + bamadj;
            string head = GetStdoutFromCommand(cmd.c_str());
            int validbam = 0;
            if (head.find("kir-mapper") != std::string::npos) {validbam = 1;}
            if (head.find("hla-mapper") != std::string::npos) {validbam = 1;}
            if (validbam != 1) {continue;}
            
            cmd = v_samtools + " samples " + bamadj;
            head = GetStdoutFromCommand(cmd.c_str());
            vector <string> fields;
            string sample_id = "";
            if (head.size() != 0)
            {
                boost::split(fields,head,boost::is_any_of("\t"));
                sample_id = fields[0];
            }
            if (sample_id == "") {continue;}
			
			pair <string,string> k;
			k = make_pair(sample_id,"KIR3DL3");
			if (copynumber.find(k) == copynumber.end()) {continue;}
			
            if (failed_samples.find(sample_id) == failed_samples.end()) {filedb.push_back(bamadj);}
            continue;
        }
		

		
		if ((! fileExists(bamadjnodup)) && (! fileExists(bamadj)))
        {
            if (! fileExists(bamunique)) {continue;}

			string cmd = v_samtools + " view -H " + bamunique;
            string head = GetStdoutFromCommand(cmd.c_str());
            int validbam = 0;
            if (head.find("kir-mapper") != std::string::npos) {validbam = 1;}
            if (head.find("hla-mapper") != std::string::npos) {validbam = 1;}
            if (validbam != 1) {continue;}
            
            cmd = v_samtools + " samples " + bamunique;
            head = GetStdoutFromCommand(cmd.c_str());
            vector <string> fields;
            string sample_id = "";
            if (head.size() != 0)
            {
                boost::split(fields,head,boost::is_any_of("\t"));
                sample_id = fields[0];
            }
            if (sample_id == "") {continue;}
			
			pair <string,string> k;
			k = make_pair(sample_id,"KIR3DL3");
			if (copynumber.find(k) == copynumber.end()) {continue;}
			
            if (failed_samples.find(sample_id) == failed_samples.end()) {filedb.push_back(bamunique);}
            continue;
        }
		
		
		
    }
    
	if (filedb.size() == 0)
    {
        v_message = "The list of BAM files is empty!";
        warnings.push_back (v_message);
        return;
    }
	
    v_message = " > Processing " + to_string(filedb.size()) + " BAM files with copy number data";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    debug_message("Loading bams - done");

    
    
    
    	
	
	
    
    


	
    for (auto &item : valid_genes)
	{
        string gene = item.first;
        if (gene == "") {continue;}

        if (v_target != "")
        {
            if (targets.find(gene) == targets.end()) {continue;}
        }
        
        if (gene == "5UPKIR") {continue;}
        if (gene == "HLA-E") {continue;}
        if (gene == "HLA-G") {continue;}

 

        debug_message("Creating output structure for this gene - start");

        cmd = "mkdir " + v_output_adjusted + gene;
        GetStdoutFromCommand(cmd.c_str());
        cmd = "mkdir " + v_output_adjusted + gene + "/reports/";
        GetStdoutFromCommand(cmd.c_str());
        cmd = "mkdir " + v_output_adjusted + gene + "/phasing/";
        GetStdoutFromCommand(cmd.c_str());
        cmd = "mkdir " + v_output_adjusted + gene + "/phasing/split_vcf";
        GetStdoutFromCommand(cmd.c_str());
        cmd = "mkdir " + v_output_adjusted + gene + "/vcf/";
        GetStdoutFromCommand(cmd.c_str());
        cmd = "mkdir " + v_output_adjusted + gene + "/vcf/log";
        GetStdoutFromCommand(cmd.c_str());
        cmd = "mkdir " + v_output_adjusted + gene + "/calls/";
        GetStdoutFromCommand(cmd.c_str());

       debug_message("Creating output structure for this gene - done");

        
        
        debug_message("Loading ignored positions - start");

		map <int,int> ignored_positions;
		string bed = v_resources + "/bed/" + gene + ".ignore.bed";
		
		ifstream ig (bed.c_str());
		for( std::string line; getline( ig, line ); )
		{
			if (line == "") {continue;}
			vector <string> data;
			boost::split(data,line,boost::is_any_of("\t"));
			int start = stoi(data[1]);
			int end = stoi(data[2]);
			for (int a = start; a <= end; a++)
			{
				ignored_positions[a] = 1;
			}
		}
		ig.close();
       debug_message("Loading ignored positions - done");

        
        	
	    debug_message("Loading exon positions - start");
    	map <int,int> exon_positions;
		bed = v_resources + "/bed/" + gene + ".CDS.bed";
		
		ifstream exon (bed.c_str());
		for( std::string line; getline( exon, line ); )
		{
			if (line == "") {continue;}
			vector <string> data;
			boost::split(data,line,boost::is_any_of("\t"));
			int start = stoi(data[1]);
			int end = stoi(data[2]);
			for (int a = (start - 1); a <= (end + 1); a++)
			{
				exon_positions[a] = 1;
			}
		}
		exon.close();
	    debug_message("Loading exon positions - done");
 
        
		
        
		int do_genotype = 1;
		string outvcfphased = v_output_adjusted + gene + "/vcf/" + gene + ".combined.trim.treated.norm.phased.vcf";
	
		if ((fileExists(outvcfphased)) && (filesize(outvcfphased.c_str()) > 0)) {
			if (v_update_call == 1) {do_genotype = 0;}
		}
		
		
		vector <string> head;
		vector <string> samples;
		map <int,string> snp_data;
		map <pair<int,string>,string> vcf_data;

		
        if (do_genotype == 1) {
        
			debug_message("Freebayes - start");
			
			v_message = " > " + gene + " : Running freebayes";
			screen_message (screen_size, 0, v_message, 1, v_quiet);
			
			vector <string> subregions;
			if (v_exome == 0) {
				for (auto sub : regions)
				{
					string subgene = sub.first.first;
					string region = sub.first.second;
					if (subgene != gene) {continue;}
					vector <string> subchr;
					boost::split(subchr,region,boost::is_any_of(":"));
					vector <string> limits;
					boost::split(limits,subchr[1],boost::is_any_of("-"));
					string subregions_txt = split_region(stoi(v_threads),subchr[0],stoi(limits[0]),stoi(limits[1]),20);
					boost::split(subregions,subregions_txt,boost::is_any_of("\n"));
					break;
				}
			}
        
			if (v_exome)
			{
				subregions.clear();
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
        
 
        
			string v_ref_gene_bam = "";
			if (v_exome == 0) {v_ref_gene_bam = v_resources + "/bam/" + gene + "_chr.bam";}
			if (v_exome == 1) {v_ref_gene_bam = v_resources + "/bam/" + gene + "_chr.cds.bam";}

			if (! fileExists(v_ref_gene_bam)) {continue;}
			string cmd = v_samtools + " samples " + v_ref_gene_bam;
			string sample_ref_list = GetStdoutFromCommand(cmd.c_str());
			vector <string> data;
			boost::split(data,sample_ref_list,boost::is_any_of("\n"));

        
        
			ofstream CP;
			string outcp = v_output_adjusted + gene + "/vcf/" + gene + "_copynumbers.txt";
			removefile(outcp,0);
			CP.open (outcp.c_str());

			for (auto item : copynumber)
			{
				if (item.first.second == gene)
				{
					CP << item.first.first << "\t" << item.second << endl;
				}
			}
			
			for (auto item : data)
			{
				if (item == "") {continue;}
				vector <string> sub;
				boost::split(sub,item,boost::is_any_of("\t"));
				if (sub[0] == "") {continue;}
				CP << sub[0] << "\t" << 1 << endl;
			}
			
			
			CP.close();
	 
        
        
        
			ThreadPool pool(stoi(v_threads));
			std::vector< std::future<int> > results;
			
			int call_count = 0;
			for (auto region : subregions)
			{
				call_count++;
				results.emplace_back(
					pool.enqueue([call_count,region,gene,filedb,&copynumber,outcp,v_ref_gene_bam,v_output_adjusted]
						{
							
							string current_region = "";
							current_region = region;
							boost::replace_all(current_region, "chrchr", "chr");
							
							string region_out = current_region;
							boost::replace_all(region_out, ":", ".");
							boost::replace_all(region_out, "-", ".");
							
							string vcfout = v_output_adjusted + gene + "/vcf/" + gene + "." + region_out + ".vcf";
							string vcflog = v_output_adjusted + gene + "/vcf/log/" + gene + "." + region_out + ".freebayes.log";
							string cmd = v_freebayes + "  --report-all-haplotype-alleles --use-best-n-alleles 8 --min-alternate-count 3 -A " + outcp + " -f " + v_ref + " -r " + current_region + " -v " + vcfout + " > " + vcflog;

					string samplelist = v_output_adjusted + gene + "/vcf/bamlist.txt";
					ofstream OUT;
					OUT.open (samplelist.c_str());
					
					for (auto item : filedb)
					{
						OUT << item << endl;
					}
					cmd.append (" -L " + samplelist);
					cmd.append (" -b " + v_ref_gene_bam);
					OUT.close();
					
					if (v_debug)
					{
						mtx_genotype.lock();
						debug_message(cmd);
						mtx_genotype.unlock();
					}
					
					
							string out = GetStdoutFromCommand(cmd.c_str());
							int failed = 0;
							if (filesize(vcfout.c_str()) == 0) {failed = 1;}
							if (out.find("Failed") != std::string::npos) {failed = 1;}
							if (out.find("Too many") != std::string::npos) {failed = 2;}
							if (failed == 1) {
								mtx_genotype.lock();
								cout << endl << "Freebayes failed. Quitting..." << endl;
								mtx_genotype.unlock();
								exit(0);
								return 1;
							}
							if (failed == 2) {
								mtx_genotype.lock();
								cout << endl << "Freebayes failed because of too many files." << endl;
								cout << "You should adjust 'ulimit -n' to a higher value. Quitting..." << endl;
								mtx_genotype.unlock();
								exit(0);
								return 1;
							}
					
							return 1;
						})
				);

			}
		
			for(auto && result: results){result.get();} // waiting for all threads
			debug_message("Freebayes - done");
 
		
			if (gene == "KIR2DS4") {
				debug_message("22bp KIR2DS4 - start");

				map <string,float> count_del;
				map <string,float> count_ins;
				map <string,int> samples;
				v_command = v_samtools + " view " + v_ref_gene_bam;
				string out = GetStdoutFromCommand(v_command);
				vector <string> lines;
				boost::split(lines,out,boost::is_any_of("\n"));
				for (auto item : lines)
				{
					if (item == "") {continue;}
					if (item.substr(0,1) == "[") {continue;}
					vector <string> fields;
					boost::split(fields,item,boost::is_any_of("\t"));
					vector <string> allele;
					boost::split(allele,fields[0],boost::is_any_of("."));
					string ref = "ref." + gene + "." + allele[1];
					samples[ref] = 1;
					string seq = "";
					if (fields.size() >= 9) {
						seq = fields[9];
					}
					if (seq.find("CCCGGAGCTCCTATGACATGTACCAT") != std::string::npos) {count_ins[ref]++;}
					if (seq.find("ATGGTACATGTCATAGGAGCTCCGGG") != std::string::npos) {count_ins[ref]++;}
					if (seq.find("TTGTCCTGCAGCTCCATCTATCCAGGG") != std::string::npos) {count_del[ref]++;}
					if (seq.find("CCCTGGATAGATGGAGCTGCAGGACAA") != std::string::npos) {count_del[ref]++;}
				}

				string samplelist = v_output_adjusted + gene + "/vcf/bamlist.txt";
				
				ifstream list (samplelist.c_str());
				for( std::string line; getline(list, line ); )
				{
					if (line == "") {continue;}
					string ref = base_name(line);
					boost::replace_all(ref, ".adjusted.nodup.bam", "");
					boost::replace_all(ref, ".adjusted.bam", "");
					boost::replace_all(ref, ".unique.bam", "");
					boost::replace_all(ref, ".unique.nodup.bam", "");
					samples[ref] = 1;
					
					v_command = v_samtools + " view " + line + " chr19:54839300-54839800";
					string out = GetStdoutFromCommand(v_command);
					vector <string> lines;
					boost::split(lines,out,boost::is_any_of("\n"));
					for (auto item : lines)
					{
						if (item == "") {continue;}
						if (item.substr(0,1) == "[") {continue;}
						vector <string> fields;
						boost::split(fields,item,boost::is_any_of("\t"));
						int flag = stoi(fields[1]);
						if (flag > 200) {continue;}
						string seq = "";
						if (fields.size() >= 9) {
							seq = fields[9];
						}
						if (seq.find("CCCGGAGCTCCTATGACATGTACCAT") != std::string::npos) {count_ins[ref]++;}
						if (seq.find("ATGGTACATGTCATAGGAGCTCCGGG") != std::string::npos) {count_ins[ref]++;}
						if (seq.find("TTGTCCTGCAGCTCCATCTATCCAGGG") != std::string::npos) {count_del[ref]++;}
						if (seq.find("CCCTGGATAGATGGAGCTGCAGGACAA") != std::string::npos) {count_del[ref]++;}
					}
					
				}
				list.close();
				
				string outvcf = v_output_adjusted + "/" + gene + "/vcf/" + gene + ".22bpIns.vcf";
				ofstream VCF;
				VCF.open (outvcf.c_str());
				
				VCF << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
				for (auto loop : samples)
				{
					VCF << "\t" << loop.first;
				}
				VCF << endl;


				VCF << "chr19\t54839510\t<22bpIns>\tC\tCCGGAGCTCCTATGACATGTACC\t.\t.\t.\tGT:DP:AD:PS";
				for (auto loop : samples)
				{
					float del = count_del[loop.first];
					float ins = count_ins[loop.first];
					float cov = del + ins;
					string depth = to_string(int(del)) + "," + to_string(int(ins));
					
					float ratio_ins = ins/cov;
					float ratio_del = del/cov;

					pair <string,string> k;
					k = make_pair(loop.first,gene);
					boost::replace_all(data[2], "+", "");
					string cn = copynumber[k];
					
					if (loop.first.find("ref.") != std::string::npos) {cn = "1";}
					
					string gen = ".";
					if (cn != "0")
					{
						if (cn == "1")
						{
							if (ratio_ins >= 0.9){gen = "1";}
							if (ratio_del >= 0.9){gen = "0";}
							if ((ratio_del >= 0.35) && (ratio_ins >= 0.35)){gen = ".";}
						}

						if (cn == "2")
						{
							if (ratio_ins >= 0.9){gen = "1/1";}
							if (ratio_del >= 0.9){gen = "0/0";}
							if ((ratio_del >= 0.25) && (ratio_ins >= 0.25)){gen = "0/1";}
						}
						
						if (cn == "3")
						{
							if (ratio_ins >= 0.9){gen = "1/1/1";}
							if (ratio_del >= 0.9){gen = "0/0/0";}
							if ((ratio_del >= 0.25) && (ratio_ins >= 0.25)){gen = "0/1/.";}
						}
						if (cn == "4")
						{
							if (ratio_ins >= 0.9){gen = "1/1/1/1";}
							if (ratio_del >= 0.9){gen = "0/0/0/0";}
							if ((ratio_del >= 0.25) && (ratio_ins >= 0.25)){gen = "0/1/./.";}
						}
						
					}

					VCF << "\t" + gen + ":" + to_string(int(cov)) + ":" + depth + ":.";
				}
				VCF << endl;

				VCF.close();
				debug_message("22bp KIR2DS4 - end");

			}

        
        
			if ((gene == "KIR3DL3")  && (v_full == 1)){
				debug_message("20bp KIR3DL3 - start");

				map <string,float> count_del;
				map <string,float> count_ins;
				map <string,int> samples;
				v_command = v_samtools + " view " + v_ref_gene_bam;
				string out = GetStdoutFromCommand(v_command);
				vector <string> lines;
				boost::split(lines,out,boost::is_any_of("\n"));
				for (auto item : lines)
				{
					if (item == "") {continue;}
					if (item.substr(0,1) == "[") {continue;}
					vector <string> fields;
					boost::split(fields,item,boost::is_any_of("\t"));
					vector <string> allele;
					boost::split(allele,fields[0],boost::is_any_of("."));
					string ref = "ref." + gene + "." + allele[1];
					samples[ref] = 1;
					string seq = "";
					if (fields.size() >= 9) {
						seq = fields[9];
					}
					if (seq.find("TGCTATTCCACCTTTCCTCAGAGTATCTT") != std::string::npos) {count_del[ref]++;}

					if (seq.find("TGCTATTCCACCTTTCCTCATGTTGTTCC") != std::string::npos) {count_ins[ref]++;}
				}

				string samplelist = v_output_adjusted + gene + "/vcf/bamlist.txt";
				
				ifstream list (samplelist.c_str());
				for( std::string line; getline(list, line ); )
				{
					if (line == "") {continue;}
					string ref = base_name(line);
					boost::replace_all(ref, ".adjusted.nodup.bam", "");
					boost::replace_all(ref, ".adjusted.bam", "");
					boost::replace_all(ref, ".unique.bam", "");
					boost::replace_all(ref, ".unique.nodup.bam", "");
					samples[ref] = 1;
					
					v_command = v_samtools + " view " + line + " chr19:54736020-54736920";
					string out = GetStdoutFromCommand(v_command);
					vector <string> lines;
					boost::split(lines,out,boost::is_any_of("\n"));
					for (auto item : lines)
					{
						if (item == "") {continue;}
						if (item.substr(0,1) == "[") {continue;}
						vector <string> fields;
						boost::split(fields,item,boost::is_any_of("\t"));
						int flag = stoi(fields[1]);
						if (flag > 200) {continue;}
						string seq = "";
						if (fields.size() >= 9) {
							seq = fields[9];
						}
						if (seq.find("TGCTATTCCACCTTTCCTCAGAGTATCTT") != std::string::npos) {count_del[ref]++;}
						if (seq.find("TGCTATTCCACCTTTCCTCATGTTGTTCC") != std::string::npos) {count_ins[ref]++;}
					}
					
				}
				list.close();
				
				string outvcf = v_output_adjusted + "/" + gene + "/vcf/" + gene + ".20bpIns.vcf";
				ofstream VCF;
				VCF.open (outvcf.c_str());
				
				VCF << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
				for (auto loop : samples)
				{
					VCF << "\t" << loop.first;
				}
				VCF << endl;


				VCF << "chr19\t54736520\t<20bpDel>\tA\tATTCCACCTTTCCTCATGTTG\t.\t.\t.\tGT:DP:AD:PS";
				for (auto loop : samples)
				{
					float del = count_del[loop.first];
					float ins = count_ins[loop.first];
					float cov = del + ins;
					string depth = to_string(int(del)) + "," + to_string(int(ins));
					
					float ratio_ins = ins/cov;
					float ratio_del = del/cov;

					pair <string,string> k;
					k = make_pair(loop.first,gene);
					boost::replace_all(data[2], "+", "");
					string cn = copynumber[k];
					
					if (loop.first.find("ref.") != std::string::npos) {cn = "1";}
					
					string gen = ".";
					if (cn != "0")
					{
						if (cn == "1")
						{
							if (ratio_ins >= 0.9){gen = "1";}
							if (ratio_del >= 0.9){gen = "0";}
							if ((ratio_del >= 0.35) && (ratio_ins >= 0.35)){gen = ".";}
						}

						if (cn == "2")
						{
							if (ratio_ins >= 0.9){gen = "1/1";}
							if (ratio_del >= 0.9){gen = "0/0";}
							if ((ratio_del >= 0.25) && (ratio_ins >= 0.25)){gen = "0/1";}
						}
						
						if (cn == "3")
						{
							if (ratio_ins >= 0.9){gen = "1/1/1";}
							if (ratio_del >= 0.9){gen = "0/0/0";}
							if ((ratio_del >= 0.25) && (ratio_ins >= 0.25)){gen = "0/1/.";}
						}
						if (cn == "4")
						{
							if (ratio_ins >= 0.9){gen = "1/1/1/1";}
							if (ratio_del >= 0.9){gen = "0/0/0/0";}
							if ((ratio_del >= 0.25) && (ratio_ins >= 0.25)){gen = "0/1/./.";}
						}
						
					}

					VCF << "\t" + gen + ":" + to_string(int(cov)) + ":" + depth + ":.";
				}
				VCF << endl;

				VCF.close();
				debug_message("20bp KIR3DL3 - end");

			}

        



		

			debug_message("Reading VCFs - start");

			v_message = " > " + gene + " : Processing VCF";
			screen_message (screen_size, 0, v_message, 1, v_quiet);


			int count = 0;
			string sampleline = "";

			for (auto region : subregions)
			{
				string current_region = "";
				current_region = region;
				boost::replace_all(current_region, "chrchr", "chr");
				
				string region_out = current_region;
				boost::replace_all(region_out, ":", ".");
				boost::replace_all(region_out, "-", ".");

				string vcfout = v_output_adjusted + gene + "/vcf/" + gene + "." + region_out + ".vcf";
				debug_message("reading " + vcfout);
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
					
					if (ignored_positions.find(stoi(data[1])) != ignored_positions.end()) {continue;}
					
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
				debug_message("closing " + vcfout);
				count++;
			}
			debug_message("Reading VCFs done ");
			
        
        
        
        
			if (gene == "KIR2DS4")
			{
				string outvcf = v_output_adjusted + "/" + gene + "/vcf/" + gene + ".22bpIns.vcf";
				ifstream vcffile (outvcf.c_str());
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
					
					if (ignored_positions.find(stoi(data[1])) != ignored_positions.end()) {continue;}
					
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
			}
			
			if (gene == "KIR3DL3")
			{
				string outvcf = v_output_adjusted + "/" + gene + "/vcf/" + gene + ".20bpIns.vcf";
				ifstream vcffile (outvcf.c_str());
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
					
					if (ignored_positions.find(stoi(data[1])) != ignored_positions.end()) {continue;}
					
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
			}


        
			debug_message("Write combined VCF - start ");
			string outvcf = v_output_adjusted + gene + "/vcf/" + gene + ".combined.vcf";
			ofstream COMBVCF;
			COMBVCF.open (outvcf.c_str());
			string logvcf = v_output_adjusted + gene + "/vcf/" + gene + ".combined.vcf.log";
			ofstream LOG;
			LOG.open (logvcf.c_str());
			LOG << "SAMPLE\tPOS\tORIGINAL_GENOTYPE\tNEW_GENOTYPE" << endl;
			
			for (auto item : head)
			{
				COMBVCF << item << endl;
			}
			COMBVCF << "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">" << endl;
			debug_message("Write combined VCF - head ");
			COMBVCF << sampleline << endl;
			debug_message("Write combined VCF - sample line ");

			for (auto item : snp_data)
			{
				int pos = item.first;
				string info = item.second;
				COMBVCF << info;
				debug_message("Write combined VCF - position " + to_string(pos));

				for (int a = 9; a < samples.size(); a++)
				{
					
					string sample = samples[a];
					debug_message("Write combined VCF - position " + to_string(pos) + " sample " + sample);
					pair <int,string> k;
					k = make_pair(pos,sample);
					string value = vcf_data[k];
					
					if (sample.find("ref.") != std::string::npos) {
						COMBVCF << "\t" << value;
						continue;
					}
					
					
					pair <string,string> j;
					j = make_pair(sample,gene);
					if (copynumber[j] == "0")  {
						COMBVCF << "\t" << value;
						continue;
					}
					
						debug_message("Write combined VCF - treating genotype start: " + value);
						if (value == "") {value = ".:.:.:.:.:.:.";}
						string newgenotype = treat_genotype(value,stoi(copynumber[j]));
						debug_message("Write combined VCF - treating genotype done" + value);
						if (newgenotype != value)
						{
							string logvalue = sample + "\t" + to_string(pos) + "\t" + value + "\t" + newgenotype;
							LOG << logvalue << endl;
						}

						COMBVCF << "\t" << newgenotype;
					continue;

				}
				COMBVCF << endl;
				
			}
			COMBVCF.close();
			LOG.close();
			debug_message("Write combined VCF - done");
		

        
			debug_message("Trimming VCF - start");
			string outvcftrim = v_output_adjusted + gene + "/vcf/" + gene + ".combined.trim.vcf";
			cmd = v_bcftools + " view --trim-alt-alleles " + outvcf + " -o " + outvcftrim;
			GetStdoutFromCommand(cmd.c_str());
			debug_message("Trimming VCF - done");


       
        
  

			debug_message("Treating VCF - start");
			map <string,int> must_keep;
			vector <string> samples_list;

			ifstream vcftmp (outvcftrim.c_str());
			for( std::string line; getline(vcftmp, line ); )
			{
				if (line == "") {continue;}
				if (line.substr(0,2) == "##")
				{
					continue;
				}
				if (line.substr(0,2) == "#C")
				{
					boost::split(samples_list,line,boost::is_any_of("\t"));
					continue;
				}
				vector <string> data;
				boost::split(data,line,boost::is_any_of("\t"));

				for (int a = 8; a < data.size(); a++){
					if (samples_list[a].find("ref.") == std::string::npos){continue;}
					vector <string> fields;
					boost::split(fields,data[a],boost::is_any_of(":"));
					if (fields[0] == ".") {continue;}
					if (fields[0] == "1"){must_keep[data[1]] = 1;}
					if (fields[0] == "2"){must_keep[data[1]] = 1;}
					if (fields[0] == "3"){must_keep[data[1]] = 1;}
					if (fields[0] == "4"){must_keep[data[1]] = 1;}
					if (fields[0] == "5"){must_keep[data[1]] = 1;}
					if (fields[0] == "6"){must_keep[data[1]] = 1;}
	 }
			}
			vcftmp.close();
    
  
        
        
			string outvcftrimtreated = v_output_adjusted + gene + "/vcf/" + gene + ".combined.trim.treated.vcf";
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

				if (must_keep.find(data[1]) == must_keep.end()) {
					if (data[4] == ".") {continue;}
					if (data[4] == "*") {continue;}
					if (data[5] == "0") {continue;}
				}
				TREATEDVCF << line << endl;
			}
			
			TREATEDVCF.close();
			vcf.close();
			
			debug_message("Treating VCF - done");


        
			string outvcftrimtreatednorm = v_output_adjusted + gene + "/vcf/" + gene + ".combined.trim.treated.norm.vcf";
			debug_message("VCF Normalization - start");
//			string ref = v_db + "/reference/hg38/" + chr_hg38[gene] + ".fasta";
			string ref = v_db + "/reference/hg38/reference.fasta";
			cmd = v_bcftools + " norm -f " + ref + " -o " + outvcftrimtreatednorm + " " + outvcftrimtreated;
			GetStdoutFromCommand(cmd.c_str());
			debug_message("VCF Normalization - done");
    



// reloading
			debug_message("Reloading VCF - start");

			head.clear();
			samples.clear();
			snp_data.clear();
			vcf_data.clear();
			sampleline = "";
			

			ifstream vcffile (outvcftrimtreatednorm.c_str());
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
				
				if (ignored_positions.find(stoi(data[1])) != ignored_positions.end()) {continue;}
				
				snp_data[stoi(data[1])] = data[0] + "\t" + data[1] + "\t" + data[2] + "\t" + data[3] + "\t" + data[4] + "\t" + data[5] + "\t" + data[6] + "\t" + data[7] + "\t" + data[8];
				
				for (int a = 9; a < data.size(); a++)
				{
					pair <int,string> k;
					k = make_pair(stoi(data[1]),samples[a]);
					vcf_data[k] = data[a];
				}
			}
			vcffile.close();

			debug_message("Reloading VCF - done");


        
        


			debug_message("Whatshap - start");
			v_message = " > " + gene + " : Phasing with WhatsHap";
			screen_message (screen_size, 0, v_message, 1, v_quiet);

			
			vector <string> phasing_cmds;
			for (int a = 9; a < samples.size(); a++)
			{
				string sample = samples[a];
				if (sample.find("ref.") != std::string::npos) {continue;}
				
				pair <string,string> k;
				k = make_pair(sample,gene);
				if (copynumber[k] == "0") {continue;}
				if (copynumber[k] == "1") {continue;}
	//            if (copynumber[k] == "3") {continue;}
				if (copynumber[k] == "4") {continue;}
				
				string vcfout = v_output_adjusted + gene + "/phasing/split_vcf/" + sample + ".vcf";
				string vcfouttrim = v_output_adjusted + gene + "/phasing/split_vcf/" + sample + ".trim.vcf";
				string vcfoutwhats = v_output_adjusted + gene + "/phasing/split_vcf/" + sample + ".trim.whats.vcf";
				string cmd = v_bcftools + " view -s " + sample + " " + outvcftrimtreatednorm + " -o " + vcfout;
				GetStdoutFromCommand(cmd.c_str());
				cmd = v_bcftools + " view --trim-alt-alleles --min-ac 1 " + vcfout + " -o " + vcfouttrim;
				GetStdoutFromCommand(cmd.c_str());
				
				string bam = "";
				for (auto item: filedb)
				{
					if (item.find(sample) != std::string::npos) {bam = item;}
				}
				
				if (copynumber[k] == "2")
				{
					cmd = v_whats + " phase --indels -r " + v_ref + " -o " + vcfoutwhats + " " + vcfouttrim + " " + bam;
					phasing_cmds.push_back(cmd);
				}

				if ((copynumber[k] == "3") && (v_nopolyphasing == 0))
				{
					cmd = v_whats + " polyphase -p 3 -r " + v_ref + " -o " + vcfoutwhats + " " + vcfouttrim + " " + bam;
					phasing_cmds.push_back(cmd);
				}
				
			}
			
			ThreadPool poolwhats(stoi(v_threads));
			std::vector< std::future<int> > results_whats;
			
			int run_count = 0;
			for (auto item : phasing_cmds)
			{
				call_count++;
				results_whats.emplace_back(
					poolwhats.enqueue([run_count,item]
						{
							GetStdoutFromCommand(item.c_str());
							return 1;
						})
				);
				
			}
			
			for(auto && result: results_whats){result.get();} // waiting for all threads
			debug_message("Whatshap - end");

			results_whats.clear();
			
			
			debug_message("Whatshap - reload VCF - start");

			map <pair<int,string>,string> original_allele_code;
			for (auto item : snp_data)
			{
				int pos = item.first;
				string info = item.second;
				vector <string> data;
				boost::split(data,info,boost::is_any_of("\t"));
				vector <string> alleles;
				string subtext = data[3] + "," + data[4];
				boost::split(alleles,subtext,boost::is_any_of(","));
				for (int a = 0; a < alleles.size();a++)
				{
					pair <int,string> k;
					k = make_pair(pos,alleles[a]);
					original_allele_code[k] = to_string(a);
				}
			}
			debug_message("Whatshap - reload VCF - done");
		
		

  
			debug_message("Writing new VCF file after whatshap - start");

			vector <string> subsamples;
			for (int a = 9; a < samples.size(); a++)
			{
				string sample = samples[a];
				
				pair <string,string> k;
				k = make_pair(sample,gene);
				if (copynumber[k] == "0") {continue;}
				if (copynumber[k] == "1") {continue;}
				
				string vcfoutwhats = v_output_adjusted + gene + "/phasing/split_vcf/" + sample + ".trim.whats.vcf";
				ifstream vcffile (vcfoutwhats.c_str());
				vector <string> subsamples;
				for( std::string line; getline( vcffile, line ); )
				{
					if (line == "") {continue;}
					if (line.substr(0,2) == "##") {continue;}
					if (line.substr(0,2) == "#C") {boost::split(subsamples,line,boost::is_any_of("\t")); continue;}
					vector <string> data;
					boost::split(data,line,boost::is_any_of("\t"));
					int pos = stoi(data[1]);
					if (ignored_positions.find(pos) != ignored_positions.end()) {continue;}
					
					vector <string> alleles;
					string subtext = data[3] + "," + data[4];
					boost::split(alleles,subtext,boost::is_any_of(","));
					
					for (int a = 9; a < data.size(); a++)
					{
						
						pair <int,string> v;
						v = make_pair(pos,subsamples[a]);
						string originaldata = vcf_data[v];
						vector <string> original_fields;
						boost::split(original_fields,originaldata,boost::is_any_of(":"));
						
						
						string info = data[a];
						vector <string> fields;
						boost::split(fields,info,boost::is_any_of(":"));
						vector <string> suballeles;
						boost::split(suballeles,fields[0],boost::is_any_of("/"));
						if (suballeles.size() == 1) {boost::split(suballeles,fields[0],boost::is_any_of("|"));}
						string currentA = "";
						if (suballeles[0] != ".") {
							currentA = alleles[stoi(suballeles[0])];
						}
						string currentB = "";
						if (suballeles[1] != ".") {
							currentB = alleles[stoi(suballeles[1])];
						}
						
						string currentC = ".";
						string currentD = ".";
						
						if ((copynumber[k] == "3") || (copynumber[k] == "4")){
							if (suballeles[2] != ".") {
								currentC = alleles[stoi(suballeles[2])];
							}
						}
						
						if (copynumber[k] == "4"){
							if (suballeles[3] != ".") {
								currentD = alleles[stoi(suballeles[2])];
							}
						}
						
						
						
						int phased = 0;
						if (info.find("|") != std::string::npos) {phased = 1;}
						
						
						pair <int,string> j;
						j = make_pair(pos,currentA);
						string oldA = original_allele_code[j];
						if (currentA == ".") {oldA = ".";}
						if (currentA == "") {oldA = ".";}
						
						pair <int,string> h;
						h = make_pair(pos,currentB);
						string oldB = original_allele_code[h];
						if (currentB == ".") {oldB = ".";}
						if (currentB == "") {oldB = ".";}
						
						pair <int,string> i;
						i = make_pair(pos,currentC);
						string oldC = original_allele_code[i];
						if (currentC == ".") {oldC = ".";}
						if (currentC == "") {oldC = ".";}

						pair <int,string> l;
						l = make_pair(pos,currentD);
						string oldD = original_allele_code[l];
						if (currentD == ".") {oldD = ".";}
						if (currentD == "") {oldD = ".";}
						
						
						string newinfo = oldA;
						if (phased == 0) {newinfo.append("/");}
						if (phased == 1) {newinfo.append("|");}
						newinfo.append(oldB);
						
						if ((copynumber[k] == "3") || (copynumber[k] == "4")){
							if (phased == 0) {newinfo.append("/");}
							if (phased == 1) {newinfo.append("|");}
							newinfo.append(oldC);
						}
						
						if (copynumber[k] == "4"){
							if (phased == 0) {newinfo.append("/");}
							if (phased == 1) {newinfo.append("|");}
							newinfo.append(oldD);
						}
						
						newinfo.append(":" + original_fields[1]);
						newinfo.append(":" + original_fields[2]);
						newinfo.append(":" + fields[3]);
						
						pair <int,string> k;
						k = make_pair(pos,subsamples[a]);
						vcf_data[k] = newinfo;
					}
				}
				vcffile.close();
			}
		


 		
			string outvcfphased = v_output_adjusted + gene + "/vcf/" + gene + ".combined.trim.treated.norm.phased.vcf";
			ofstream phasedVCF;
			phasedVCF.open (outvcfphased.c_str());
			for (auto item : head)
			{
				phasedVCF << item << endl;
			}
			phasedVCF << sampleline << endl;
			
			for (auto item : snp_data)
			{
				int pos = item.first;
				string info = item.second;
				phasedVCF << info;
				
				for (int a = 9; a < samples.size(); a++)
				{
					string sample = samples[a];
					pair <int,string> k;
					k = make_pair(pos,sample);
					string value = vcf_data[k];
					phasedVCF << "\t" << value;
				}
				phasedVCF << endl;
			}
			phasedVCF.close();
			

                
		
			if (gene == "KIR3DP1")
			{
				cmd = "mv " + outvcfphased + " " + outvcfphased + ".tmp";
				GetStdoutFromCommand(cmd.c_str());
				string perlfile = v_resources + "/scripts/KIR3DP1_E2.pl";
				cmd = "perl " + perlfile + " " + outvcfphased + ".tmp > " + outvcfphased;
				GetStdoutFromCommand(cmd.c_str());
				
				
				ifstream newfile;
				newfile.open(outvcfphased);
				
				for( std::string line; getline( newfile, line ); )
				{
					if (line == "") {continue;}
					if (line.substr(0,2) == "##"){continue;}
					if (line.substr(0,2) == "#C"){continue;}
					vector <string> data;
					boost::split(data,line,boost::is_any_of("\t"));
				
					if (data[1] == "63062")
					{
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
				}
				newfile.close();
				removefile(outvcfphased + ".tmp",v_debug);
			}

		
        
			debug_message("Writing new VCF file after whatshap - end");

        
		} // do while not updating the VCF

        
        if (v_exome == 1) {
            v_message = " > " + gene + " : Checking SNP compatibility for CDS";
            screen_message (screen_size, 0, v_message, 0, v_quiet);
        }
        else {
            v_message = " > " + gene + " : Checking full SNP compatibility";
            screen_message (screen_size, 0, v_message, 0, v_quiet);
        }

        
        int call_higher_cn = 1;

        string nullallele = gene + "*null";
        string unresolved = gene + "*unresolved";
        
        head.clear();
        samples.clear();
        snp_data.clear();
        map <pair<int,string>,string> sample_snp;
        map <pair<string,string>,string> sample_ps;
        map <pair<int,string>,string> ref_data;
        map <string,int> ref_names;
        map <string,int> sample_names;
        map <pair<string,string>,string> ps_list;
        map <string,string> ps_positions;
        string sampleline = "";
        string chr;
        
        
        head.clear();
        samples.clear();
        snp_data.clear();
        ref_data.clear();
        ps_list.clear();
        
        sampleline = "";
 

        
        debug_message("Loading data for SNP comparison...start");
        ifstream vcffilein (outvcfphased.c_str());
        for( std::string line; getline( vcffilein, line ); )
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
            
            
            if (ignored_positions.find(stoi(data[1])) != ignored_positions.end()) {continue;}

            
            snp_data[stoi(data[1])] = data[0] + "\t" + data[1] + "\t" + data[2] + "\t" + data[3] + "\t" + data[4] + "\t" + data[5] + "\t" + data[6] + "\t" + data[7] + "\t" + data[8];
            if (chr == "") {chr = data[0];}
            
            for (int a = 9; a < data.size(); a++)
            {
                pair <int,string> k;
                k = make_pair(stoi(data[1]),samples[a]);
                
                vector <string> fields;
                boost::split(fields,data[a],boost::is_any_of(":"));
                
                pair <string,string> j;
                j = make_pair(samples[a],fields[3]);
                
                if (samples[a].find("ref.") != std::string::npos) {
                    ref_data[k] = fields[0];
                    ref_names[samples[a]] = 1;
                }
                else {
                    sample_snp[k] = fields[0];
                    sample_ps[j].append("," + fields[3]);
                    sample_names[samples[a]] = 1;
                    
                    string ps = fields[3];
                    if (fields[3] == ".") {ps = data[1];}
                    pair <string,string> l;
                    l = make_pair(samples[a],ps);
                    ps_list[l].append("," + data[1]);
                    ps_positions[samples[a]].append(","+ps);
                }
            }
        }
        vcffilein.close();
        debug_message("Loading data for SNP comparison...end");
        
        debug_message("Loading combinations...start");
        vector <string> main_combination_list;
        map <string,int> used;
        for (auto itemA : ref_names)
        {
            for (auto itemB : ref_names)
            {
                if (used.find(itemB.first) != used.end()) {continue;}
                string comb = itemA.first + "\t" + itemB.first;
                main_combination_list.push_back(comb);
            }
            used[itemA.first] = 1;
        }
        
        
        vector <string> main_combination_list_3cp;
        if (highcp[gene] > 2) {
            string combfile = "";
            if (v_exome == 1) {combfile = v_resources + "/lists/" + gene + ".cds.3alleles.txt";}
            if (v_exome == 0) {combfile = v_resources + "/lists/" + gene + ".full.3alleles.txt";}
            if (! fileExists(combfile)) {call_higher_cn = 0;}
            if (fileExists(combfile)) {
                ifstream comb (combfile .c_str());
                for( std::string line; getline( comb, line ); )
                {
                    if (line == "") {continue;}
                    main_combination_list_3cp.push_back(line);
                }
                comb.close();
            }
        }
        
        vector <string> main_combination_list_4cp;

        
        debug_message("Loading combinations...end");
        
        map <string,string> final_results;
            

        ThreadPool poolsamples(stoi(v_threads));
        std::vector< std::future<int> > results_genotype;
        int loop = 0;
        float done = 0;



        vector <string> ordered_samples;
        for (auto item : sample_names)
        {
            string sample = item.first;
            pair <string,string> k;
            k = make_pair(sample,gene);
            string cn = copynumber[k];
            if (cn == "0") {ordered_samples.push_back(sample);}
        }
        for (auto item : sample_names)
        {
            string sample = item.first;
            pair <string,string> k;
            k = make_pair(sample,gene);
            string cn = copynumber[k];
            if (cn == "1") {ordered_samples.push_back(sample);}
        }
        for (auto item : sample_names)
        {
            string sample = item.first;
            pair <string,string> k;
            k = make_pair(sample,gene);
            string cn = copynumber[k];
            if (cn == "2") {ordered_samples.push_back(sample);}
        }
        for (auto item : sample_names)
        {
            string sample = item.first;
            pair <string,string> k;
            k = make_pair(sample,gene);
            string cn = copynumber[k];
            if (cn == "3") {ordered_samples.push_back(sample);}
        }
        for (auto item : sample_names)
        {
            string sample = item.first;
            pair <string,string> k;
            k = make_pair(sample,gene);
            string cn = copynumber[k];
            if (cn == "4") {ordered_samples.push_back(sample);}
        }


        for (auto item : ordered_samples)
//        for (auto item : sample_names)
        {
//            string sample = item.first;
            string sample = item;
            loop++;
            
            results_genotype.emplace_back(
                poolsamples.enqueue([loop,sample, &copynumber, gene, nullallele, ref_names, &ps_list, &ps_positions, &sample_snp, &ref_data, &main_combination_list, &main_combination_list_3cp, &main_combination_list_4cp, v_output_adjusted, call_higher_cn, unresolved, &final_results, chr, &done, &filedb,&snp_data]
                    {
 
            string check_out_file = v_output_adjusted + gene + "/reports/" + sample + "." + gene + ".txt";
                
            pair <string,string> k;
            k = make_pair(sample,gene);
            string cn = copynumber[k];
                
            mtx_genotype.lock();
            debug_message("Starting sample " + sample);
            mtx_genotype.unlock();
            
                
            vector <string> selected_alleles;
            map <int,string> allele_dif;
            if ((cn != "0") && (cn != "4")) {
                for (auto refs : ref_names)
                {
                    string ref = refs.first;
                    int tested = 0;
                    int valid = 0;
                    
                    for (auto data : snp_data)
                    {
                        int pos = data.first;

                        pair <int,string> k;
                        k = make_pair(pos,ref);
                        string reference = ref_data[k];
                        if (reference == ".") {continue;}
                        k = make_pair(pos,sample);
                        string genotype = sample_snp[k];
                        boost::replace_all(genotype, "|", "/");
                        vector <string> alleles;
                        boost::split(alleles,genotype,boost::is_any_of("/"));
                        
                        float ok = 0;
                        float miss = 0;
                        for (auto test : alleles)
                        {
                            if (test == reference) {ok++;}
                            if (test == ".") {miss++;}
                        }
                        if (ok > 0) {tested++; valid++; continue;}
                        if ((ok == 0) && (miss > 0)){tested++; valid++; continue;}
                        tested++;
                    }
                    int dif = tested - valid;
                    allele_dif[dif].append("," + ref);
                }
                
                for (auto data : allele_dif)
                {
                    string list = data.second;
                    vector <string> sub;
                    boost::split(sub,list,boost::is_any_of(","));
                    for (auto item : sub)
                    {
                        if (item == "") {continue;}
                        selected_alleles.push_back(item);
                    }
                    if (selected_alleles.size() >= 25) {break;}
                }
            }
                
            debug_message("End allele selection for " + sample);

            map <float,string> results_ratio;
            map <float,string> results_dif;
            
            if (cn == "0") {
                string res = sample + "\t" + cn + "\t" + chr + "\t" + nullallele + "\t" + nullallele + "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
                results_ratio[1] = res;
                results_dif[0] = res;
            }

            
            if (cn == "4") {
                string res = sample + "\t" + cn + "\t" + chr + "\t" + unresolved + "\t" + unresolved + "\t" + unresolved + "\t" + unresolved + "\tNA\tNA\tNA\tNA\tNA\tNA";
                results_ratio[1] = res;
                results_dif[0] = res;
            }
            
            if (cn == "1")
            {
                
                for (auto comb : ref_names)
                {
                    string alleleA = comb.first;
                    
                    if (std::find(selected_alleles.begin(), selected_alleles.end(), alleleA) == selected_alleles.end())
                    {
                        continue;
                    }
                    
                    float tested = 0;
                    float miss = 0;
                    float valid = 0;
                    string miss_list = "";
                    string error_list = "";
                    

                    string ps_positions_str = ps_positions[sample];
                    map <string,int> ps_positions_map_nonredunt;
                    vector <string> tmp;
                    boost::split(tmp,ps_positions_str,boost::is_any_of(","));

                    for (auto item : tmp)
                    {
                        if (item == "") {continue;}
                        ps_positions_map_nonredunt[item] = 1;
                    }
                    ps_positions_str = "";
                    for (auto item : ps_positions_map_nonredunt)          
                    {
                        ps_positions_str.append("," + item.first);
                    }
                    ps_positions_str = ps_positions_str.substr(1);          
                    
 //                   string ps_positions_str = ps_positions[sample];
 //                   ps_positions_str = ps_positions_str.substr(1);
                    vector <string> ps_positions_list;
                    boost::split(ps_positions_list,ps_positions_str,boost::is_any_of(","));
                    
                    for (auto ps : ps_positions_list)
                    //for (auto ps : ps_list)
                    {
                        string subsample = sample;
                        //string subsample = ps.first.first;
                        //if (subsample != sample) {continue;}
                        //string pos = ps.second;
                        
                        pair <string,string> t;
                        t = make_pair(subsample,ps);
                        string pos = ps_list[t];

                        vector <string> positions;
                        string substringpos = pos.substr(1);
                        boost::split(positions,substringpos,boost::is_any_of(","));

                        
                        string h1s = "";
                        string h1r = "";
                        for (auto linked : positions)
                        {
                            int position = stoi(linked);
                            pair <int,string> k;
                            k = make_pair(position,sample);
                            string snp = sample_snp[k];
                            h1s.append("," + snp);
                            
                            k = make_pair(position,alleleA);
                            string nucA = ref_data[k];
                            h1r.append("," + nucA);
                        }
                        
                        vector <string> h1sample;
                        string subh1s = h1s.substr(1);
                        string subh1r = h1r.substr(1);

                        boost::split(h1sample,subh1s,boost::is_any_of(","));
                        vector <string> h1ref;
                        boost::split(h1ref,subh1r,boost::is_any_of(","));
                        
                        if (h1sample.size() == 1)
                        {
                            if (h1ref[0] == ".")     {continue;}
                            if (h1sample[0] == ".")    {miss++; miss_list.append(";" + positions[0]);continue;}
                            if (h1sample[0] != ".") {tested ++;}
                            
                            if (h1ref[0] == h1sample[0]){valid++;continue;}
                            error_list.append(";" + positions[0]);
                        }
                        
                        
                    }
                    float ratio = (valid / tested);
                    float dif = valid - tested;
                    
                    if (miss_list == "") {miss_list = ",none";}
                    miss_list = miss_list.substr(1);
                    
                    if (error_list == "") {error_list = ",none";}
                    error_list = error_list.substr(1);
                    
                    string res = sample + "\t" + cn + "\t" + chr + "\t" + alleleA + "\t" + nullallele + "\tNA\tNA\t" + to_string(int(tested)) + "\t" + to_string(int(valid)) + "\t" + error_list + "\t" + to_string(int(miss)) + "\t" + miss_list + "\t" + to_string(ratio) + "\n";
                    results_ratio[ratio].append(res);
                    results_dif[dif].append(res);
                    //break;
                }// end looop combinations
                
            }


            
            if (cn == "2")
            {
                
                for (auto comb : main_combination_list)
                {
                    vector <string> data;
                    boost::split(data,comb,boost::is_any_of("\t"));
                    string alleleA = data[0];
                    string alleleB = data[1];
                    
                    if (std::find(selected_alleles.begin(), selected_alleles.end(), alleleA) == selected_alleles.end())
                    {
                        continue;
                    }
                    if (std::find(selected_alleles.begin(), selected_alleles.end(), alleleB) == selected_alleles.end())
                    {
                        continue;
                    }
                    
                    float tested = 0;
                    float miss = 0;
                    float valid = 0;
                    string miss_list = "";
                    string error_list = "";



                    string ps_positions_str = ps_positions[sample];
                    map <string,int> ps_positions_map_nonredunt;
                    vector <string> tmp;
                    boost::split(tmp,ps_positions_str,boost::is_any_of(","));

                    for (auto item : tmp)
                    {
                        if (item == "") {continue;}
                        ps_positions_map_nonredunt[item] = 1;
                    }
                    ps_positions_str = "";
                    for (auto item : ps_positions_map_nonredunt)          
                    {
                        ps_positions_str.append("," + item.first);
                    }
                    ps_positions_str = ps_positions_str.substr(1);          


 //                   string ps_positions_str = ps_positions[sample];
 //                   ps_positions_str = ps_positions_str.substr(1);
                    vector <string> ps_positions_list;
                    boost::split(ps_positions_list,ps_positions_str,boost::is_any_of(","));
                    
                    for (auto ps : ps_positions_list)
                    //for (auto ps : ps_list)
                    {
                        string subsample = sample;
                        //string subsample = ps.first.first;
                        //if (subsample != sample) {continue;}
                        //string pos = ps.second;




                        
                        pair <string,string> t;
                        t = make_pair(subsample,ps);
                        string pos = ps_list[t];
                                    
                                    
                                   
                        vector <string> positions;
                        string substrpos = pos.substr(1);
                        boost::split(positions,substrpos,boost::is_any_of(","));
                        
                        
                        string h1s = "";
                        string h2s = "";
                        string h1r = "";
                        string h2r = "";
                        for (auto linked : positions)
                        {
                            int position = stoi(linked);
                            pair <int,string> k;
                            k = make_pair(position,sample);
                            string snp = sample_snp[k];
                            vector <string> alleles;
                            if (snp == ".") {snp = "./.";}
                            boost::replace_all(snp, "|", "/");
                            boost::split(alleles,snp,boost::is_any_of("/"));
                            h1s.append("," + alleles[0]);
                            h2s.append("," + alleles[1]);
                            
                            k = make_pair(position,alleleA);
                            string nucA = ref_data[k];
                            k = make_pair(position,alleleB);
                            string nucB = ref_data[k];
                            h1r.append("," + nucA);
                            h2r.append("," + nucB);
                        }
                        
                        vector <string> h1sample;
                        string subh1s = h1s.substr(1);
                        string subh2s = h2s.substr(1);
                        string subh1r = h1r.substr(1);
                        string subh2r = h2r.substr(1);

                        boost::split(h1sample,subh1s,boost::is_any_of(","));
                        vector <string> h2sample;
                        boost::split(h2sample,subh2s,boost::is_any_of(","));
                        vector <string> h1ref;
                        boost::split(h1ref,subh1r,boost::is_any_of(","));
                        vector <string> h2ref;
                        boost::split(h2ref,subh2r,boost::is_any_of(","));
                        
                        
                        if (h1sample.size() == 1)
                        {
                            if ((h1ref[0] == ".") && (h2ref[0] == "."))    {continue;}
                            if ((h1sample[0] == ".") && (h2sample[0] == "."))    {miss++; miss_list.append(";" + positions[0]);continue;}
                            if ((h1sample[0] == ".") || (h2sample[0] == "."))    {miss++; miss_list.append(";" + positions[0]);}
                            if ((h1sample[0] != ".") || (h2sample[0] != "."))    {tested ++;}
                            
                            if ((h1ref[0] == h1sample[0]) && (h2ref[0] == h2sample[0])){valid++;continue;}
                            if ((h1ref[0] == h2sample[0]) && (h2ref[0] == h1sample[0])){valid++;continue;}

                            if (((h1ref[0] == ".") || (h2ref[0] == ".")) && ((h1sample[0] == ".") || (h2sample[0] == "."))){valid++;continue;}
                            if (h1ref[0] == ".") {    if ((h2ref[0] == h1sample[0]) || (h2ref[0] == h2sample[0])){valid++;continue;}     }
                            if (h2ref[0] == ".") {    if ((h1ref[0] == h1sample[0]) || (h1ref[0] == h2sample[0])){valid++;continue;}     }
                            if (h1sample[0] == ".") {    if ((h2sample[0] == h1ref[0]) || (h2sample[0] == h2ref[0])){valid++;continue;}     }
                            if (h2sample[0] == ".") {    if ((h1sample[0] == h1ref[0]) || (h1sample[0] == h2ref[0])){valid++;continue;}     }
                            
                            error_list.append(";" + positions[0]);
                        }
                        
                        if (h1sample.size() > 1)
                        {
                            float validfor = 0;
                            float validrev = 0;
                            string error_list_for = "";
                            string error_list_rev = "";
                            
                            for (int a = 0; a < h1sample.size();a++)
                            {
                                if ((h1ref[a] == ".") && (h2ref[a] == "."))    {continue;}
                                
                                if ((h1sample[a] == ".") && (h2sample[a] == "."))    {miss++; miss_list.append(";" + positions[a]);continue;}
                                
                                if ((h1sample[a] == ".") || (h2sample[a] == "."))    {miss++; miss_list.append(";" + positions[a]);}
                                
                                if ((h1sample[a] != ".") || (h2sample[a] != "."))    {tested ++;}
                                
                                if ((h1ref[a] != ".") && (h2ref[a] != ".")) {
                                    if ((h1ref[a] == h1sample[a]) && (h2ref[a] == h2sample[a])){validfor++;}
                                        else {error_list_for.append(";" + positions[a]);}
                                        
                                        
                                    
                                    if ((h1ref[a] == h2sample[a]) && (h2ref[a] == h1sample[a])){validrev++;}
                                        else {error_list_rev.append(";" + positions[a]);}
                                    

                                    continue;
                                }
                                
                                if ((h1ref[a] != ".") && (h2ref[a] == ".")) {
                                    if ((h1ref[a] == h1sample[a]) || (h1ref[a] == h2sample[a])){validfor++;validrev++;continue;}
                                    else {error_list_for.append(";" + positions[a]);error_list_rev.append(";" + positions[a]); continue;}
                                }

                                if ((h1ref[a] == ".") && (h2ref[a] != ".")) {
                                    if ((h2ref[a] == h1sample[a]) || (h2ref[a] == h2sample[a])){validfor++;validrev++;continue;}
                                    else {error_list_for.append(";" + positions[a]);error_list_rev.append(";" + positions[a]);continue;}
                                }
                                
                            }
                            
                            if (validfor >= validrev)
                            {
                                valid = valid + validfor;
                                error_list.append(error_list_for);
                            }
                            if (validfor < validrev)
                            {
                                valid = valid + validrev;
                                error_list.append(error_list_rev);
                            }
                        }
        
                    }
                    float ratio = (valid / tested);
                    float dif = valid - tested;

                    if (miss_list == "") {miss_list = ",none";}
                    miss_list = miss_list.substr(1);

                    if (error_list == "") {error_list = ",none";}
                    error_list = error_list.substr(1);

                    string res = sample + "\t" + cn + "\t" + chr + "\t" + alleleA + "\t" + alleleB + "\tNA\tNA\t" + to_string(int(tested)) + "\t" + to_string(int(valid)) + "\t" + error_list + "\t" + to_string(int(miss)) + "\t" + miss_list + "\t" + to_string(ratio) + "\n";
                    results_ratio[ratio].append(res);
                    results_dif[dif].append(res);
                    //break;
                }// end looop combinations
                
            } // end 2 copies
            

            
            if (cn == "3")
            {
                
                if (call_higher_cn == 0) {
                    string res = sample + "\t" + cn + "\t" + chr + "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
                    results_ratio[1] = res;
                    results_dif[0] = res;
                }
                
                if (call_higher_cn == 1) {
                    for (auto comb : main_combination_list_3cp)
                    {
                        vector <string> data;
                        boost::split(data,comb,boost::is_any_of(","));
                        string alleleA = data[0];
                        string alleleB = data[1];
                        string alleleC = data[2];
                        
                        if (std::find(selected_alleles.begin(), selected_alleles.end(), alleleA) == selected_alleles.end())
                        {
                            continue;
                        }
                        if (std::find(selected_alleles.begin(), selected_alleles.end(), alleleB) == selected_alleles.end())
                        {
                            continue;
                        }
                        if (std::find(selected_alleles.begin(), selected_alleles.end(), alleleC) == selected_alleles.end())
                        {
                            continue;
                        }
                        
                        float tested = 0;
                        float miss = 0;
                        float valid = 0;
                        string miss_list = "";
                        string error_list = "";


                        string ps_positions_str = ps_positions[sample];
                        map <string,int> ps_positions_map_nonredunt;
                        vector <string> tmp;
                        boost::split(tmp,ps_positions_str,boost::is_any_of(","));

                        for (auto item : tmp)
                        {
                            if (item == "") {continue;}
                            ps_positions_map_nonredunt[item] = 1;
                        }
                        ps_positions_str = "";
                        for (auto item : ps_positions_map_nonredunt)          
                        {
                            ps_positions_str.append("," + item.first);
                        }
                        ps_positions_str = ps_positions_str.substr(1);          


                        
//                        string ps_positions_str = ps_positions[sample];
//                        ps_positions_str = ps_positions_str.substr(1);
                        vector <string> ps_positions_list;
                        boost::split(ps_positions_list,ps_positions_str,boost::is_any_of(","));
                    
                    for (auto ps : ps_positions_list)
                    //for (auto ps : ps_list)
                    {
                        string subsample = sample;
                        //string subsample = ps.first.first;
                        //if (subsample != sample) {continue;}
                        //string pos = ps.second;
                        
                        pair <string,string> t;
                        t = make_pair(subsample,ps);
                        string pos = ps_list[t];


                            vector <string> positions;
                            string substrpos = pos.substr(1);
                            boost::split(positions,substrpos,boost::is_any_of(","));
                            
                            string h1s = "";
                            string h2s = "";
                            string h3s = "";
                            string h1r = "";
                            string h2r = "";
                            string h3r = "";
                            for (auto linked : positions)
                            {
                                int position = stoi(linked);
                                pair <int,string> k;
                                k = make_pair(position,sample);
                                string snp = sample_snp[k];
                                vector <string> alleles;
                                if (snp == ".") {snp = "././.";}
                                boost::replace_all(snp, "|", "/");
                                boost::split(alleles,snp,boost::is_any_of("/"));
                                h1s.append("," + alleles[0]);
                                h2s.append("," + alleles[1]);
                                h3s.append("," + alleles[2]);
                                
                                k = make_pair(position,alleleA);
                                string nucA = ref_data[k];
                                k = make_pair(position,alleleB);
                                string nucB = ref_data[k];
                                k = make_pair(position,alleleC);
                                string nucC = ref_data[k];
                                h1r.append("," + nucA);
                                h2r.append("," + nucB);
                                h3r.append("," + nucC);
                            }
                            
                            vector <string> h1sample;
                            
                            string subh1s = h1s.substr(1);
                            string subh2s = h2s.substr(1);
                            string subh3s = h3s.substr(1);
                            string subh1r = h1r.substr(1);
                            string subh2r = h2r.substr(1);
                            string subh3r = h3r.substr(1);



                            boost::split(h1sample,subh1s,boost::is_any_of(","));
                            vector <string> h2sample;
                            boost::split(h2sample,subh2s,boost::is_any_of(","));
                            vector <string> h3sample;
                            boost::split(h3sample,subh3s,boost::is_any_of(","));
                            vector <string> h1ref;
                            boost::split(h1ref,subh1r,boost::is_any_of(","));
                            vector <string> h2ref;
                            boost::split(h2ref,subh2r,boost::is_any_of(","));
                            vector <string> h3ref;
                            boost::split(h3ref,subh3r,boost::is_any_of(","));
                            
                            
                            if (h1sample.size() == 1)
                            {
                                if (((h1ref[0] == ".") && (h2ref[0] == ".")) && (h3ref[0] == "."))    {continue;}
                                if (((h1sample[0] == ".") && (h2sample[0] == ".")) && (h3sample[0] == "."))    {miss++; miss_list.append(";" + positions[0]);continue;}
                                if (((h1sample[0] == ".") || (h2sample[0] == ".")) || (h3sample[0] == "."))    {miss++; miss_list.append(";" + positions[0]);}
                                if (((h1sample[0] != ".") || (h2sample[0] != ".")) || (h3sample[0] != "."))    {tested ++;}

                                
                                if (((h1ref[0] == h1sample[0]) && (h2ref[0] == h2sample[0])) && (h3ref[0] == h3sample[0])){valid++;continue;}
                                if (((h1ref[0] == h1sample[0]) && (h2ref[0] == h3sample[0])) && (h3ref[0] == h2sample[0])){valid++;continue;}
        
                                if (((h1ref[0] == h2sample[0]) && (h2ref[0] == h1sample[0])) && (h3ref[0] == h3sample[0])){valid++;continue;}
                                if (((h1ref[0] == h2sample[0]) && (h2ref[0] == h3sample[0])) && (h3ref[0] == h1sample[0])){valid++;continue;}
        
                                if (((h1ref[0] == h3sample[0]) && (h2ref[0] == h1sample[0])) && (h3ref[0] == h2sample[0])){valid++;continue;}
                                if (((h1ref[0] == h3sample[0]) && (h2ref[0] == h2sample[0])) && (h3ref[0] == h1sample[0])){valid++;continue;}
                                
                                
                                int countmissref = 0;
                                if (h1ref[0] == ".") {countmissref++;}
                                if (h2ref[0] == ".") {countmissref++;}
                                if (h3ref[0] == ".") {countmissref++;}
                                int countmisssample = 0;
                                if (h1sample[0] == ".") {countmisssample++;}
                                if (h2sample[0] == ".") {countmisssample++;}
                                if (h3sample[0] == ".") {countmisssample++;}
                                
                                if ((countmissref == 2) || (countmisssample == 2)) {valid++;continue;}
                                
                                
            
                                map <string,int> sample_alleles;
                                map <string,int> ref_alleles;
                                
                                if (h1sample[0] != ".") {sample_alleles[h1sample[0]]++;}
                                if (h2sample[0] != ".") {sample_alleles[h2sample[0]]++;}
                                if (h3sample[0] != ".") {sample_alleles[h3sample[0]]++;}
                                
                                if (h1ref[0] != ".") {ref_alleles[h1ref[0]]++;}
                                if (h2ref[0] != ".") {ref_alleles[h2ref[0]]++;}
                                if (h3ref[0] != ".") {ref_alleles[h3ref[0]]++;}
                                
                                int overlap = 0;
                                for (auto item : sample_alleles)
                                {
                                    if (ref_alleles.find(item.first) != ref_alleles.end()) {overlap++;}
                                }
                                
                            
                                if ((countmissref == 0) && (countmisssample == 0)) {overlap = 0;}
                                
                                if (overlap > 0) {valid++;continue;    }
                                
                                error_list.append(";" + positions[0]);
                            }
                            
                            
                            if (h1sample.size() > 1)
                            {
                                for (int a = 0; a < h1sample.size();a++)
                                {
                                    
                                    if (((h1ref[a] == ".") && (h2ref[a] == ".")) && (h3ref[a] == "."))    {continue;}
                                    if (((h1sample[a] == ".") && (h2sample[a] == ".")) && (h3sample[a] == "."))    {miss++; miss_list.append(";" + positions[a]);continue;}
                                    if (((h1sample[a] == ".") || (h2sample[a] == ".")) || (h3sample[a] == "."))    {miss++; miss_list.append(";" + positions[a]);}
                                    if (((h1sample[a] != ".") || (h2sample[a] != ".")) || (h3sample[a] != "."))    {tested ++;}
        
                                    if (((h1ref[a] == h1sample[a]) && (h2ref[a] == h2sample[a])) && (h3ref[a] == h3sample[a])){valid++;continue;}
                                    if (((h1ref[a] == h1sample[a]) && (h2ref[a] == h3sample[a])) && (h3ref[a] == h2sample[a])){valid++;continue;}
                                    
                                    if (((h1ref[a] == h2sample[a]) && (h2ref[a] == h1sample[a])) && (h3ref[a] == h3sample[a])){valid++;continue;}
                                    if (((h1ref[a] == h2sample[a]) && (h2ref[a] == h3sample[a])) && (h3ref[a] == h1sample[a])){valid++;continue;}
                                    
                                    if (((h1ref[a] == h3sample[a]) && (h2ref[a] == h1sample[a])) && (h3ref[a] == h2sample[a])){valid++;continue;}
                                    if (((h1ref[a] == h3sample[a]) && (h2ref[a] == h2sample[a])) && (h3ref[a] == h1sample[a])){valid++;continue;}
        
                                }
                            }
                            
                        }
                        float ratio = (valid / tested);
                        float dif = valid - tested;
                        
                        if (miss_list == "") {miss_list = ",none";}
                        miss_list = miss_list.substr(1);
                        
                        if (error_list == "") {error_list = ",none";}
                        error_list = error_list.substr(1);
                        
                        string res = sample + "\t" + cn + "\t" + chr + "\t" + alleleA + "\t" + alleleB + "\t" + alleleC + "\tNA\t" + to_string(int(tested)) + "\t" + to_string(int(valid)) + "\t" + error_list + "\t" + to_string(int(miss)) + "\t" + miss_list + "\t" + to_string(ratio) + "\n";
                        results_ratio[ratio].append(res);
                        results_dif[dif].append(res);

                        //break;
                    }// end looop combinations
               } // end 3 copies higher cn
            } // end 3 copies
            
            
                    
            mtx_genotype.lock();
            debug_message("Starting report for sample " + sample);
            mtx_genotype.unlock();
            
            ofstream report;
            string outrep = v_output_adjusted + gene + "/reports/" + sample + "." + gene + ".txt";
 //               cout << outrep << endl;
            report.open (outrep.c_str());
            
            report << "Sample\tCopy_number\tChr\tAllele_A\tAllele_B\tAllele_C\tAllele_D\tTested_genotypes\tValid_genotypes\tError_list\tMissed_genotypes\tMissed_list\tRatio" << endl;
            
            string bestcall = "";



            for (auto iter = results_dif.rbegin(); iter != results_dif.rend(); ++iter) {
                string res = iter->second;
                boost::replace_all(res, "ref.", "");
                report << res;
                if (bestcall == "") {bestcall = res;}
            }
            report.close();


            
            mtx_genotype.lock();
            debug_message("Ending report for sample " + sample);
            mtx_genotype.unlock();
            
            
            vector <string> lines;
            boost::split(lines,bestcall,boost::is_any_of("\n"));
            map <string,string> calls;
            map <string,string> calls_miss;
            for (auto item : lines)
            {
                if (item == "") {continue;}
                vector <string> data;
                boost::split(data,item,boost::is_any_of("\t"));
                string t_alleleA = data[3];
                string t_alleleB = data[4];
                string t_alleleC = data[5];
                string t_alleleD = data[6];
                string t_tested = data[7];
                string t_valid = data[8];
                string t_error_list = data[9];
                string t_missed = data[10];
                string t_missed_list = data[11];
                string t_ratio = data[12];
                
                
                if (v_exome == 1)
                {
                    if (((cn == "1") || (cn == "2")) || (cn == "3")){
                        if (t_alleleA != "NA") {
                            vector <string> fields;
                            boost::split(fields,t_alleleA,boost::is_any_of("."));
                            string last = fields[1].substr(fields[1].length() - 1);
                            string newallele = fields[0] + "*" + fields[1].substr(0,5);
                            if (((last == "N") || (last == "S")) || (last == "L")) {newallele.append(last);}
                            t_alleleA = newallele;
                        }
                    }

                    if ((cn == "2") || (cn == "3")){
                        if (t_alleleA != "NA") {
                            vector <string> fields;
                            boost::split(fields,t_alleleB,boost::is_any_of("."));
                            string last = fields[1].substr(fields[1].length() - 1);
                            string newallele = fields[0] + "*" + fields[1].substr(0,5);
                            if (((last == "N") || (last == "S")) || (last == "L")) {newallele.append(last);}
                            t_alleleB = newallele;
                        }
                    }

                    if (cn == "3"){
                        if (t_alleleC != "NA") {
                            vector <string> fields;
                            boost::split(fields,t_alleleC,boost::is_any_of("."));
                            string last = fields[1].substr(fields[1].length() - 1);
                            string newallele = fields[0] + "*" + fields[1].substr(0,5);
                            if (((last == "N") || (last == "S")) || (last == "L")) {newallele.append(last);}
                            t_alleleC = newallele;
                        }
                    }
                }
                
                
                string call = "";
                if (cn == "0") {call = t_alleleA + "+" + t_alleleB;}
                if (cn == "1") {call = t_alleleA + "+" + t_alleleB;}
                if (cn == "2") {call = t_alleleA + "+" + t_alleleB;}
                if (cn == "3") {call = t_alleleA + "+" + t_alleleB + "+" + t_alleleC;}
                if (cn == "4") {call = t_alleleA + "+" + t_alleleB + "+" + t_alleleC + "+" + t_alleleD;}
                
                calls[call] = t_ratio;
                calls_miss[call] = t_missed;
            }
                        
            string call_line = "";
            int count_calls = 0;
            float ratio = 0;
            string miss_line = "";
            for (auto item : calls)
            {
                call_line.append(";" + item.first);
                ratio = 0;
                if (cn == "0"){ratio = 1;}
                if (cn == "4"){ratio = 0;}
                if (item.second != "NA"){
                    if (cn == "1"){ratio = stof(item.second);}
                    if (cn == "2"){ratio = stof(item.second);}
                    if (cn == "3"){ratio = stof(item.second);}
                }
                count_calls++;
                miss_line.append(";" + calls_miss[item.first]);
            }
            
            string line = sample + "\t" + cn + "\t" + call_line.substr(1) + "\t" + to_string(ratio) + "\t" + miss_line.substr(1);
            
            if (count_calls > 5) {
                if (cn == "1") {line = sample + "\t" + cn + "\t" + unresolved + "\tNA\tNA";}
                if (cn == "2") {line = sample + "\t" + cn + "\t" + unresolved + "+" + unresolved + "\tNA\tNA";}
                if (cn == "3") {line = sample + "\t" + cn + "\t" + unresolved + "+" + unresolved + "+" + unresolved + "\tNA\tNA";}
            }
            

            mtx_genotype.lock();
            final_results[sample] = line;
            done++;
                float progress = (done / float(filedb.size())) * 100;
                if (v_exome == 1) {
                    v_message = " > " + gene + " : Checking SNP compatibility for CDS - " + to_string(int(progress)) + "%";
                    screen_message (screen_size, 0, v_message, 0, v_quiet);

                }
                else {
                    v_message = " > " + gene + " : Checking full SNP compatibility - " + to_string(int(progress)) + "%";
                    screen_message (screen_size, 0, v_message, 0, v_quiet);
                }
                
                
            debug_message("Ending sample " + sample);
            mtx_genotype.unlock();

            lines.clear();
            calls.clear();
            calls_miss.clear();
            results_ratio.clear();
            results_dif.clear();
            selected_alleles.clear();
            allele_dif.clear();
            
            return 1;
                    })
            );// end thread
            

 //           if (results_genotype.size() == 2000)
 //           {
 //               for (auto && result: results_genotype){result.get();} // waiting for all threads
 //               results_genotype.clear();
 //           }


            //break;
        } //end sample
        for (auto && result: results_genotype){result.get();} // waiting for all threads
        
        if (v_exome == 1) {
            v_message = " > " + gene + " : Checking SNP compatibility for CDS - done";
            screen_message (screen_size, 0, v_message, 1, v_quiet);
        }
        else {
            v_message = " > " + gene + " : Checking full SNP compatibility - done";
            screen_message (screen_size, 0, v_message, 1, v_quiet);
        }

		
        ofstream report;
        string outrep = v_output_adjusted + gene + "/calls/" + gene + ".calls.txt";
        report.open (outrep.c_str());
        report << "Sample\tCopy_number\tCalls\tRatio\tMissings" << endl;
        for (auto item : final_results)
        {
            report << item.second << endl;
        }
        report.close();
        
        
        
        // break;
        
	} // final valid_genes


    v_message = "Computation done.";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
}


