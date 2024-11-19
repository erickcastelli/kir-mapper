
//  kir-mapper
//
//  Created by Erick C. Castelli
//  2024 GeMBio.Unesp.
//  erick.castelli@unesp.br


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

#include <sys/ioctl.h> //for screen_size
#include <stdio.h> //for screen_size
#include <unistd.h> //for screen_size

#include "external.hpp"
#include "preselect.hpp"
#include "functions.hpp"
#include "map_dna.hpp"
#include "setup.hpp"
#include "ncopy.hpp"
#include "genotype.hpp"
#include "combine.hpp"
#include "haplotypes.hpp"


using namespace std;

auto clock_start = std::chrono::steady_clock::now();


// max threads
int v_concurentThreadsSupported = std::thread::hardware_concurrency() / 2;
string v_threads = std::to_string(v_concurentThreadsSupported);


// memory
unsigned long v_memory = getTotalSystemMemory() / long(1024) / long(1024) / long(1024);


//Program identification
string Program_name = "kir-mapper";
string Program_company = "GeMBio.Unesp";
string Program_version = "1.0";
string Program_author = "Erick C. Castelli";
string Program_date = "November 20th 2024";
string Program_website = "https://github.com/erickcastelli/kir-mapper";


map <string,int> samtools_versions;
map <string,int> bwa_versions;
map <string,int> whats_versions;
map <string,int> freebayes_versions;
map <string,int> bcf_versions;
map <string,int> R_versions;
map <string,int> star_versions;
map <string,int> minimap_versions;


// Internal variables
string ostype = "";
int screen_size = 80;
string v_message = "";
string v_system_out = "";
vector<string> warnings;
string configfile = "";
vector<string> v_system_out_list;
string v_command;
vector <string> selected_reads;
int motif_size = 30;
int v_buffer = 1000000;
string v_genes = ""; // list of genes avaliable
vector<string> v_gene_list;
int v_mm_max = 20; //-n BWA ALN
int v_mm_open = 1; //-o BWA ALN
int v_size = 50;
double v_tolerance = 0.05; // tolerance
double v_mtrim_error = 0.08f; // error limit (mtrim)
//int keep_sam = 0;
int minselect = 30;
string v_bed = "";
string v_map_type = "";
int downsampling = 30;
//int v_rnaseq = 0;
int v_lowmem = 0;
int v_callint = 0;
//int v_forceindex = 0;
int v_exome = 0;
int v_exons = 0;
int v_cds = 1;
int v_rna = 0;
int v_full = 0;
int v_nanopore = 0;
string v_thresholds;
string v_target;
string v_list;
string v_reference_name = "";
string v_sample_list = "";

map <string,int> position_hg38;
map <string,string> chr_hg38;
map <string,string> chr_size_hg38;
map <string,int> gene_opt_start;
map <string,int> gene_opt_end;
map <string,string> gene_type;
map <string,int> ref_size;

map <pair<string,string>,int> reads_mapped_to_allele_A;
map <pair<string,string>,int> reads_mapped_to_allele_B;
map <pair<string,string>,int> reads_possibly_mapped;
map <string,int> reads_mapped_nm;
map <string,string> reads_mapped_gene;

double read_mean_size;
map <string,double> read_count_gene;

map <pair<string,string>, int> sequence_list_r1;
map <pair<string,string>, int> sequence_list_r2;
map <string, int> sequence_size_r1;
map <string, int> sequence_size_r2;

map <pair<string,string>, int> sequence_address;

int v_skiptyping = 0;
int v_typeall = 0;
int v_quiet = 0;
int v_useconfig = 0;
int v_skip_unmapped = 0;
int v_nopolyphasing = 0;
int v_recall = 0;
int v_haploforce = 0;

//Main parameters
string v_output = "";
string v_db = "";
string v_dbprepare = "";
string v_r0 = "";
string v_r1 = "";
string v_r2 = "";
string v_bam = "";
string v_mapout = "";
string v_sample = "";
string homedir = "";
int usechr = 0;
string v_tag = "";
int v_force_two_copies = 0;
int v_update_call = 0;

// For developers
int v_debug = 0; // defines the debug mode
int v_verbose = 0; // defines the verbose mode
int v_bypass_version_check = 0; //bypass bwa and samtools version
int v_replicates = 20;

string v_bwa = "";
string v_samtools = "";
string v_whats = "";
string v_freebayes = "";
string v_bcftools = "";
string v_picard = "";
string v_star = "";
string v_shapeit = "";
string v_minimap = "";

map <string,string> capture_bias_region;
map <string,float> exome_depth;
map <pair<string,string>,string> genes_avaliable_in_sample;

void main_help(void)
{
    screen_message (screen_size, 0, "", 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);

    v_message = "Program:  " + Program_name;
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "Version:  " + Program_version + ", " + Program_date;
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    screen_message (screen_size, 0, "Contact:  Erick C. Castelli <erick.castelli@unesp.br>", 1, v_quiet);

    screen_message (screen_size, 0, "", 1, v_quiet);
    screen_message (screen_size, 0, "Usage:    kir-mapper <command> [options]", 1, v_quiet);

    screen_message (screen_size, 0, "", 1, v_quiet);
    
    screen_message (screen_size, 0, "Commands:  setup         configure kir-mapper", 1, v_quiet);
    screen_message (screen_size, 0, "           map           map/align sequences (WGS, WES, Amplicons)", 1, v_quiet);
    screen_message (screen_size, 0, "           ncopy         detect copy numbers", 1, v_quiet);
    screen_message (screen_size, 0, "           genotype      call SNPs and alleles", 1, v_quiet);
    screen_message (screen_size, 0, "           haplotype     estimate haplotypes and call alleles", 1, v_quiet);
    screen_message (screen_size, 0, "           group         combine multiple kir-mapper map and ncopy runs", 1, v_quiet);
    screen_message (screen_size, 0, "           join          join variants into a single VCF for plink", 1, v_quiet);
    screen_message (screen_size, 0, "           select        preselect KIR-like sequences", 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);
    return;
}




void load_config (void)
{
    if (v_useconfig == 0) {
        ifstream config(configfile.c_str());
        if (config)
        {
            for( std::string line; getline( config, line ); )
            {
                line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                vector<string> v_set;
                boost::split(v_set,line,boost::is_any_of("="));
                if (v_set[0] == "bwa") {v_bwa = v_set[1];}
                if (v_set[0] == "samtools") {v_samtools = v_set[1];}
                if (v_set[0] == "db") {v_db = v_set[1] + "/";}
                if (v_set[0] == "whatshap") {v_whats = v_set[1];}
                if (v_set[0] == "freebayes") {v_freebayes = v_set[1];}
                if (v_set[0] == "bcftools") {v_bcftools = v_set[1];}
                if (v_set[0] == "minimap") {v_minimap = v_set[1];}
                if (v_set[0] == "picard") {
					v_picard = v_set[1]; 
					if (v_picard == "") {v_picard = "DISABLED";}
				}
                if (v_set[0] == "star") {v_star = v_set[1];}
                if (v_set[0] == "shapeit4") {v_shapeit = v_set[1];}

            }
        }
        config.close();
    }
    return;
}




int main(int argc, const char * argv[]) {

    
    
    #if defined (__linux__)
        ostype = "linux";
    #endif
    #if defined (__APPLE__)
        ostype = "mac";
    #endif
    
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    //screen_size = w.ws_col;
    
    
    if (stoi(v_threads) < 1) {v_threads = "2";}
    if (stoi(v_threads) > 20) {v_threads = "20";}
    
    
    // tested program versions
    samtools_versions["1.18"] = 1;
    samtools_versions["1.19.2"] = 1;

    bwa_versions["0.7.17"] = 1;
    bwa_versions["0.7.16"] = 1;
    
    whats_versions["2.0"] = 1;
    whats_versions["2.1"] = 1;
    whats_versions["2.2"] = 1;
    whats_versions["2.3"] = 1;

    freebayes_versions["v1.3.6"] = 1;

    bcf_versions["1.13"] = 1;
    bcf_versions["1.19"] = 1;
	
	minimap_versions["2.24"] = 1;

    R_versions["3.6.3"] = 1;
    R_versions["4.1.1"] = 1;
    R_versions["4.2.2"] = 1;
    R_versions["4.3.1"] = 1;

    star_versions["2.7.11a"] = 1;
    star_versions["2.7.11b"] = 1;
    star_versions["2.7.10b"] = 1;

    homedir = getpwuid(getuid())->pw_dir;
    configfile = homedir + "/.kir-mapper";
    
    
    if (fileExists(configfile))
    {
        load_config();
    }
    else
    {
        cout << endl << endl;
        cout << endl << endl;
        cout << "Apparently, this is the first time you have run kir-mapper." << endl;
        cout << "Starting configuration ..." << endl;
        main_setup();
        return 0;
    }
    
    
    if (! fileExists(v_bwa)) {cout << endl << "You should run 'kir-mapper setup'. BWA not detected." << endl << endl; main_setup(); return 0; }

    if (! fileExists(v_samtools)) {cout << endl << "You should run 'kir-mapper setup'. Samtools not detected." << endl << endl; main_setup(); return 0;}

    if (! fileExists(v_bcftools)) {cout << endl << "You should run 'kir-mapper setup'. Bcftools not detected." << endl << endl; main_setup(); return 0;}
    
    if (! fileExists(v_db)) {cout << endl << "You should run 'kir-mapper setup'. Database not detected." << endl << endl; main_setup(); return 0;}

    if (! fileExists(v_star)) {cout << endl << "You should run 'kir-mapper setup'. STAR not detected." << endl << endl; main_setup(); return 0;}

    if (! fileExists(v_minimap)) {cout << endl << "You should run 'kir-mapper setup'. Minimap2 not detected." << endl << endl; main_setup(); return 0;}

    if (v_whats != "DISABLED") {
        if (! fileExists(v_whats)) {cout << endl << "You should run 'kir-mapper setup'. Whatshap not detected." << endl << endl; main_setup(); return 0;}
    }
  
    if (v_picard != "DISABLED") {
        if (! fileExists(v_picard)) {cout << endl << "You should run 'kir-mapper setup'. Picard tools not detected." << endl << endl; main_setup(); return 0;}
    }
    
    if (! fileExists(v_freebayes)) {cout << endl << "You should run 'kir-mapper setup'. Freebayes not detected." << endl << endl; main_setup(); return 0;}

    
    int command_ok = 1;
    
    if (argc > 2)
    {
       int a;
       for (a = 2; a < argc; a++)
       {
           string str = argv[a];
           
           if (str == "--quiet")
           {
               v_quiet= 1;
               continue;
           }

           if (str == "--force")
           {
               v_haploforce = 1;
               continue;
           }
		   
		   if ((str == "--update_calls") || (str == "--update_call"))
           {
               v_update_call = 1;
               continue;
           }


           else if (str == "--exome")
           {
               v_exome = 1;
               continue;
           }
           
           else if (str == "--cds")
           {
               v_cds = 1;
               continue;
           }

           else if (str == "--wgs")
           {
               v_exome = 0;
               continue;
           }

           else if (str == "--recall")
           {
               v_recall = 1;
               continue;
           }

           else if (str == "--full")
           {
               v_exome = 0;
               v_full = 1;
               continue;
           }
           
		   
		    else if (str == "--nanopore")
           {
               v_exome = 0;
               v_full = 1;
			   v_nanopore = 1;
			   v_skiptyping = 1;
			   v_tolerance = 0.1;
               continue;
           }
           
           else if (str == "--exons")
           {
               v_exons = 1;
               v_cds = 0;
               continue;
           }

           else if (str == "--telomeric")
           {
               v_target= "KIR2DL4,KIR3DL1,KIR3DS1,KIR2DL5AB,KIR2DS3,KIR2DS5,KIR2DS1,KIR2DS4,KIR3DL2";
               continue;
           }
        
           else if (str == "--centromeric")
           {
               v_target= "KIR3DL3,KIR2DS2,KIR2DL2,KIR2DL3,KIR2DL5AB,KIR2DS3,KIR2DP1,KIR2DL1,KIR3DP1";
               continue;
           }
 

           else if (str == "--nopolyphase")
           {
               v_nopolyphasing = 1;
               continue;
           }
           
           else if (str == "--noconfig")
           {
               v_useconfig = 1;
               continue;
           }

           else if (str == "--low-mem")
           {
               v_lowmem = 1;
               continue;
           }

            else if (str == "--debug")
           {
               v_debug = 1;
               continue;
           }

           
           else if (str == "--verbose")
           {
               v_verbose = 1;
               continue;
           }
           
           else if (str == "--bypass-version-check")
           {
               v_bypass_version_check = 1;
               continue;
           }
           

           else if (str == "--skip-adjust")
           {
               v_skiptyping = 1;
               continue;
           }
           
           else if (str == "--rna")
           {
               v_rna = 1;
               v_exome = 0;
               continue;
           }
           
           else if (str == "--skip-unmapped")
           {
               v_skip_unmapped = 1;
               continue;
           }

           else if (str == "--two_copies")
           {
               v_force_two_copies = 1;
               continue;
           }



           
           else if (str =="-r0")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_r0 = argv[b];
               if (! fileExists(v_r0))
               {
                   v_message = "Invalid r0: " + v_r0;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }

           
           else if (str == "-r1")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_r1 = argv[b];
               if (! fileExists(v_r1))
               {
                   v_message = "Invalid r1: " + v_r1;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }


           else if (str == "-r2")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_r2 = argv[b];
               if (! fileExists(v_r2))
               {
                   v_message = "Invalid r2: " + v_r2;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }

           else if (str == "-bam")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_bam = argv[b];
               if (! fileExists(v_bam))
               {
                   v_message = "Invalid bam: " + v_bam;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }
           
           else if (str == "-map")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_mapout = argv[b];
               if (! fileExists(v_mapout))
               {
                   v_message = "Invalid map folder: " + v_mapout;
                   v_mapout = "";
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }


           else if (str == "-bed")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_bed = argv[b];
               if (! fileExists(v_bed))
               {
                   v_message = "Invalid bed: " + v_bed;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }
           
           else if (str == "-threads")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_threads  = argv[b];
               if (std::stoi(v_threads) < 1) {v_threads = "2";}
			   if (std::stoi(v_threads) > 20) {v_threads = "20";}
               continue;
           }

           else if (str == "-reference")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_reference_name  = argv[b];
               int subtest = 0;
               if (((v_reference_name == "KIR3DL3") || (v_reference_name == "5UPKIR")) || ((v_reference_name == "HLA-G") || (v_reference_name == "HLA-E")))
               {subtest = 1;}
               if (subtest == 0)
               {
                   v_message = "Invalid reference: " + v_reference_name;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }
           
           else if (str == "-db")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_db  = argv[b];
               char last_ch = v_db.back();
               char ch = '/';
               if (last_ch != ch)
               {
                   v_db = v_db + "/";
               }
               v_dbprepare = v_db;
               
               if (! fileExists(v_db))
               {
                   v_message = "Invalid db: " + v_db;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               
               continue;
           }
         
           else if (str == "-buffer")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_buffer  = stoi(argv[b]);
               if (v_buffer < 10000) {v_buffer = 10000;}
               continue;
           }

           else if (str == "-downsample")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               downsampling  = stoi(argv[b]);
               downsampling = stoi(str.substr(11));
               if (downsampling < 5) {downsampling = 5;}
               continue;
           }
        
           else if (str == "-output")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_output  = argv[b];
               char last_ch = v_output.back();
               char ch = '/';
               if (last_ch != ch)
               {
                   v_output = v_output + "/";
               }
               continue;
           }
           
           else if (str == "-sample")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_sample  = argv[b];
               if (v_sample != "") {v_sample = v_sample + "_";}
               continue;
           }

           else if (str == "-list")
           {
              int b = a + 1;
              if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
              v_list  = argv[b];
               continue;
           }
           
           else if (str == "-samples")
           {
              int b = a + 1;
              if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
              v_sample_list  = argv[b];
               if (! fileExists(v_sample_list))
               {
                   v_message = "Invalid file:" + v_sample_list;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }


           
           else if (str == "-samtools")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_samtools  = argv[b];
               if (! fileExists(v_samtools))
               {
                   v_message = "Invalid samtools: " + v_samtools;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }

           else if (str == "-whatshap")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_whats  = argv[b];
               if (! fileExists(v_whats))
               {
                   v_message = "Invalid WhatsHap: " + v_whats;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }
           
           
           else if (str == "-bwa")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_bwa  = argv[b];
               if (! fileExists(v_bwa))
               {
                   v_message = "Invalid bwa: " + v_bwa;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }
		   
		   else if (str == "-minimap")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_minimap  = argv[b];
               if (! fileExists(v_bwa))
               {
                   v_message = "Invalid minimap: " + v_minimap;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }


 
           
           else if (str == "-size")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_size  = stoi(argv[b]);
               if (v_size < 30) {v_size = 30;}
               continue;
           }

           else if (str == "-replicates")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_replicates  = stoi(argv[b]);
               continue;
           }

           else if (str == "-threshold")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_thresholds = argv[b];
               continue;
           }

           else if (str == "-tag")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_tag  = argv[b];
               continue;
           }
           
           else if (str == "-target")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_target  = argv[b];
               continue;
           }
           
           else if (str == "-tolerance")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_tolerance  = stod(argv[b]);
               continue;
           }

           else if (str == "-error")
           {
               int b = a + 1;
               if (b > (argc-1)) {continue;}
               string test = argv[b];
               if (test.substr(0,1) == "-") {continue;}
               v_mtrim_error = stod(argv[b]);
               continue;
           }
           
           else
           {
                if (str.substr(0,1) == "-") {
                    v_message = "Unknown option: " + str;
                    warnings.push_back(v_message);
                    command_ok = 0;
                }
           }
       }
    }


   
    
    if (argc == 1)
    {
        main_help();
 
        if (warnings.size() > 0)
        {
            v_message = "Warning:  " + warnings[0];
            screen_message (screen_size, 0, v_message, 1, v_quiet);
        }
        if (warnings.size() > 1)
        {
            for(int a = 1; a < warnings.size(); a++)
                v_message = "          " + warnings[a];
                screen_message (screen_size, 0, v_message, 1, v_quiet);
        }
        screen_message (screen_size, 0, "", 1, v_quiet);
  
        return (0);
    }
    
    
    
    if (command_ok == 0) {v_r1 = ""; v_r0 = ""; v_bam = ""; v_output = ""; v_sample="";}
    command_ok = 0;
    
    
 
    if (strcmp(argv[1],"map") == 0)
    {
        main_dna_map();
        command_ok = 1;
    }
    
    if (strcmp(argv[1],"select") == 0)
    {
        main_preselect();
        command_ok = 1;
    }
    
    if (strcmp(argv[1],"ncopy") == 0)
    {
        main_ncopy();
        command_ok = 1;
    }
    
    if (strcmp(argv[1],"genotype") == 0)
    {
        main_genotype();
        command_ok = 1;
    }
    
    if (strcmp(argv[1],"genotypes") == 0)
    {
        main_genotype();
        command_ok = 1;
    }

    if (strcmp(argv[1],"setup") == 0)
    {
        main_setup();
        command_ok = 1;
    }

    if (strcmp(argv[1],"group") == 0)
    {
        main_combine();
        command_ok = 1;
    }

    if (strcmp(argv[1],"combine") == 0)
    {
        main_combine();
        command_ok = 1;
    }
    
    if (strcmp(argv[1],"haplotypes") == 0)
    {
        main_haplotypes();
        command_ok = 1;
    }
    if (strcmp(argv[1],"haplotype") == 0)
    {
        main_haplotypes();
        command_ok = 1;
    }
	

    if (strcmp(argv[1],"join") == 0)
    {
        main_join();
        command_ok = 1;
    }
  
    if (command_ok == 0)
    {
        string str = argv[1];
        v_message = "Unknown command: " + str;
        warnings.push_back(v_message);
        main_help();
        if (warnings.size() > 0)
        {
            screen_message (screen_size, 0, "", 1, v_quiet);
            v_message = "Warning:  " + warnings[0];
            screen_message (screen_size, 0, v_message, 1, v_quiet);
        }
        if (warnings.size() > 1)
        {
            for (int a = 1; a < warnings.size(); a++) {
                v_message = "          " + warnings[a];
                screen_message (screen_size, 0, v_message, 1, v_quiet);
            }
        }
        screen_message (screen_size, 0, "", 1, v_quiet);
        return (0);
    }

    
    
    
    auto clock_end = std::chrono::steady_clock::now();
    auto diff = clock_end - clock_start;
    v_message = "Elapsed time: " + to_string(((std::chrono::duration <double, std::milli> (diff).count())/1000)) + " s";
    warnings.push_back(v_message);
    
    
    if (warnings.size() > 0)
    {
        screen_message (screen_size, 0, "", 1, v_quiet);
        v_message = "Warning:  " + warnings[0];
        screen_message (screen_size, 0, v_message, 1, v_quiet);
    }
    
    if (warnings.size() > 1)
    {
        for (int a = 1; a < warnings.size(); a++) {
            v_message = "          " + warnings[a];
            screen_message (screen_size, 0, v_message, 1, v_quiet);
        }
    }
    screen_message (screen_size, 0, "", 1, v_quiet);

    return 0;
}
