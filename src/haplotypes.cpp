//  kir-mapper
//
//  Created by Erick C. Castelli
//  2022 GeMBio.Unesp.
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

mutex haplo_mtx;

void main_haplotypes() {
	
    
    int v_check = 0;
    if (v_output != "") {v_check = 1;}
    if (! fileExists(v_shapeit)) {
        v_message = "Shapeit4 must be available for kir-mapper haplotypes.";
        warnings.push_back (v_message);
        v_message = "Please run setup again!";
        warnings.push_back (v_message);
        return;
    }
 
    if (((v_db == "") || (v_check == 0)))
    {
        screen_message (screen_size, 0, "", 1, v_quiet);
        v_message = "Program:   " + Program_name + "::haplotypes";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        v_message = "Version:   " + Program_version + ", " + Program_date;
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        v_message = "Usage:     kir-mapper haplotype -output kir-mapper_output_folder <options>";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Mandatory options:", 1, v_quiet);
        screen_message (screen_size, 2, "-output          output folder (same as map and ncopy)", 1, v_quiet);
        screen_message (screen_size, 2, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Other options:", 1, v_quiet);
        screen_message (screen_size, 2, "-db              path to the kir-mapper database", 1, v_quiet);
        screen_message (screen_size, 2, "-threads         number of threads [" + v_threads + "]", 1, v_quiet);
        screen_message (screen_size, 2, "-replicates      number of replicates [" + to_string(v_replicates) + "]", 1, v_quiet);
        screen_message (screen_size, 2, "-tag             tag to differentiate multiple haplotype runs", 1, v_quiet);
        screen_message (screen_size, 2, "-target          list of target genes (e.g., KIR2DL4,KIR3DL3)", 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        screen_message (screen_size, 2, "--cds           include only the CDS (no 3'UTR) - default", 1, v_quiet);
        screen_message (screen_size, 2, "--exons         include exons and 3'UTR - need genotype --full", 1, v_quiet);
//        screen_message (screen_size, 2, "--full          full genotype, include introns and UTRs", 1, v_quiet);
//        screen_message (screen_size, 2, "--recall        only recall alleles, skip shapeit4", 1, v_quiet);
        screen_message (screen_size, 2, "--telomeric     limit target to telomeric genes", 1, v_quiet);
        screen_message (screen_size, 2, "--centromeric   limit target to centromeric genes", 1, v_quiet); 
//        screen_message (screen_size, 2, "--short         do not include 2DP1 and 3DP1", 1, v_quiet); 


        screen_message (screen_size, 2, "--quiet         quiet mode", 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        return;
    }


    if (v_full == 1)
    {
        v_exons = 0;
        v_cds = 0;
        v_full = 1;
    }

    if (v_exons == 1)
    {
        v_exons = 1;
        v_cds = 0;
        v_full = 0;
    }




    map <string,int> fake_gene_positions;
    fake_gene_positions["KIR3DL3"] = 50000;
    fake_gene_positions["KIR2DS2"] = 100000;
    fake_gene_positions["KIR2DL2"] = 150000;
    fake_gene_positions["KIR2DL3"] = 200000;
    fake_gene_positions["KIR2DS3"] = 250000;
    fake_gene_positions["KIR2DP1"] = 300000;
    fake_gene_positions["KIR2DL1"] = 350000;
    fake_gene_positions["KIR3DP1"] = 400000;

    fake_gene_positions["KIR2DL4"] = 600000;
    fake_gene_positions["KIR3DL1"] = 650000;
    fake_gene_positions["KIR3DS1"] = 700000;
    fake_gene_positions["KIR2DL5AB"] = 750000;
    fake_gene_positions["KIR2DS5"] = 800000;
    fake_gene_positions["KIR2DS1"] = 850000;
    fake_gene_positions["KIR2DS4"] = 900000;
    fake_gene_positions["KIR3DL2"] = 950000;

    vector <string> list_of_genes_expected_positions;
    list_of_genes_expected_positions.push_back("KIR3DL3");
    list_of_genes_expected_positions.push_back("KIR2DS2");
    list_of_genes_expected_positions.push_back("KIR2DL2");
    list_of_genes_expected_positions.push_back("KIR2DL3");
    list_of_genes_expected_positions.push_back("KIR2DS3");
    list_of_genes_expected_positions.push_back("KIR2DP1");
    list_of_genes_expected_positions.push_back("KIR2DL1");
    list_of_genes_expected_positions.push_back("KIR3DP1");

    list_of_genes_expected_positions.push_back("KIR2DL4");
    list_of_genes_expected_positions.push_back("KIR3DL1");
    list_of_genes_expected_positions.push_back("KIR3DS1");
    list_of_genes_expected_positions.push_back("KIR2DL5AB");
    list_of_genes_expected_positions.push_back("KIR2DS5");
    list_of_genes_expected_positions.push_back("KIR2DS1");
    list_of_genes_expected_positions.push_back("KIR2DS4");
    list_of_genes_expected_positions.push_back("KIR3DL2");


    if (! fileExists(v_output))
    {
        v_message = "You must indicate a valid kir-mapper output folder";
        warnings.push_back (v_message);
        return;
    }




    screen_message (screen_size, 0, "", 1, v_quiet);
    v_message = "Program:   " + Program_name + "::haplotype";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "Version:   " + Program_version + ", " + Program_date;
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "Cores:     using " + v_threads + " thread(s)";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    if (v_exons == 1) {
        v_message = "Mode:      only exons and 3'UTR";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
    }

    if ((v_exons == 0) && (v_cds == 1)) {
        v_message = "Mode:      only CDS, no 3'UTR";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
    }

    if (v_full== 1) {
        v_message = "Mode:      full genotype";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
    }




    string master_output = "";
    if (v_tag == ""){master_output = v_output + "/haplotypes";}
    if (v_tag != ""){master_output = v_output + "/haplotypes_" + v_tag;}
    v_command = "mkdir " + master_output;
    string out = GetStdoutFromCommand(v_command);
    if (! fileExists(master_output)) {cout << "Something went wrong!" << endl;}


    if (v_full == 1){master_output.append("/full/");}
    if (v_cds == 1){master_output.append("/cds/");}
    if (v_exons == 1){master_output.append("/exons/");}


    v_command = "mkdir " + master_output;
    out = GetStdoutFromCommand(v_command);
    if (! fileExists(master_output)) {cout << "Something went wrong!" << endl;}




    map <pair<int,string>,string> vcf_data;
    map <pair<int,string>,string> ref_data;

    map <string,string> snp_info;
    map <int,string> snp_correct_position;
    map <int,string> fakepos_gene;
    map <pair<string,string>,string> ps_fake_pos;
    map <string,string> chr_reference;




    string v_source_genotype = "";
    if (v_cds == 1) {v_source_genotype = v_output + "/genotype/cds/";}
    if ((v_full == 1) || (v_exons == 1)) {v_source_genotype = v_output + "/genotype/full/";}
    if (! fileExists(v_source_genotype)) {screen_message(screen_size, 0, "You need a kir-mapper output folder with genotypes (function genotype)",1,0); return;}




    string v_cn = v_output + "/ncopy/copy_numbers.txt";
    if (! fileExists(v_cn)) {screen_message(screen_size, 0, "You need a kir-mapper output folder with copy numbers (function ncopy)",1,0); return;}
    screen_message(screen_size, 1, " > Loading copy numbers ...",1,v_debug); 
    ifstream cn (v_cn.c_str());
    string line = "";
    getline( cn, line );
    map <pair<string,string>,string> copy_numbers;
    for( std::string line; getline( cn, line ); )
    {
        if (line == "") {continue;}
        if (line.substr(0,2) == "##") {continue;}
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));
        pair <string,string> k;
        k = make_pair(data[1],data[0]);
        copy_numbers[k]= data[2];
    }
    cn.close();




    map <string,int> target_genes;
    
    if (v_target != "") {
        vector <string> target_split;
        boost::split(target_split,v_target,boost::is_any_of(","));
        for (auto item : target_split)
        {
            if (item == "") {continue;}
            if (fake_gene_positions.find(item) == fake_gene_positions.end()) {continue;}
            target_genes[item] = 1;
        }
        if (target_genes.size() == 0){cout << "Invalid list of targets!" << endl; return;}
    }

    






    map <string,int> valid_genes;
    map <string,int> sample_list;
    map <string,int> gene_end_position;
    map <string,int> ref_list;

    for (auto gene : fake_gene_positions)
    {
        

        if (v_target != "")
        {
            if (target_genes.find(gene.first) == target_genes.end()) {continue;}
        }
        
        string test_file = v_source_genotype + "/" + gene.first + "/vcf/" + gene.first + ".combined.trim.treated.norm.phased.vcf";
        if (! fileExists(test_file)) {test_file = v_source_genotype + "/" + gene.first + "/vcf/" + gene.first + ".combined.trim.treated.norm.vcf";}
        
        if (! fileExists(test_file)) {
            screen_message(screen_size, 1, " > Skipping " + gene.first + ": no VCF file.",1,0); 
            continue;
        }


        valid_genes[gene.first] = 1;
        screen_message(screen_size, 1, " > Loading " + gene.first + " ...",1,v_debug); 


        ifstream vcf (test_file.c_str());
        vector <string> samples;
        int first = 0;
        for( std::string line; getline( vcf, line ); )
        {
            if (line == "") {continue;}
            if (line.substr(0,2) == "##") {continue;}
            if (line.substr(0,2) == "#C") {
                boost::split(samples,line,boost::is_any_of("\t"));
                continue;}
            vector <string> data;
            boost::split(data,line,boost::is_any_of("\t"));

            chr_reference[gene.first] = data[0];

            if (first == 0) {first = stoi(data[1]);}

            int fake_pos = (stoi(data[1]) - first) + fake_gene_positions[gene.first] + 200;
            snp_correct_position[fake_pos] = data[0] + "\t" + data[1];
            snp_info[data[0] + "\t" + data[1]] = data[2] + "\t" + data[3] + "\t" + data[4];
            fakepos_gene[fake_pos] = gene.first;
            gene_end_position[gene.first] = fake_pos;

            pair <string,string> k;
            k = make_pair(gene.first, data[1]);
            ps_fake_pos[k] = to_string(fake_pos);


            for (int a = 9; a < data.size(); a++)
            {
                
                pair <int,string> k;
                k = make_pair(fake_pos,samples[a]);
                if (samples[a].find("ref.") != std::string::npos) {
                    ref_data[k] = data[a]; 
                    ref_list[samples[a]]++;
                    continue;
                }
                sample_list[samples[a]]++;
                vcf_data[k] = data[a];
            }
        }
        vcf.close();
    }


    if (sample_list.size() < 100)
    {
        cout << "The minimum sample size required to this step is 100. Quitting..." << endl;
        return;    
    }



    for (auto gene : fake_gene_positions)
    {
        int pos = gene.second;
        
        for (auto sample : sample_list)
        {
            
            pair <string,string> k;
            k = make_pair(gene.first,sample.first);
            string cn = copy_numbers[k];

            pair <int,string> j;
            j = make_pair(pos,sample.first);
            
            if (cn == "NA"){vcf_data[j] = "./.:.";}
            if (cn == "0"){vcf_data[j] = "1/1:.";}
            if (cn == "1"){vcf_data[j] = "1|0:" + to_string(pos);}
            if (cn == "2"){vcf_data[j] = "0/0:.";}
            if (cn == "3"){vcf_data[j] = "./.:.";}
            if (cn == "4"){vcf_data[j] = "./.:.";}

            snp_correct_position[pos] = "fakechr\t" + to_string(pos);
            snp_info["fakechr\t" + to_string(pos)] = gene.first + "\t<PRESENCE>\t<ABSENSE>";
            fakepos_gene[pos] = gene.first;
        }
   
    }












    screen_message(screen_size, 1, " > Writing fake diploid VCF ...",1,v_debug); 

     
    ofstream fakevcf;
    string file = master_output + "/fake_diploid_for_shapeit4.vcf";;
    fakevcf.open (file.c_str());
    fakevcf << "##fileformat=VCFv4.2" << endl;
    fakevcf << "##contig=<ID=fakechr,length=3000000>" << endl;
    fakevcf << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
    fakevcf << "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">" << endl;
    

    string lineout = "";
    lineout.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (auto sample : sample_list)
    {
        lineout.append("\t" + sample.first);
    }
    fakevcf << lineout << endl;

    string valid_gene = "";
    int first_pos = 0;
    for (auto pos : snp_correct_position)
    {
        
        lineout = "";
        string realpos = snp_correct_position[pos.first];
        lineout.append("fakechr\t" + to_string(pos.first) + "\t");
        lineout.append(snp_info[realpos]);
        lineout.append("\t.\t.\t.\tGT:PS");

        string gene = fakepos_gene[pos.first];
        if (gene != valid_gene)
        {
            valid_gene = gene;
            first_pos = pos.first;
        }

        float missing = 0;
        for (auto sample : sample_list)
        {
            pair <string,string> k;
            k = make_pair(gene,sample.first);
            string cn = copy_numbers[k];

            pair <int,string> j;
            j = make_pair(pos.first,sample.first);
            string info = vcf_data[j];

            
            vector <string> data;
            boost::split(data,info,boost::is_any_of(":"));

            string outgen = "";
            
            if (first_pos == pos.first)
            {
                outgen = data[0] + ":" + data[1];
            }



            if (first_pos != pos.first) 
            {
                if (cn == "0") {outgen = "0/0:.";}

                if (cn == "1") {
                    outgen.append("0|");
                    outgen.append(data[0]);
                    outgen.append(":" + to_string(first_pos));

                    if (data[0].find(".") != std::string::npos) {
                        missing++;
                    }

                }
                if (cn == "2") {
                    outgen.append(data[0]);
                    outgen.append(":");
    
                    string ps = ".";
                    if (data[3] != ".") {
                        pair <string,string> k;
                        k = make_pair(gene, data[3]);
                        ps = ps_fake_pos[k];
                    }
                    outgen.append(ps);

                    if (data[0].find(".") != std::string::npos) {
                        missing++;
                    }
                }

                if (cn == "3") {
                    
                    vector <string> sub;
                    boost::replace_all(data[0], "|", "/");
                    boost::split(sub,data[0],boost::is_any_of("/"));

                    map <string,int> alleles;
                    for (auto item : sub) {alleles[item]++;}
                    sub.clear();
                    for (auto item: alleles) {sub.push_back(item.first);}
                    if (sub.size() == 1){
                        outgen.append(sub[0] + "/" + sub[0] + ":.");
                    }
                    if (sub.size() == 2){
                        outgen.append(sub[0] + "/" + sub[1] + ":.");
                    }

                    if (sub.size() > 2) {outgen.append("./.:.");missing++;}                    
                    if (sub.size() == 0) {outgen.append("./.:.");missing++;}                    
                    
                    
                }
                if (cn == "4") {
                    outgen.append("./.:."); missing++;
                }
                if (cn == "NA") {
                    outgen.append("./.:."); missing++;
                }
            }
            lineout.append("\t" + outgen);
        }

        fakevcf << lineout << endl;
           
    }
    fakevcf.close();




    string outnorm = master_output + "/fake_diploid_for_shapeit4.bi.vcf";
    v_command = v_bcftools + " norm -m-any -o " + outnorm + " " + file;
    out = GetStdoutFromCommand(v_command);


    v_command = "bgzip -f " + outnorm;
    out = GetStdoutFromCommand(v_command);

    v_command = "tabix -p vcf " + outnorm + ".gz";
    out = GetStdoutFromCommand(v_command);












    screen_message(screen_size, 1, " > Phasing with shapeit4 ...",1,v_debug); 

    if (v_recall == 1){screen_message(screen_size, 1, " > Phasing with shapeit4 ... skipping ...",1,v_debug); }



    string shapeit_dir = master_output + "/shapeit4";
    v_command = "mkdir " + shapeit_dir;
    out = GetStdoutFromCommand(v_command);


    if (v_recall == 0) {

        vector <string> cmd_list;
        for (int a = 1; a <= v_replicates; a++)
        {
            string outshapeit = shapeit_dir + "/shapeit_run" + to_string(a) + ".vcf";
            float seed = a + 1;
            v_command = v_shapeit + " --input " + outnorm + ".gz --region fakechr --output " + outshapeit + " --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,1b,1p,10m --sequencing --use-PS 0.0001 --seed " + to_string(int(seed));
            cmd_list.push_back(v_command);
        }

        ThreadPool poolshapeit(stoi(v_threads));
        std::vector< std::future<int> > results_shapeit;
        for (auto item : cmd_list)
        {
            results_shapeit.emplace_back(
                    poolshapeit.enqueue([item]
                        {
                            string out = GetStdoutFromCommand(item);
                            return 1;
                        })
                );

        }
        for(auto && result: results_shapeit){result.get();} // waiting for all threads
    
    }
    
    











    map <string,string> best_pair_haplo;
    map <string,float> best_pair_ratio;

    if (v_recall == 0) {

        screen_message(screen_size, 1, " > Comparing haplotypes ...",1,v_debug); 
        map <pair<string,string>,int> haplodata;
        vector <string> snps_from_shapeit;

        for (int a = 1; a <= v_replicates; a++)
        {

            string outshapeit = shapeit_dir + "/shapeit_run" + to_string(a) + ".vcf";
            
            ifstream vcf (outshapeit.c_str());
            vector <string> sampledata;
            map <string,string> vector_h1;
            map <string,string> vector_h2;
            for( std::string line; getline( vcf, line ); )
            {
                if (line == "") {continue;}
                if (line.substr(0,2) == "##"){continue;}
                if (line.substr(0,2) == "#C") {boost::split(sampledata,line,boost::is_any_of("\t")); continue;}

                vector <string> data;
                boost::split(data,line,boost::is_any_of("\t"));

                if (a == 1)
                {
                    string info = data[0] + "\t" + data[1] + "\t" + data[2] + "\t" + data[3]  + "\t" + data[4];
                    snps_from_shapeit.push_back(info);
                }


                for (int a = 9; a < data.size(); a++)
                {
                    string genotype = data[a];
                    string sample = sampledata[a];
                    vector <string> alleles;
                    boost::split(alleles,genotype,boost::is_any_of("|"));
                    vector_h1[sample].append("," + alleles[0]);
                    vector_h2[sample].append("," + alleles[1]);
                }

            }
            vcf.close();


            for (auto item : vector_h1)
            {
                string sample = item.first;
                string haplotype_h1 = item.second;
                string haplotype_h2 = vector_h2[sample];
                vector <string> haplotypes;
                haplotypes.push_back(haplotype_h1);
                haplotypes.push_back(haplotype_h2);
                std::sort(haplotypes.begin(), haplotypes.end());
                string haplo = haplotypes[0] + "\n" + haplotypes[1];

                pair <string,string> j;
                j = make_pair(sample,haplo);

                haplodata[j]++;

            }
        }



        
        for (auto item : haplodata)
        {
            string sample = item.first.first;
            int hits = item.second;
            string haplo = item.first.second;
            float ratio = float(hits) / float(v_replicates);

            if (best_pair_ratio.find(sample) == best_pair_ratio.end())
            {
                best_pair_ratio[sample] = ratio;
                best_pair_haplo[sample] = haplo;
                continue;
            }

            if (best_pair_ratio[sample] < ratio)
            {
                best_pair_ratio[sample] = ratio;
                best_pair_haplo[sample] = haplo;
                continue;
            }

        }
   

   






        screen_message(screen_size, 1, " > Writing best haplotypes ...",1,v_debug); 

        map <pair<int,string>,string> final_genotypes;
        for (auto sample : sample_list)
        {
            vector <string> haplos;
            boost::split(haplos,best_pair_haplo[sample.first],boost::is_any_of("\n"));
            vector <string> alleles_h1;
            vector <string> alleles_h2;
            boost::split(alleles_h1,haplos[0],boost::is_any_of(","));
            boost::split(alleles_h2,haplos[1],boost::is_any_of(","));
            
            for (int a = 1; a <= snps_from_shapeit.size(); a++) {
                
                string gen = alleles_h1[a] + "|" + alleles_h2[a];

                pair <int,string> k;
                k = make_pair(a,sample.first);

                final_genotypes[k] = gen;
            }
        }




        ofstream bestvcf;
        file = master_output + "/fake_diploid_for_shapeit4.phased.vcf";

        bestvcf.open (file.c_str());
        bestvcf << "##fileformat=VCFv4.2" << endl;
        bestvcf << "##contig=<ID=fakechr,length=3000000>" << endl;
        bestvcf << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
        bestvcf << "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">" << endl;
        

        lineout = "";
        lineout.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for (auto sample : sample_list)
        {
            lineout.append("\t" + sample.first);
        }
        bestvcf << lineout << endl;

        int loop = 1;
        for (auto pos : snps_from_shapeit)
        {
        
            lineout = "";
            lineout.append(pos);
            lineout.append("\t.\t.\t.\tGT");

            for (auto sample : sample_list)
            {
                pair <int,string> k;
                k = make_pair(loop,sample.first);
                lineout.append("\t" + final_genotypes[k]);
            }
            bestvcf << lineout << endl;
            loop++;
        }
        bestvcf.close();




        

    }

        string outmulti = master_output + "/fake_diploid_for_shapeit4.phased.multi.vcf";
        v_command = v_bcftools + " norm -m+any -o " + outmulti + " " + file;
        if (v_recall == 0) {out = GetStdoutFromCommand(v_command);}
        if (! fileExists(outmulti )) {cout << "Something went wrong!" << endl << endl; return;}


 
    
        map <pair<string,string>,string> absence_vector;
        vector <string> sampledata;
        ifstream shapeitvcf (outmulti.c_str());
        for( std::string line; getline( shapeitvcf, line ); )
        {
            if (line == "") {continue;}
            if (line.substr(0,2) == "##"){continue;}
            if (line.substr(0,2) == "#C") {boost::split(sampledata,line,boost::is_any_of("\t")); continue;}

            vector <string> data;
            boost::split(data,line,boost::is_any_of("\t"));

            if (data[3] != "<PRESENCE>") {continue;}
            string gene = data[2];

            for (int a = 9; a < data.size(); a++)
            {
                pair <string,string> j;
                j = make_pair(gene,sampledata[a]);
                absence_vector[j] = data[a];
            }

        }
        shapeitvcf.close();






        screen_message(screen_size, 1, " > Writing final phased VCFs ...",1,v_debug); 
        map <pair<string,int>,string> fakedata;

        
        string sampleline = "";
        map <string,int> chrlist;

        shapeitvcf.open (outmulti.c_str());
        for( std::string line; getline( shapeitvcf, line ); )
        {
            if (line == "") {continue;}
            if (line.substr(0,2) == "##") {continue;}
            if (line.substr(0,2) == "#C") {sampleline = line; continue;}
            vector <string> data;
            boost::split(data,line,boost::is_any_of("\t"));
            string realpos = snp_correct_position[stoi(data[1])];
            if (realpos == "") {continue;}
            vector <string> sub;
            boost::split(sub,realpos,boost::is_any_of("\t"));
    

            pair <string,int> j;
            j = make_pair(sub[0],stoi(sub[1]));
            chrlist[sub[0]] = 1;

            fakedata[j] = line;

        }
        shapeitvcf.close();


        for (auto chr : chrlist)
        {

            if (chr.first == "fakechr") {continue;}
            string outchr = master_output + "/fake_diploid_phased." + chr.first + ".vcf";

            ofstream chrvcf;
            chrvcf.open (outchr.c_str());

            chrvcf << "##fileformat=VCFv4.2" << endl;
            chrvcf << "##reference=hg38DH.fa" << endl;
            chrvcf << "##contig=<ID=chr19,length=58617616>" << endl;
            chrvcf << "##contig=<ID=chr19_KI270890v1_alt,length=184499>" << endl;
            chrvcf << "##contig=<ID=chr19_KI270921v1_alt,length=282224>" << endl;
            chrvcf << "##contig=<ID=chr19_KI270923v1_alt,length=189352>" << endl;
            chrvcf << "##contig=<ID=chr19_KI270938v1_alt,length=1066800>" << endl;
            chrvcf << "##contig=<ID=chr6,length=170805979>" << endl;
            chrvcf << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
            chrvcf << sampleline << endl;

            for (auto pos : fakedata)
            {
                
                if (pos.first.first != chr.first) {continue;}
                string lineout = "";

                lineout.append(chr.first + "\t" + to_string(pos.first.second));
                vector <string> sub;
                boost::split(sub,pos.second,boost::is_any_of("\t"));
                for (int a = 2; a < sub.size();a++)
                {
                    lineout.append("\t" + sub[a]);
                }
                chrvcf << lineout << endl;

            }
            chrvcf.close();
        }






	ThreadPool poolfasta(stoi(v_threads));
	std::vector< std::future<int> > results_fasta;
    int count = 0;

    screen_message(screen_size, 1, " > Writing fasta files ...",1,v_debug); 
    for (auto gene : fake_gene_positions)
    {
        if (valid_genes.find(gene.first) == valid_genes.end()) {continue;}
        count++;

        string ref = chr_reference[gene.first];
        string currentgene = gene.first;
        
        string outchr = master_output + "/fake_diploid_phased." + ref + ".vcf";
        string reffas = v_db + "/reference/hg38/" + ref + ".fasta";
        
        string bed =  "";
        if (v_exons == 1) {bed =  v_db + "/haplotypes/bed/" + gene.first + ".exons.bed";}
        if (v_full == 1) {bed =  v_db + "/haplotypes/bed/" + gene.first + ".exons.bed";}
        if (v_cds == 1) {bed =  v_db + "/haplotypes/bed/" + gene.first + ".CDS.bed";}

        
       string outfas = master_output + "/" + gene.first + ".fas";

       string cmd = "";
       results_fasta.emplace_back(
				poolfasta.enqueue([count,cmd, currentgene, outfas, bed, reffas, outchr, master_output]
					{
                       

                        fasta_exons (bed, outchr, reffas, currentgene, outfas, master_output);


                        if (currentgene == "KIR2DS5")
                        {
                            string mvcmd = "mv " + outfas + " " + outfas + ".tmp";
                            string outcmd = GetStdoutFromCommand(mvcmd);

                            ofstream outrev;
                            outrev.open (outfas.c_str());

                            string inp = outfas + ".tmp";
                            ifstream transeq (inp.c_str());
                            for( std::string line; getline( transeq, line ); )
                            {
                                if (line == "") {continue;}
                                if (line.substr(0,1) == ">") {outrev << line << endl; continue;}
                                string revseq = reverse_and_complement(line);
                                outrev << revseq << endl;
                            }
                            transeq.close();
                            outrev.close();
                            removefile(outfas + ".tmp",v_debug);

                        }


                        if (currentgene == "KIR3DP1")
                        {
                            string mvcmd = "mv " + outfas + " " + outfas + ".tmp";
                            string outcmd = GetStdoutFromCommand(mvcmd);

                            ofstream outrev;
                            outrev.open (outfas.c_str());

                            string inp = outfas + ".tmp";
                            ifstream transeq (inp.c_str());
                            for( std::string line; getline( transeq, line ); )
                            {
                                if (line == "") {continue;}
                                if (line.substr(0,1) == ">") {outrev << line << endl; continue;}
                                
                                
                                string seq = line;

                                if (seq.find("DEL:KIR3DP1:EXON2") != std::string::npos) {
                                    size_t start;
                                    start = seq.find("CATGGCGTGTGTTG");
                                    size_t end;
                                    end = seq.find("GTGGTCAGGA");
                                    if (end == std::string::npos) {end = seq.find("GTGGTCAGAA");}

                                    if ((start != std::string::npos) && (end != std::string::npos))     
                                    {
                                        string subA = seq.substr(0,start+14);
                                        string subB = seq.substr(end);
                                        seq = subA + subB;
                                    }                               
                                }
                                
                                outrev << seq << endl;
                            }
                            transeq.close();
                            outrev.close();
                            removefile(outfas + ".tmp",v_debug);

                        }



						return 1;
					})
			);

    }
    for(auto && result: results_fasta){result.get();} // waiting for all threads









    
    screen_message(screen_size, 1, " > Comparing fasta files with the database ...",1,v_debug); 

    map <pair<string,string>,string> final_results;
    map <string,int> sample_names;
    for (auto gene : fake_gene_positions)
    {
        if (valid_genes.find(gene.first) == valid_genes.end()) {continue;}
        string currentgene = gene.first;
        string imgt = v_db + "/haplotypes/imgt_db/" + currentgene + "_exons.fasta";
        map <string,string> known_sequences_name;
        map <string,string> known_sequences_seq;
        ifstream imgtdata (imgt.c_str());
        string id = "";
        for( std::string line; getline( imgtdata, line ); )
        {
            if (line == "") {continue;}
            if (line.substr(0,1) == ">")
            {
                vector <string> data;
                boost::split(data,line,boost::is_any_of(" "));
                id = data[1];
                continue;
            }
            known_sequences_name[id].append(line);
        }
        imgtdata.close();

        for (auto item : known_sequences_name)
        {
            known_sequences_seq[item.second].append("," + item.first);
        }





        string outfas = master_output + "/" + gene.first + ".fas";
        ifstream ourfas (outfas.c_str());
        map <string,int> ourseqs;
        map <pair<string,string>,string> ourdata;
  
        id = "";
        string vec = "";
        for( std::string line; getline( ourfas, line ); )
        {
            if (line == "") {continue;}
            if (line.substr(0,1) == ">")
            {
                vector <string> data;
                string subid = line.substr(1);
                boost::split(data,subid,boost::is_any_of("_"));
                id = data[0];
                vec = data[1];
                continue;
            }
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            ourseqs[line]++;

            pair <string,string> j;
            j = make_pair(id,vec);
            ourdata[j] = line;
            
        }
        ourfas.close();



        string sequences_names = master_output + "/" + gene.first + ".names.fas";
        ofstream seqname;
        seqname.open (sequences_names.c_str());
        int newseq = 1;
        map <string,int> versions;
        map <string,string> converted_seq_names;
        for (auto item : ourseqs)
        {
            string seq = item.first;
            string name = "";
            map <string,int> used;
            used.clear();

            if (known_sequences_seq.find(seq) != known_sequences_seq.end()) 
            {
                string subname = known_sequences_seq[seq];
                if (v_full == 0) {
                    vector <string> data;
                    boost::split(data,subname,boost::is_any_of("*"));
                    if (data[1].size() > 5) {data[1] = data[1].substr(0,5);}
                    subname = data[0] + "*" + data[1];
                }
                if (used.find(subname) == used.end()) {
                    name = "," + subname;
                    used[subname] = 1;
                }
                
            }

            if (name == "")
            {
                for (auto item : known_sequences_seq)
                {
                    if (seq.find(item.first) != std::string::npos) {
                        
                        string subname = known_sequences_seq[item.first];
                        if (v_full == 0) {
                            vector <string> data;
                            boost::split(data,subname,boost::is_any_of("*"));
                            if (data[1].size() > 5) {data[1] = data[1].substr(0,5);}
                            subname = data[0] + "*" + data[1];
                        }
                        if (used.find(subname) == used.end()) {
                            name = "," + subname;
                            used[subname] = 1;
                        }
                    }
                }
            }

            if (name == "")
            {
                for (auto item : known_sequences_seq)
                {
                    if (item.first.find(seq) != std::string::npos) {
                        string subname = known_sequences_seq[item.first];
                        if (v_full == 0) {
                            vector <string> data;
                            boost::split(data,subname,boost::is_any_of("*"));
                            if (data[1].size() > 5) {data[1] = data[1].substr(0,5);}
                            subname = data[0] + "*" + data[1];
                        }
                        if (used.find(subname) == used.end()) {
                            name = "," + subname;
                            used[subname] = 1;
                        }
                     }
                }
            }

            if (name == "") {name = "," + gene.first + ".new." + to_string(newseq); newseq++;}
            boost::replace_all(name, ",,", ",");

            name = name.substr(1);
            versions[name]++;
            if (versions[name] > 1) {
                name.append("(version " + to_string(versions[name]) + ")");
            }
            seqname << ">" << name << endl;           
            seqname << seq << endl;
            converted_seq_names[seq] = name;
        }

        seqname.close();

        ofstream outdb;
        string outdbfile = master_output + "/" + gene.first + ".db.txt";
        outdb.open (outdbfile.c_str());
        outdb << "Sample\tvector\tcopy_number\twarning\tsequence_name\tsequence" << endl;

        for (auto item : ourdata)
        {
            string sample = item.first.first;
            string vec = item.first.second;
            string seq = item.second;
            string name = converted_seq_names[seq];
            string warn = "none";

            pair <string,string> k;
            k = make_pair(gene.first,sample);
            string cn = copy_numbers[k];


            pair <string,string> j;
            j = make_pair(gene.first,sample);
            string absence = absence_vector[j];

            if (cn == "0") {warn = "zero copies for this gene";seq = "null_allele"; name = gene.first + "*null";}
            if (cn == "1") 
            {
                warn = "one copy for this gene. One of the haplotypes is *null";
                if ((absence == "0|1") && (vec == "h2")){seq = "null_allele"; name = gene.first + "*null";}
                if ((absence == "1|0") && (vec == "h1")){seq = "null_allele"; name = gene.first + "*null";}
            }
            if (cn == "3") {warn = "three copies for this gene, reduced to two.";}
            if (cn == "4") {warn = "four copies for this gene. Everything was imputed. Do not consider this result.";}

            outdb << sample << "\t" << vec << "\t" << cn << "\t" << warn << "\t" << name << "\t" << seq << endl;

            pair <string,string> l;
            l = make_pair(sample + "\t" + vec, gene.first);
            final_results[l] = name;
            sample_names[sample] = 1;

        }
     
        outdb.close();
    }





    ofstream outdb;
    string outdbfile = master_output + "/merged_one_line_per_chromosome.db.txt";
    outdb.open (outdbfile.c_str());
    outdb << "Sample\tvector\twarning\thaplotype_P_value";
    for (auto gene : list_of_genes_expected_positions)
    {
        if (valid_genes.find(gene) == valid_genes.end()) {continue;}
        outdb << "\t" << gene;
    }
    outdb << endl;


    for (auto item : sample_names)
    {
        
        outdb << item.first << "\th1\t";

        string warn = "";

        for (auto gene : list_of_genes_expected_positions)
        {
            if (valid_genes.find(gene) == valid_genes.end()) {continue;}
            pair <string,string> k;
            k = make_pair(gene,item.first);
            string cn = copy_numbers[k];
            if (cn == "3") {warn.append(", 3 copies of " + gene);}
            if (cn == "4") {warn.append(", 4 copies of " + gene);}
        }
        if (warn != "") {warn = warn.substr(2);}
        if (warn == "") {warn = "none";}
        outdb << warn;
        outdb << "\t" << best_pair_ratio[item.first];


        for (auto gene : list_of_genes_expected_positions)
        {
            if (valid_genes.find(gene) == valid_genes.end()) {continue;}
            pair <string,string> l;
            l = make_pair(item.first + "\th1", gene);
            outdb << "\t" << final_results[l];
        }
        outdb << endl;

        outdb << item.first << "\th2\t" << warn;
        outdb << "\t" << best_pair_ratio[item.first];
        for (auto gene : list_of_genes_expected_positions)
        {
            if (valid_genes.find(gene) == valid_genes.end()) {continue;}
            pair <string,string> l;
            l = make_pair(item.first + "\th2", gene);
            outdb << "\t" << final_results[l];
        }
        outdb << endl;

    }
 
    outdb.close();





    outdbfile = master_output + "/merged_two_chromosomes_perl_line.db.txt";
    outdb.open (outdbfile.c_str());
    outdb << "Sample\twarning\thaplotype_P_value";
    for (auto gene : list_of_genes_expected_positions)
    {
        if (valid_genes.find(gene) == valid_genes.end()) {continue;}
        outdb << "\t" << gene << " (h1)";
        outdb << "\t" << gene << " (h2)";
    }
    outdb << endl;


    for (auto item : sample_names)
    {
        
        outdb << item.first << "\t"; 

        string warn = "";

        for (auto gene : list_of_genes_expected_positions)
        {
            if (valid_genes.find(gene) == valid_genes.end()) {continue;}
            pair <string,string> k;
            k = make_pair(gene,item.first);
            string cn = copy_numbers[k];
            if (cn == "3") {warn.append(", 3 copies of " + gene);}
            if (cn == "4") {warn.append(", 4 copies of " + gene);}
        }
        if (warn != "") {warn = warn.substr(2);}
        if (warn == "") {warn = "none";}
        outdb << warn;
        outdb << "\t" << best_pair_ratio[item.first];


        for (auto gene : list_of_genes_expected_positions)
        {
            if (valid_genes.find(gene) == valid_genes.end()) {continue;}
            pair <string,string> l;
            l = make_pair(item.first + "\th1", gene);
            outdb << "\t" << final_results[l];

            pair <string,string> k;
            k = make_pair(item.first + "\th2", gene);
            outdb << "\t" << final_results[k];

        }

        outdb << endl;
 
    }

 
    outdb.close();



    

}