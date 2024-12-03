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

#include "setup.hpp"
#include "external.hpp"
#include "functions.hpp"


void main_setup() {
    
    screen_message (screen_size, 0, "", 1, v_quiet);
    v_message = "Program:   " + Program_name + "::setup";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "Version:   " + Program_version + ", " + Program_date;
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);


    v_message = "kir-mapper needs some third-party software.";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);
    v_message = "Here is the list and their tested versions.";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    
    string list = "  Samtools, version(s)";
    for (auto item : samtools_versions)
    {
        list.append (" " + item.first);
    }
    screen_message (screen_size, 0, list, 1, v_quiet);



    list = "  Bcftools, version(s)";
    for (auto item : bcf_versions)
    {
        list.append (" " + item.first);
    }
    screen_message (screen_size, 0, list, 1, v_quiet);
    


    list = "  BWA, version(s)";
    for (auto item : bwa_versions)
    {
        list.append (" " + item.first);
    }
    screen_message (screen_size, 0, list, 1, v_quiet);



    list = "  Freebayes, version(s)";
    for (auto item : freebayes_versions)
    {
        list.append (" " + item.first);
    }
    screen_message (screen_size, 0, list, 1, v_quiet);
    


    list = "  WhatsHap, version(s)";
    for (auto item : whats_versions)
    {
        list.append (" " + item.first);
    }
    screen_message (screen_size, 0, list, 1, v_quiet);

    

    list = "  R, version(s)";
    for (auto item : R_versions)
    {
        list.append (" " + item.first);
    }
    screen_message (screen_size, 0, list, 1, v_quiet);
    

    
    list = "  Picard tools, any version";
    screen_message (screen_size, 0, list, 1, v_quiet);

    
    list = "  STAR RNA-seq aligner, version(s)";
    for (auto item : star_versions)
    {
        list.append (" " + item.first);
    }
    screen_message (screen_size, 0, list, 1, v_quiet);


    list = "  Shapeit4";
    screen_message (screen_size, 0, list, 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);


    screen_message (screen_size, 0, "", 1, v_quiet);
    cout << "Shall we proceed with the setup process? yes (y) or close/cancel (c): ";
    string res = "";
    cin >> res;
    res.erase(std::remove(res.begin(), res.end(), '\n'), res.end());
    if ((res != "y") && (res != "Y")) 
    { 
        screen_message (screen_size, 0, "", 1, v_quiet);
        return;
    }

   



    
 // ****** R ******

    screen_message (screen_size, 0, "", 1, v_quiet);
    string command_R = "which R";
    cout << "**** Configuring R and libraries ***" << endl;
    string outR = GetStdoutFromCommand(command_R);
    outR.erase(std::remove(outR.begin(), outR.end(), '\n'), outR.end());
    if (outR == "") {
        cout << "You need R to continue. Please install R and run 'kir-mapper setup' again." << endl; 
        screen_message (screen_size, 0, "", 1, v_quiet);
        return;
    }
    
    string Rfile = "kir-mapper_setup.R";
    ofstream myfileR;
    myfileR.open (Rfile);
    myfileR << "options(warn=-1)" << endl;
	myfileR << "list.of.packages <- c('ggplot2','plotly','htmlwidgets','stringr','shiny','shinyjs','purrr','readr','dplyr','tidyr','DT','forcats','shinyalert','pacman')" << endl;
	myfileR << "new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,'Package'])]" << endl;

    myfileR << "if(length(new.packages) == 0){cat ('All R packages are installed!\\n')}" << endl;
    myfileR << "if(length(new.packages) > 0){" << endl;
    myfileR << "cat ('List of missing packages:\\n')" << endl;
    myfileR << "for (i in new.packages){" << endl;
    myfileR << "cat(paste(i,'\\n'))}" << endl;
    myfileR << "cat ('Trying to install...\\n')" << endl;
    myfileR << "install.packages(new.packages, repos = 'http://cran.us.r-project.org')" << endl;
    myfileR << "}" << endl;

    myfileR.close();
    command_R = "Rscript kir-mapper_setup.R";
    outR = GetStdoutFromCommand(command_R);
    cout << outR << endl << endl;
    removefile(Rfile.c_str(),v_debug);



    
// ********** SAMTOOLS  ************
    screen_message (screen_size, 0, "", 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);
    string program = "Samtools";
    string command = "which samtools";
    string samtools = "";
    cout << "**** Configuring " << program << " ***" << endl;
    cout << "Versions < 1.18 won't work." << endl;
    string item = GetStdoutFromCommand(command);
    item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
    
samtools_start:
    if (item != "")
    {
        cout << program << " path (automatic detection): " << item << endl;
        cout << "Use this one (y), or change (c), or quit (q): ";
        string res = "";
        cin >> res;
        res.erase(std::remove(res.begin(), res.end(), '\n'), res.end());
        int valid = 0;
        if ((res == "y") || (res == "Y")) {valid = 1; samtools = item;}
        if ((res == "c") || (res == "C")) {valid = 1; item = "";}
        if ((res == "q") || (res == "Q")) {return;}
        if (valid == 0) {goto samtools_start;}
    }

samtools:
    if (item == "")
    {
        while (! fileExists(item)) {
            cout << "Please write the " << program << " path, or drag it below (or q to quit): " << endl;
            cin >> item;
            item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
            if ((item == "q") || (item == "Q")) {cout << endl << endl; return;}
            if (! fileExists(item))
            {
                cout << "Invalid file or program." << endl;
            }
            else {cout << program << " path: " << item << endl;}
        }
    }
    samtools = item;

    item = GetStdoutFromCommand(samtools + " --version");
    vector <string> data;
    boost::split(data,item,boost::is_any_of("\n"));
    vector <string> sub;
    boost::split(sub,data[0],boost::is_any_of(" "));
    if (samtools_versions.find(sub[1]) != samtools_versions.end())
        {
            cout << "This version (" << sub[1] << ") is compatible with kir-mapper" << endl;
        }
    else {
        cout << "This version has not been tested with kir-mapper\nContinue anyway (y), quit (q), or return (r)?";
        cin >> item;
        item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
        int valid = 0;
        if ((item == "q") || (item == "Q")) {valid = 1; cout << endl << endl; return; }
        if ((item == "r") || (item == "R")) {valid = 1;item = ""; goto samtools;}
        if ((item == "y") || (item == "Y")) {valid = 1;}
        if (valid == 0) {goto samtools;}
    }

    

    
    
    
    
    
    
// ********** BWA  ************
    program = "BWA";
    command = "which bwa";
    string bwa = "";
    cout << endl;
    cout << endl;
    cout << "**** Configuring " << program << " ***" << endl;
    item = GetStdoutFromCommand(command);
    item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
 
bwa_start:
    if (item != "")
        {
            cout << program << " path (automatic detection): " << item << endl;
            cout << "Use this one (y), or change (c): ";
            string res = "";
            cin >> res;
            int valid = 0;
            if ((res == "y") || (res == "Y")) {valid = 1; bwa = item;}
            if ((res == "c") || (res == "C")) {valid = 1;item = "";}
            if (valid == 0) {goto bwa_start;}
        }
    
bwa:
    if (item == "")
        {
            while (! fileExists(item)) {
                cout << "Please write the " << program << " path, or drag it bellow (or q to quit): " << endl;
                cin >> item;
                item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
                if ((item == "q") || (item == "Q"))  {cout << endl << endl; return;}
                if (! fileExists(item))
                    {
                        cout << "Invalid file or program." << endl;
                    }
                else {cout << program << " path: " << item << endl << endl;}
            }
        }
    
    bwa = item;
    item = GetStdoutFromCommand(bwa);
    data.clear();
    boost::split(data,item,boost::is_any_of("\n"));
    sub.clear();
    boost::split(sub,data[2],boost::is_any_of(" "));
    boost::split(data,sub[1],boost::is_any_of("-"));
    if (bwa_versions.find(data[0]) != bwa_versions.end())
        {
            cout << "This version (" << data[0] << ") is compatible with kir-mapper" << endl;
        }
    else {
        cout << "This version has not been tested with kir-mapper\nContinue anyway (y), quit (q), or return (r)? ";
        cin >> item;
        item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
        int valid = 0;
        if ((item == "q") || (item == "Q")) {valid = 1; cout << endl << endl; return; }
        if ((item == "r") || (item == "R")) {valid = 1; item = ""; goto bwa;}
        if ((item == "y") || (item == "Y")) {valid = 1;}
        if (valid == 0) {goto bwa;}
    }

    

    
    
    
// ********** BCFTOOLS  ************
    program = "Bcftools";
    command = "which bcftools";
    string bcf = "";
    cout << endl;
    cout << endl;
    cout << "**** Configuring " << program << " ***" << endl;
    item = GetStdoutFromCommand(command);
    item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
 
bcf_start:
    if (item != "")
        {
            cout << program << " path (automatic detection): " << item << endl;
            cout << "Use this one (y), or change (c): ";
            string res = "";
            cin >> res;
            int valid = 0;
            if ((res == "y") || (res == "Y")) {valid = 1; bcf = item;}
            if ((res == "c") || (res == "C")) {valid = 1;item = "";}
            if (valid == 0) {goto bcf_start;}
        }
    
bcf:
    if (item == "")
        {
            while (! fileExists(item)) {
                cout << "Please write the " << program << " path, or drag it bellow (or q to quit): " << endl;
                cin >> item;
                item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
                if ((item == "q") || (item == "Q"))  {cout << endl << endl; return;}
                if (! fileExists(item))
                    {
                        cout << "Invalid file or program." << endl;
                    }
                else {cout << program << " path: " << item << endl << endl;}
            }
        }
    
    bcf = item;
    item = GetStdoutFromCommand(bcf + " --version");
    data.clear();
    boost::split(data,item,boost::is_any_of("\n"));
    sub.clear();
    boost::split(sub,data[0],boost::is_any_of(" "));
    if (bcf_versions.find(sub[1]) != bcf_versions.end())
    {
        cout << "This version (" << data[0] << ") is compatible with kir-mapper" << endl;
    }
    else {
        cout << "This version has not been tested with kir-mapper\nContinue anyway (y), quit (q), or return (r)? ";
        cin >> item;
        item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
        int valid = 0;
        if ((item == "q") || (item == "Q")) {valid = 1; cout << endl << endl; return; }
        if ((item == "r") || (item == "R")) {valid = 1; item = ""; goto bcf;}
        if ((item == "y") || (item == "Y")) {valid = 1;}
        if (valid == 0) {goto bcf;}
    }

    
    
    
    
    
    
// ********** Database  ************
    
    item = "";
    program = "Database";
    string db = "";
    cout << endl;
    cout << endl;
    cout << "**** Configuring " << program << " ***" << endl;
    
    if (item == "")
        {
            while (! fileExists(item)) {
                cout << "Please write the " << program << " path, or drag it below (or q to quit): " << endl;
                cin >> item;
                item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
                if ((item == "q") || (item == "Q"))  {cout << endl << endl; return;}
                if (! fileExists(item))
                    {
                        cout << "Invalid file or folder." << endl;
                    }
                else {cout << program << " path: " << item << endl << endl;}
            }
        }
    db = item;
    
    

    
    
    
    
    
    
string version = "";
// ********** Freebayes  ************
program = "Freebayes";
command = "which freebayes";
string freebayes = "";
cout << endl;
cout << endl;
cout << "**** Configuring " << program << " ***" << endl;

item = GetStdoutFromCommand(command);
try
{
    item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
}
    catch (const std::exception& e)
{ item = "";}

freebayes_start:
    if (item != "")
        {
            cout << program << " path (automatic detection): " << item << endl;
            cout << "Use this one (y), or change (c), or quit (q): ";
            string res = "";
            cin >> res;
            int valid = 0;
            if ((res == "y") || (res == "Y")) {valid = 1; freebayes = item;}
            if ((res == "c") || (res == "C")) {valid = 1;item = "";}
            if ((res == "q") || (res == "Q")) {cout << endl << endl; return; }
            if (valid == 0) {goto freebayes_start;}
        }
    
freebayes:
    if (item == "")
        {
            while (! fileExists(item)) {
                cout << "Please write the " << program << " path, or drag it below (or q to quit): " << endl;
                cin >> item;
                item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
                if ((item == "q") || (item == "Q")) {cout << "Quiting!" << endl; return;}
                if (! fileExists(item))
                    {
                        cout << "Invalid file or program." << endl;
                    }
                else {cout << program << " path: " << item << endl << endl;}
            }
        }
    
    freebayes = item;
    
    try
    {
        item = GetStdoutFromCommand(freebayes);
    }
        catch (const std::exception& e)
    { item = ""; freebayes = "";}
    
    data.clear();
    
    try
    {
        boost::split(data,item,boost::is_any_of("\n"));
    }
        catch (const std::exception& e)
    { data.clear();}
    
    
    version = "";
    for (auto value : data)
    {
        vector <string> version_data;
        boost::split(version_data,value,boost::is_any_of(":"));
        if (version_data[0] == "version")
        {
            version = version_data[1];
            version.erase(remove(version.begin(), version.end(), ' '), version.end());
        }
    }
    
    if (freebayes_versions.find(version) != freebayes_versions.end())
        {
            cout << "This version (" << version << ") is compatible with kir-mapper" << endl;
        }
    else {
        cout << "This version has not been tested with kir-mapper\nContinue anyway (y), quit (q), or return (r)? ";
        cin >> item;
        item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
        int valid = 0;
        if ((res == "q") || (res == "Q")) {valid = 1; cout << "Quiting!" << endl; return; }
        if ((res == "r") || (res == "R")) {valid = 1; item = ""; goto freebayes;}
        if ((res == "y") || (res == "Y")) {valid = 1;}
        if (valid == 0) {goto freebayes;}
    }
 
    
    
    
    
    
    
// ********** WHATSHAP  ************
    
    program = "WhatsHap";
    command = "which whatshap";
    string whats = "";
    cout << endl;
    cout << endl;
    cout << "**** Configuring " << program << " ***" << endl;
    cout << program << " is mandatory for kir-mapper genotype." << endl;
    cout << "Shall we configure " + program + "? Yes (y), ignore it (i), or quit (q): ";

    cin >> res;
    res.erase(std::remove(res.begin(), res.end(), '\n'), res.end());
    if ((res == "i") || (res == "i")) {whats = "DISABLED";}
    if ((res == "q") || (res == "Q")) {cout << endl << endl; return; }

    if (((res == "y") || (res == "Y")) && (whats != "DISABLED"))
    {
    
        item = GetStdoutFromCommand(command);
        try
        {
            item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
        }
            catch (const std::exception& e)
        { item = "";}

    whats_start:
        if (item != "")
            {
                cout << program << " path (automatic detection): " << item << endl;
                cout << "Use this one (y), or change (c), or quit (q): ";
                string res = "";
                cin >> res;
                int valid = 0;
                if ((res == "y") || (res == "Y")) {valid = 1; whats = item;}
                if ((res == "c") || (res == "C")) {valid = 1; item = "";}
                if ((res == "q") || (res == "Q")) {cout << endl << endl; return; }
                if (valid == 0) {goto whats_start;}
            }
        
    whats:
        if (item == "")
            {
                while (! fileExists(item)) {
                    cout << "Please write the " << program << " path, or drag it below (or q to quit): " << endl;
                    cin >> item;
                    item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
                    if ((item == "q") || (item == "Q")) {cout << endl << endl; return;}
                    if (! fileExists(item))
                        {
                            cout << "Invalid file or program." << endl;
                        }
                    else {cout << program << " path: " << item << endl << endl;}
                }
            }
        
        whats = item;
        
        item = GetStdoutFromCommand(whats + " --version");
        data.clear();
        boost::replace_all(item , "\n", "");

 //       boost::split(data,item,boost::is_any_of("\n"));
 //       sub.clear();
 //       boost::split(sub,data[0],boost::is_any_of(" "));
        
        if (whats_versions.find(item) != whats_versions.end())
            {
                cout << "This version (" << sub[1] << ") is compatible with kir-mapper" << endl;
            }
        else {
            cout << "This version has not been tested with kir-mapper\nContinue anyway (y), quit (q), or return (r)? ";
            cin >> item;
            item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
            int valid = 0;
            if ((item == "q") || (item == "Q")) {valid = 1; cout << "Quiting!" << endl; return; }
            if ((item == "r") || (item == "R")) {valid = 1; item = ""; goto whats;}
            if ((item == "y") || (item == "Y")) {valid = 1;}
            if (valid == 0) {goto whats;}
        }
    }

    
    
    
    
// ********** STAR  ************
    program = "STAR";
    command = "which STAR";
    string star = "";
    cout << endl;
    cout << endl;
    cout << "**** Configuring " << program << " ***" << endl;
    item = GetStdoutFromCommand(command);
    item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
 
star_start:
    if (item != "")
        {
            cout << program << " path (automatic detection): " << item << endl;
            cout << "Use this one (y), or change (c): ";
            string res = "";
            cin >> res;
            int valid = 0;
            if ((res == "y") || (res == "Y")) {valid = 1; star = item;}
            if ((res == "c") || (res == "C")) {valid = 1;item = "";}
            if (valid == 0) {goto star;}
        }
    
star:
    if (item == "")
        {
            while (! fileExists(item)) {
                cout << "Please write the " << program << " path, or drag it bellow (or q to quit): " << endl;
                cin >> item;
                item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
                if ((item == "q") || (item == "Q"))  {cout << endl << endl; return;}
                if (! fileExists(item))
                    {
                        cout << "Invalid file or program." << endl;
                    }
                else {cout << program << " path: " << item << endl << endl;}
            }
        }
    
    star = item;
    item = GetStdoutFromCommand(star + " --version");
    data.clear();
    boost::split(data,item,boost::is_any_of("\n"));
    if (star_versions.find(data[0]) != star_versions.end())
    {
        cout << "This version (" << data[0] << ") is compatible with kir-mapper" << endl;
    }
    else {
        cout << "This version has not been tested with kir-mapper\nContinue anyway (y), quit (q), or return (r)? ";
        cin >> item;
        item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
        int valid = 0;
        if ((item == "q") || (item == "Q")) {valid = 1; cout << endl << endl; return; }
        if ((item == "r") || (item == "R")) {valid = 1; item = ""; goto star;}
        if ((item == "y") || (item == "Y")) {valid = 1;}
        if (valid == 0) {goto star;}
    }

    









// ********** SHAPEIT  ************
    
    program = "Shapeit4";
    command = "which shapeit4";
    string shapeit = "";
    cout << endl;
    cout << endl;
    cout << "**** Configuring " << program << " ***" << endl;
    cout << program << " is mandatory for kir-mapper haplotypes." << endl;
    cout << "Shall we configure " + program + "? Yes (y), ignore it (i), or quit (q): ";

    cin >> res;
    res.erase(std::remove(res.begin(), res.end(), '\n'), res.end());
    if ((res == "i") || (res == "i")) {shapeit = "DISABLED";}
    if ((res == "q") || (res == "Q")) {cout << endl << endl; return; }

    if (((res == "y") || (res == "Y")) && (shapeit  != "DISABLED"))
    {
    
        item = GetStdoutFromCommand(command);
        try
        {
            item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
        }
            catch (const std::exception& e)
        { item = "";}

    shapeit_start:
        if (item != "")
            {
                cout << program << " path (automatic detection): " << item << endl;
                cout << "Use this one (y), or change (c), or quit (q): ";
                string res = "";
                cin >> res;
                int valid = 0;
                if ((res == "y") || (res == "Y")) {valid = 1; shapeit = item;}
                if ((res == "c") || (res == "C")) {valid = 1; item = "";}
                if ((res == "q") || (res == "Q")) {cout << endl << endl; return; }
                if (valid == 0) {goto shapeit_start;}
            }
        
    shapeit:
        if (item == "")
            {
                while (! fileExists(item)) {
                    cout << "Please write the " << program << " path, or drag it below (or q to quit): " << endl;
                    cin >> item;
                    item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
                    if ((item == "q") || (item == "Q")) {cout << endl << endl; return;}
                    if (! fileExists(item))
                        {
                            cout << "Invalid file or program." << endl; goto shapeit;
                        }
                    else {cout << program << " path: " << item << endl << endl;}
                }
            }
        
        shapeit = item;

    }

    
    


    









    
    // ********** PICARD  ************
        
        program = "Picard Tools";
        command = "which picard";
        string picard = "";
        cout << endl;
        cout << endl;
        cout << "**** Configuring " << program << " ***" << endl;
        cout << program << " is optional for kir-mapper map." << endl;
        cout << "Shall we configure " + program + "? Yes (y), ignore it (i), or quit (q): ";

        cin >> res;
        res.erase(std::remove(res.begin(), res.end(), '\n'), res.end());
        if ((res == "i") || (res == "i")) {picard = "DISABLED";}
        if ((res == "q") || (res == "Q")) {cout << endl << endl; return; }
        if (((res == "y") || (res == "Y")) && (picard != "DISABLED"))
        {
        
            item = GetStdoutFromCommand(command);
            try
            {
                item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
            }
                catch (const std::exception& e)
            { item = "";}
			
			if (item == "")
			{
				command = "which picard.jar";
				item = GetStdoutFromCommand(command);
			}

        picard_start:
            if (item != "")
                {
                    cout << program << " path (automatic detection): " << item << endl;
                    cout << "Use this one (y), or change (c), or quit (q): ";
                    string res = "";
                    cin >> res;
                    int valid = 0;
                    if ((res == "y") || (res == "Y")) {valid = 1; picard = item;}
                    if ((res == "c") || (res == "C")) {valid = 1; item = "";}
                    if ((res == "q") || (res == "Q")) {cout << endl << endl; return; }
                    if (valid == 0) {goto picard_start;}
                }
            
        picard:
            if (item == "")
                {
                    while (! fileExists(item)) {
                        cout << "Please write the " << program << " path, or drag it bellow (or q to quit): " << endl;
                        cin >> item;
                        item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
                        if ((item == "q") || (item == "Q")) {cout << endl << endl; return;}
                        if (! fileExists(item))
                            {
                                cout << "Invalid file or program." << endl;
								
                            }
                        else {cout << program << " path: " << item << endl << endl;}
                    }
                }
            
            picard = item;
            
        }

        





    
    
    cout << endl << endl;
    cout << "**** Configuration list ***" << endl;
    cout << "Database:  " << db << endl;
    cout << "Samtools:  " << samtools << endl;
    cout << "Bcftools:  " << bcf << endl;
    cout << "BWA:       " << bwa << endl;
    cout << "WhatsHap:  " << whats << endl;
    cout << "Freebayes: " << freebayes << endl;
    cout << "STAR:      " << star << endl;
    cout << "Picard:    " << picard << endl;
    cout << "shapeit4   " << shapeit << endl;
    cout << endl;
    
write:
    cout << "Confirm writing the configuration file? (y or n) ";
    cin >> item;
    item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
    int valid = 0;
    if (item == "n") {valid = 1; cout << endl << endl; return; }
    if (item == "y") {valid = 1;}
    if (valid == 0) {goto write;}
    
    
    string homedir = getpwuid(getuid())->pw_dir;
    string configfile = homedir + "/.kir-mapper";
    
    ofstream myfile;
    myfile.open (configfile);
    myfile << "db=" + db + "\n";
    myfile << "samtools=" + samtools + "\n";
    myfile << "bcftools=" + bcf + "\n";
    myfile << "bwa=" + bwa + "\n";
    myfile << "whatshap=" + whats + "\n";
    myfile << "freebayes=" + freebayes + "\n";
    myfile << "picard=" + picard + "\n";
    myfile << "star=" + star + "\n";
    myfile << "shapeit4=" + shapeit + "\n";
    myfile.close();
    cout << "\n";
    return;
}
