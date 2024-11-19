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


void main_join() {
    
    
    int v_check = 0;
    if (v_output != "") {v_check = 1;}
 
    if (((v_db == "") || (v_check == 0)))
    {
        screen_message (screen_size, 0, "", 1, v_quiet);
        v_message = "Program:   " + Program_name + "::join";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        v_message = "Version:   " + Program_version + ", " + Program_date;
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        v_message = "Usage:     kir-mapper join -output map_output_folder <options>";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Mandatory options:", 1, v_quiet);
        screen_message (screen_size, 2, "-output      output folder (same as map and ncopy)", 1, v_quiet);
        screen_message (screen_size, 2, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Other options:", 1, v_quiet);
        screen_message (screen_size, 2, "-db          path to the kir-mapper database", 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        screen_message (screen_size, 2, "--full      full genotype, not only exons", 1, v_quiet);
        screen_message (screen_size, 2, "--quiet     quiet mode", 1, v_quiet);
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


    if (v_output == "") {
        v_message = "You must indicate a valid map output folder";
        warnings.push_back (v_message);
        precheck = 0;
    }

    if (! fileExists(v_output))
    {
        v_message = "You must indicate a valid map output folder";
        warnings.push_back (v_message);
        precheck = 0;
    }
    
    string script = v_db + "/ncopy/join.pl";
    if (! fileExists(script)) {v_message = "There is something wrong with this database";
        warnings.push_back (v_message);
        precheck = 0;
    }

    if (precheck == 0){ v_db = ""; main_join();return; }
    
    screen_message (screen_size, 0, "", 1, v_quiet);
    v_message = "Program:   " + Program_name + "::join";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "Version:   " + Program_version + ", " + Program_date;
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    if (v_exome == 1) {
        v_message = "Mode:      only exons";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
    }
    if (v_exome == 0) {
        v_message = "Mode:      full genotype";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
    }

    string cmd = "";
    if (v_exome == 1)
    {
        cmd = "perl " + script + " -o " + v_output + " -m cds";
    }
    if (v_exome == 0)
    {
        cmd = "perl " + script + " -o " + v_output + " -m full";
    }
    GetStdoutFromCommand(cmd);
    v_message = "Computation done.";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    
}


void main_combine() {
    
    
    int v_check = 0;
    if (v_output != "") {v_check = 1;}
 
    if (((v_db == "") || (v_check == 0)))
    {
        screen_message (screen_size, 0, "", 1, v_quiet);
        v_message = "Program:   " + Program_name + "::group";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        v_message = "Version:   " + Program_version + ", " + Program_date;
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        v_message = "Usage:     kir-mapper group -list list_of_kir-mapper_outputs.txt -output output_folder <options>";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Mandatory options:", 1, v_quiet);
        screen_message (screen_size, 2, "-output      output folder", 1, v_quiet);
        screen_message (screen_size, 2, "-list        txt file with all kir-mapper outputs to join", 1, v_quiet);
        screen_message (screen_size, 2, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Other options:", 1, v_quiet);
        screen_message (screen_size, 2, "-db          path to the kir-mapper database", 1, v_quiet);
        screen_message (screen_size, 2, "--quiet     quiet mode", 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        return;
    }
    
    
    int precheck = 1;
    
    if (v_output == "") {
        v_message = "You must indicate a valid map output folder";
        warnings.push_back (v_message);
        precheck = 0;
    }

    if (! fileExists(v_output))
    {
        v_command = "mkdir " + v_output;
        GetStdoutFromCommand(v_command);
        if (! fileExists(v_output)) {
            v_message = "You must indicate a valid map output folder";
            warnings.push_back (v_message);
            precheck = 0;
        }
    }
    
    if (! fileExists(v_list))
    {
        v_message = "You must indicate a TXT file with the list of outputs to combine.";
        warnings.push_back (v_message);
        precheck = 0;
    }
  
    if (precheck == 0){ v_db = ""; main_combine();return; }
    
    screen_message (screen_size, 0, "", 1, v_quiet);
    v_message = "Program:   " + Program_name + "::combine";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "Version:   " + Program_version + ", " + Program_date;
    screen_message (screen_size, 0, v_message, 1, v_quiet);

 
    vector <string> outputs;
    ifstream file(v_list.c_str());
    int valid = 1;
    for( std::string line; getline( file, line ); )
    {
        if (line == "") {continue;}
        outputs.push_back(line);
        if (! fileExists(line))
        {
            valid = 0;
        }

        if (! fileExists(line + "/map/"))
        {
            valid = 0;
        }
        
        if (! fileExists(line + "/ncopy/"))
        {
            valid = 0;
        }
    }
    file.close();
    
    if (valid == 0)
    {
        v_message = "Some outputs at the provided list are not valid.";
        warnings.push_back (v_message);
    }
  
    if (valid == 0){ v_db = ""; main_combine();return; }
    
    v_command = "mkdir " + v_output + "/map/";
    GetStdoutFromCommand(v_command);
    
    v_command = "mkdir " + v_output + "/ncopy/";
    GetStdoutFromCommand(v_command);

   
    v_message = "Copying map files ...";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    
    
    for (auto item : outputs)
    {
        v_message = "Copying map files from " + item + " ...";
        screen_message (screen_size, 0, v_message, 1, v_quiet);

        for (const auto & entry : fs::directory_iterator(item + "/map"))
        {
            string sample =base_name(entry.path());
            string origin = entry.path();
            if (fileExists(origin + "/" + sample + ".kir-mapper.log"))
            {
                v_command = "cp -r -u " + item + "/map/" + sample + " " + v_output + "/map/";
                GetStdoutFromCommand(v_command);
            }
                
        }
    }
  
  
    vector <string> data;
    v_message = "Copying plots ...";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    for (auto item : outputs)
    {
        string name = base_name(item);
        v_command = "mkdir " + v_output + "/ncopy/" + name;
        GetStdoutFromCommand(v_command);
        
        v_command = "cp -r " + item + "/ncopy/ " + v_output + "/ncopy/" + name;
        GetStdoutFromCommand(v_command);
        
        string cn = item + "/ncopy/copy_numbers.txt";
        
        ifstream file(cn.c_str());
        for( std::string line; getline( file, line );)
        {
            if (line.find("COPY_NUMBER") != std::string::npos) {continue;}
            data.push_back(line);
        }
        file.close();
    }
    
    string cn = v_output + "/ncopy/copy_numbers.txt";
    ofstream CP;
    CP.open (cn.c_str());
    CP << "SAMPLE\tGENE\tCOPY_NUMBER" << endl;
    for (auto item : data)
    {
        CP << item << endl;
    }
    CP.close();
             
    v_message = "Computation done.";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    
}
