//
//  exome_bias.hpp
//  kir-mapper
//
//  Created by Erick Castelli on Jan/10/24.
//  Copyright © 2024 GeMBio.Unesp. All rights reserved.
//

#ifndef exome_bias_hpp
#define exome_bias_hpp

#include <stdio.h>
void main_exome_bias();
string calculate_capture_bias (std::vector<string> files, string gene, string v_out, string fileout, int number_of_tests);
string calculate_depth (std::vector<string> files);
#endif /* exome_bias_hpp */
