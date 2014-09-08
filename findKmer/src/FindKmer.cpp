/*
 * Copyright (c) 2014 Kalen A. Brown, August C. Thies, Gavin Conant, Xiang Wang,
 * Michela Becchi and University of Missouri in Columbia.
 * All rights reserved
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *    1. Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *    2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *    3. The name of the author or the University may not be used
 *       to endorse or promote products derived from this source code
 *       without specific prior written permission.
 *    4. Conditions of any other entities that contributed to this are also
 *       met. If a copyright notice is present from another entity, it must
 *       be maintained in redistributions of the source code.
 *    5. You notify the author and give your intentions.
 *       Notification can be given to kab8c8 at mail dot missouri dot edu
 *
 * THIS INTELLECTUAL PROPERTY (WHICH MAY INCLUDE BUT IS NOT LIMITED TO SOFTWARE,
 * FIRMWARE, VHDL, etc) IS PROVIDED BY  THE AUTHOR AND THE UNIVERSITY
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR THE UNIVERSITY
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS INTELLECTUAL PROPERTY, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * */

//============================================================================
// Name        : findKmer.cpp
// Author      : Kalen Brown and Gus Thies
// Version     : This version is currently in transition to object orientated and iterative code
// Copyright   : Do not copy
// Description : This code uses k which is a ones count for the length of a sequence to be found in DNA.
// Expects     : The parse file is a DNA file of any file extension to be formatted in this way:
//                 All valid bases are in upper case (A, C, T or G) not lower case.
//                 Sequences are broken by any character other than '\n', A, C, T or G
//                 '>' designates the beginning of an ID and newlines ('\n') designates the end of an ID
//                 Newlines ('\n') are otherwise ignored completely
//                 No kmer sequence will contain the letter N at any position, thus it breaks a sequence
// Bugs        : Known bugs include memory leaks.
// TODO        : Passing by value actually is worse than pass by reference for single elements of some types, weed out that case (esp 64 bit sys = 64 bit pointer)
// TODO        : A choice could be made in development to either minimize storage of tree in memory by having different nodes for branch and leaf,
// TODO        : Changing the findKmer function to take in a single letter and a size k and have it only increment the counter when it reaches depth == its k.
// TODO        : This would be easily done if findKmer was a class that we could instantiate n classes and they would retain their values.
// Fixed bugs  : Verified that the tree is actually creating the correct number of nodes. # of nodes != 4^k; see estimate ram usage.
// Compile     : g++  -o "findKmer" [-O3 seems to work OK but unknown benefit]
//============================================================================
//API creation header files.
#include "FindKmer.h"
#include "CmdLineParser.h" //header file for object oriented command line parser. use "" for local files.
#include "MemMgt.h" //header file for memory management tracking total allocated/deallocated bytes.

FindKmer::FindKmer() {
	//can wait and only initialize and store pointers to new objects when they are needed...
	initalizeConfiguration(&config_);
}

/*
 * parses command line arguments and displays configuration to the user.
 */
FindKmer::FindKmer(int argc, char** argv) {
	CmdLineParser *tempParser = new CmdLineParser();
	tempParser->usage();
	tempParser->parse_arguments(argc,argv);
	tempParser->print_conf(argc);
	tempParser->estimate_RAM_usage();
	tempParser->getConfiguration(&config_);
	delete tempParser;
}


FindKmer::~FindKmer() {
	//free the file name strings.
	if (config_.out_file_name != NULL) {

		fprintf(stdout,
				"Your file can be found in the current directory as: \n    %s\n",
				config_.out_file_name);
		delete (config_.out_file_name);
	}
	if (config_.sequence_file_name != NULL) {
		delete (config_.sequence_file_name);
	}
	//free(sequence_file);
	//free(out_file);

	//If we try to use fclose on a null pointer then the program crashes.
	if (config_.out_file_pointer) {
		if (fclose(config_.out_file_pointer) == EOF) {
			fprintf(stderr,
					"Out file close error! This is not expected and might mean the data was not written to the file properly before the close.\n");
		}
	}

	if (fclose(config_.sequence_file_pointer) == EOF) {
		fprintf(stderr,
				"Sequence file close error! This is likely ok though.\n");
	}
}

const configuration& FindKmer::getConfig() const {
	return config_;
}

/*
 * TODO Requires Lots of error checking! leaving it for now.
 */
void FindKmer::setConfig(const configuration& config) {
	config_ = config;
}
