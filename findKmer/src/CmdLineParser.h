/*
 * ============================================================================
 * cmdLineParser.h of findKmer project.
 *
 *  Created on: Aug 26, 2014
 *      Author: Kalen A. Brown
 *
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
 * ============================================================================
 */

#ifndef CMDLINEPARSER_H_
#define CMDLINEPARSER_H_
#include "FindKmer.h"
#include "MemMgt.h"

/*
 * This class is required to build the configuration struct from command line arguments.
 * FindKmer will not work without this configuration unless you do it manually.
 */
class CmdLineParser {

public:
	//instantiates the parser
	CmdLineParser();

	//parser de-allocator
	~CmdLineParser();

	//prints command line arguments that are available.
	void usage();

	int parse_arguments(int argc, char **argv);

	//print the current configuration
	void print_conf(int argc);

	//Check to see if the file can be opened in the given mode.
	void check_file(const char* filename, const char* mode);

	//Estimates RAM usage given the current configuration.
	unsigned long int estimate_RAM_usage();

	//Simple Getters
	int getK() const;
	const FILE* getOutFilePointer() const;
	const FILE* getSequenceFilePointer() const;
	int getSuppressOutputEnable() const;
	long double getThreshold() const;
	int getThresholdEnable() const;
	unsigned long int getTotalAllocatedBytes() const;
	const char* getSequenceFileName() const;
	char* getOutFileName() const;
	int getConfiguration(configuration *config);
	bool isValidConfiguration() const;
	bool isConfigurationCopied() const;

private:
	const char *sequence_file_name_; //holds the string representation of the file name.
	FILE *sequence_file_pointer_; //holds the FILE pointer to the file itself
	char *out_file_name_;		//holds the string representation of the file name.
	FILE *out_file_pointer_;  //holds the FILE pointer to the file itself
	int k_; //holds the length of k for the size of the sequence to be recorded.
	int suppressOutputEnable_; //Suppress identifier printing and getchar(); breaks.
	long double zThreshold_; //holds the minimum Z score value to print to outfile
	int zThresholdEnable_; //The z threshold enable set to 1 OR GREATER causes outfile to only contain sequences with z score above z threshold.
	memMgt * memMgr_;
	bool validConfiguration_;
	bool configurationCopied_;

	void set_default_conf_();
};

#endif /* CMDLINEPARSER_H_ */
