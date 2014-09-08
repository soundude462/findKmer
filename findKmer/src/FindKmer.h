/*
 * ============================================================================
 * FINDKMER.h of findKmer project.
 *
 *  Created on: Aug 27, 2014
 *      Author: Kalen Brown 
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

#ifndef FINDKMER_H_
#define FINDKMER_H_


//previously included headers taken from original findKmer.cpp file.
using namespace std;
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h> //for strcmp(string1,string2) string comparison returns a 0 if they are the same.
#include <stdlib.h> //malloc is in this.
#include <fstream> // basic file operations




/*
 * Below are some defaults you can setup at compile time.
 * Any combination of command line arguments can override these.
 * Default sequence file names must be a string literal with file extension.
 * These macros may bypass the error checking scheme!
 * The default suppress output value set to 1 OR GREATER will bypass DNA ID printing and getchar() breaks
 * The default z threshold enable set to 1 OR GREATER causes outfile to only contain sequences with z score above z threshold.
 */
//#define DEFAULT_SEQUENCE_FILE_NAME "homo_sapiensupstream.fas"
//#define DEFAULT_SEQUENCE_FILE_NAME "Full_homo_sapiens.fa"
#define DEFAULT_SEQUENCE_FILE_NAME "test.txt"
//#define DEFAULT_SEQUENCE_FILE_NAME "shortend_test_Homo_sapiens_1_and_2.fa"

#define DEFAULT_K_VALUE 7
#define DEFAULT_SUPPRESS_OUTPUT_VALUE 0
#define DEFAULT_Z_THRESHOLD_ENABLE 0
#define DEFAULT_Z_THRESHOLD 1000

/*
 * This will appear as the first line of the output file.
 * We put the shannon entropy after the sequence because it is common to the sequence on both the genome and the upstream.
 * The frequency is not constant but the hope was that it would make combining files easier by grouping common things and uncommon things.
 */
#define OUT_FILE_COLUMN_HEADERS "Sequence, Shannon Entropy h, Shannon Entropy H, Frequency, Z score"

/*
 * Commenting out the last x will NOT make the tree in memory and will NOT create a valid histogram.
 * It is only to save memory and count base occurrences.
 */
#define COUNT_BASES_ONLY(x) x

//debugging
#define DEBUG(x) x
#define DEBUG_TREE_CREATE(x) //x
#define DEBUG_HISTO_AND_FREE_RECURSIVE(x) //x
#define DEBUG_SHIFT_AND_INSERT(x) //x
#define DEBUG_STATISTICS(x) //x
#define DEBUG_FREE(x) //x

//Holds the statistics for a base
struct statistics_t {
	unsigned int Count; //number of time this base was encountered in the entire file.
	long double Probability; //probability that this base will be encountered out of all bases in the file.
};

struct configuration {
	const char *sequence_file_name; //holds the string representation of the file name.
	FILE *sequence_file_pointer; //holds the FILE pointer to the file itself
	const char *out_file_name;		//holds the string representation of the file name.
	FILE *out_file_pointer;  //holds the FILE pointer to the file itself
	int k; //holds the length of k for the size of the sequence to be recorded.
	int suppressOutputEnable; //Suppress identifier printing and getchar(); breaks.
	long double zThreshold; //holds the minimum Z score value to print to outfile
	int zThresholdEnable; //The z threshold enable set to 1 OR GREATER causes outfile to only contain sequences with z score above z threshold.
};

inline void initalizeConfiguration(configuration * config){
	config->sequence_file_name = NULL;
	config->sequence_file_pointer = NULL;
	config->out_file_name = NULL;
	config->out_file_pointer = NULL;
	config->k = 0;
	config->suppressOutputEnable = -1;
	config->zThresholdEnable = -1;
	config->zThreshold = -1;
}
/*
 * returns EXIT_SUCCESS on success.
 * Checks for NULL and checks for k range of 1 to 20.
 * File pointers and names != NULL.
 */
inline int checkConfiguration(configuration * config){
	if(config->sequence_file_name != NULL &&
	config->sequence_file_pointer != NULL &&
	config->out_file_name != NULL &&
	config->out_file_pointer != NULL &&
	config->k > 0 && config->k <= 20){
		return EXIT_SUCCESS;
	}else{
		return EXIT_FAILURE;
	}
}

/* Data structure for a tree.
 * http://msdn.microsoft.com/en-us/library/s3f49ktz.aspx
 * holds the ranges of each data type.
 * The number of bases in the human genome is
 * 2,858,658,142 bases which is less than a signed int.
 * Unsigned int has a range of 4,294,967,295
 * Unsigned char is used for the base even though the
 * program treats it as an int elsewhere to save space.
 */
struct node_t {
	unsigned char base; //the letter that this node holds. if depth = k then it is the last letter in the kmer sequence.
	struct node_t *nextNodePtr[4]; //pointer to next node
	unsigned int frequency; //number of times that this "sequence" was encountered in the whole file.
};

//Not yet used: Extra functions taken from Dr. Becchi's stdinc.h file.
//Inline functions are like, but safer than, macros.
/* max and min */
inline int max(int x, int y) { return x > y ? x : y; }
inline double max(double x, double y) { return x > y ? x : y; }
inline int min(int x, int y) { return x < y ? x : y; }
inline double min(double x, double y) { return x < y ? x : y; }

/* warnings and errors */
inline void warning(const char* p) { fprintf(stderr,"Warning:%s \n",p); }
inline void fatal(const char* string) {fprintf(stderr,"Fatal:%s\n",string); exit(EXIT_FAILURE); }

class FindKmer {
public:
	FindKmer();
	FindKmer(int argc, char **argv);
	virtual ~FindKmer();
	const configuration& getConfig() const;
	void setConfig(const configuration& config);

private:
	configuration config_;
};


#endif /* FINDKMER_H_ */
