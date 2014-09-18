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
 *
 * */

//============================================================================
// Name        : findKmer.cpp
// Author      : Kalen Brown and Gus Thies
// Version     :
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
 * The default z threshold enable set to 1 OR GREATER causes outfile to only contain sequences with abs(z) score above z threshold.
 */
//#define DEFAULT_SEQUENCE_FILE_NAME "homo_sapiensupstream.fas"
//#define DEFAULT_SEQUENCE_FILE_NAME "Full_homo_sapiens.fa"
#define DEFAULT_SEQUENCE_FILE_NAME "test.txt"

#define DEFAULT_K_VALUE 7
#define OUT_FILE_COLUMN_HEADERS "Sequence, Shannon Entropy h, Shannon Entropy H, Frequency, Z score"
#define DEFAULT_SUPPRESS_OUTPUT_VALUE 0
#define DEFAULT_Z_THRESHOLD_ENABLE 0
#define DEFAULT_Z_THRESHOLD 1000

//debugging
#define DEBUG(x) //x
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

/* structure definition for configuration of file names, pointers, and length of k.
 * For enables: 0 == false, > 1 is true, < 1 means none supplied from user.
 */
static struct conf {
	const char *sequence_file; //holds the string representation of the file name.
	FILE *sequence_file_pointer; //holds the FILE pointer to the file itself
	char *out_file;		//holds the string representation of the file name.
	FILE *out_file_pointer;  //holds the FILE pointer to the file itself
	int k; //holds the length of k for the size of the sequence to be recorded.
	int suppressOutputEnable; //Suppress identifier printing and getchar(); breaks.
	long double zThreshold; //holds the minimum Z score value to print to outfile
	int zThresholdEnable; //The z threshold enable set to 1 OR GREATER causes outfile to only contain sequences with z score above z threshold.
} config; /* Config is a GLOBAL VARIABLE for configuration of file names, pointers, and length of k.*/

//Global variable that needs to be localized.
unsigned long long int nodeCounter = 0; //number of nodes created in memory.

extern int recurse_factorial(int i) {
	if (i > 1)
		return (i * recurse_factorial(i - 1));
	else
		return (1);
}
extern long double float_factorial(int i) {
	int j;
	long double val = 1.0;

	if (i > 1) {
		for (j = 2; j <= i; j++)
			val = val * (long double) j;
		return (val);
	} else
		return (1.0);
}
int n_choose_k(int n, int k)
//Note that this function assumes that n choose k can be represented as the 32-bit int
		{
	int i, retval = 1, large_denom, small_denom;

	if (k > (n - k)) {
		large_denom = k;
		small_denom = n - k;
	} else {
		large_denom = n - k;
		small_denom = k;
	}
	for (i = n; i > large_denom; i--)
		retval *= i;

	retval = retval / recurse_factorial(small_denom);
	return (retval);
}
long double float_n_choose_k(int n, unsigned int k) {
	int i;
	long double fretval = 1.0, large_denom, small_denom;

	if (k > (n - k)) {
		large_denom = k;
		small_denom = n - k;
	} else {
		large_denom = n - k;
		small_denom = k;
	}
	for (i = n; i > large_denom; i--) {
		fretval = fretval * (long double) i;
	}

	fretval = fretval / float_factorial(small_denom);
	return (fretval);

}
/*
 * This is a test to see if we can do the normal approximation test or not.
 * The mean must be greater than or equal to five
 * N*Q must be greater than or equal to five
 */
bool normal_approx_check(unsigned long long n, long double p, long double q) {
	bool pass = true;

	if ((n * p >= 5) && (n * q >= 5)) {
		pass = true;
	} else {
		pass = false;
	}
	return pass;
}
void *allocate_array(int size, size_t element_size) {
	void *mem = malloc(size * element_size);
	if (!mem) {
		fprintf(stderr, "allocate_array():: memory allocation failed\n");
		exit(EXIT_FAILURE);
	}
	return mem;
}
void *reallocate_array(void *array, int size, size_t element_size) {
	void *new_array = realloc(array, element_size * size);
	if (!new_array) {
		fprintf(stderr, "reallocate_array():: memory reallocation failed\n");
		exit(EXIT_FAILURE);
	}

	return new_array;
}
void deallocate_array(void** array) {
	free(*array);
	*array = NULL;
}
//for testing purposes.
void random_array(int sizeOfArray) {
	//Initializing array to test code. this will come from the pre processed line
	int *array = (int*) allocate_array(sizeOfArray, sizeof(int));
	int* currentArrayPosition = array;
	srand(time(NULL)); //seeding the random function.
	for (int i = 0; i < sizeOfArray; i++) {
		*currentArrayPosition = (rand() % 4);
		currentArrayPosition++;
	}
	deallocate_array((void**) &array);
}
/* check that the given file can be read/written */
void check_file(const char *filename, const char *mode) {
	FILE *file = fopen(filename, mode);
	if (file == NULL) {
		fprintf(stderr,
				"Unable to open file %s in %s mode\nFile MUST be in current directory.\n",
				filename, mode);
		exit(EXIT_FAILURE);
	} else
		fclose(file);
}
/* initialize the configuration
 * Set to null or an invalid value to determine default or user defined.
 */
void init_conf() {
	config.sequence_file = NULL;
	config.sequence_file_pointer = NULL;
	config.out_file = NULL;
	config.out_file_pointer = NULL;
	config.k = 0;
	config.suppressOutputEnable = -1;
	config.zThresholdEnable = -1;
	config.zThreshold = -1;
}
/* This function fills in any gaps in the configuration file.*/
void set_default_conf() {

	if (!config.sequence_file) {
		config.sequence_file = DEFAULT_SEQUENCE_FILE_NAME;
	}

	if (!config.k) {
		config.k = DEFAULT_K_VALUE;
	}

	//double check default and user defined K value.
	if (!config.k) {
		fprintf(stdout, "k must be greater than zero. Ending program.\n\n");
		exit(EXIT_FAILURE);
	}
	//if lessthan zero, then the user did not specify.
	if (config.suppressOutputEnable < 0) {
		config.suppressOutputEnable = DEFAULT_SUPPRESS_OUTPUT_VALUE;
	}

	if (config.zThresholdEnable < 0) {
		config.zThresholdEnable = DEFAULT_Z_THRESHOLD_ENABLE;
		config.zThreshold = DEFAULT_Z_THRESHOLD;
	}

	if (!config.out_file) {
		const char* nameOfFile = "mer_Historam_Of_";
		const char* outFileExension = ".csv";
		const char* zScoreFiltered = "zScoreFiltered";

		if (config.zThresholdEnable == 0) {
			config.out_file = (char*) allocate_array(
					strlen("999") + strlen(nameOfFile)
							+ strlen(config.sequence_file)
							+ strlen(outFileExension), sizeof(char));
			sprintf(config.out_file, "%d%s%s%s", config.k, nameOfFile,
					config.sequence_file, outFileExension);
		} else {
			config.out_file = (char*) allocate_array(
					strlen("999") + strlen(nameOfFile)
							+ strlen(config.sequence_file)
							+ strlen(outFileExension) + strlen(zScoreFiltered),
					sizeof(char));
			sprintf(config.out_file, "%d%s%s%s%s", config.k, nameOfFile,
					config.sequence_file, zScoreFiltered, outFileExension);
		}
	}

}
/* print the configuration */
void print_conf(int argc) {
	fprintf(stdout, "\nATTEMPTING CONFIGURATION: \n");
	set_default_conf();
	if (config.sequence_file)
		fprintf(stdout, "- sequence_file file: %s\n", config.sequence_file);
	if (config.out_file)
		fprintf(stdout, "- export file: %s\n", config.out_file);
	if (config.k)
		fprintf(stdout, "- k size: %d\n", config.k);

	fprintf(stdout, "- %s\n",
			config.suppressOutputEnable > 0 ?
					"Suppressing file read output and breaks." :
					"Showing DNA Sequence identifier and allowing breaks.");

	fprintf(stdout, "- Z score filtering is %s",
			config.zThresholdEnable ? "enabled" : "disabled");

	if (config.zThresholdEnable > 0) {
		fprintf(stdout, "\n    with threshold of %LG", config.zThreshold);
	}
	fprintf(stdout, ".\n");

	//if suppressOutputEnable is false and no command line arguments have been given:
	if (config.suppressOutputEnable == 0 && argc < 2) {
		fprintf(stdout, "Press enter to proceed with this configuration.");
		getchar();
	}

	/* Double check configuration */
	if (config.k < 0 || config.k > 20) {
		fprintf(stderr,
				"%d is not a valid value for k. Please select a number greater than zero\n",
				config.k);
		exit(EXIT_FAILURE);
	}

	if ((config.sequence_file_pointer = fopen(config.sequence_file, "r"))
			!= NULL) {
		//fprintf(stdout, "Sequence file opened properly\n");
	} else {
		fprintf(stderr, "Sequence file failed to open\n\n");
		exit(EXIT_FAILURE);
	}

	if ((config.out_file_pointer = fopen(config.out_file, "w")) != NULL) {
		//fprintf(stdout, "Out file opened properly\n");
		fprintf(config.out_file_pointer, OUT_FILE_COLUMN_HEADERS);
	} else {
		fprintf(stderr,
				"Out file failed to open\nFile MUST be in current directory.\n");
		exit(EXIT_FAILURE);
	}

	fprintf(stdout, "Sequence file and out file opened properly\n");

	fprintf(stdout, "\n");
}
/* usage */
static void usage() {
	fprintf(stdout, "\n");
	fprintf(stdout, "Usage: findKmer [options]\n");
	fprintf(stdout, "             [--parse|-p <sequence_file.txt>] \n"
			"               File with DNA sequence data.\n"
			"               File must be in current directory.\n"
			"               Parser follows .fas and .fa formats\n"
			"                Default is %s.\n\n",
	DEFAULT_SEQUENCE_FILE_NAME);
	fprintf(stdout, "             [--export|-e  <out_file.csv>] \n"
			"               File to output histogram data to.\n"
			"                Default output file name is dynamic.\n\n");
	fprintf(stdout, "             [--ksize|-k  <k>] \n"
			"               Size of sequence for histogram.\n"
			"                Default is %d.\n\n",
	DEFAULT_K_VALUE);
	fprintf(stdout, "             [--quiet|-q  < 0 for FALSE | 1 for TRUE >] \n"
			"               Suppress file read output and breaks.\n"
			"                Default is %s.\n\n",
	DEFAULT_SUPPRESS_OUTPUT_VALUE ? "true" : "false");

	long double tempzThreshold = DEFAULT_Z_THRESHOLD;
	fprintf(stdout, "             [--zthreshold|-z  < Threshold_for_Z >] \n"
			"               Suppress sequences with Z scores < threshold.\n"
			"                Default is %s with a value of %LG.\n\n",
	DEFAULT_Z_THRESHOLD_ENABLE ? "enabled" : "disabled", tempzThreshold);
	fprintf(stdout, "\n");
}
int parse_arguments(int argc, char **argv) {
	int i = 1;
	if (argc < 2) {
		return 1;

	} else {

		while (i < argc) {
			if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
				exit(1);
			} else if (strcmp(argv[i], "-e") == 0
					|| strcmp(argv[i], "--export") == 0) {
				i++;
				if (i == argc) {
					fprintf(stderr, "Export file name missing.\n");
					return 0;
				}

				if (argv[i] != NULL) {
					check_file(argv[i], "w");
				}
				config.out_file = argv[i];
			} else if (strcmp(argv[i], "-p") == 0
					|| strcmp(argv[i], "--parse") == 0) {
				i++;
				if (i == argc) {
					fprintf(stderr, "Sequence data file name missing.\n");
					return 0;
				}

				if (argv[i] != NULL) {
					check_file(argv[i], "r");
				}

				config.sequence_file = argv[i];
			} else if (strcmp(argv[i], "-k") == 0
					|| strcmp(argv[i], "--ksize") == 0) {
				i++;
				if (i == argc) {
					fprintf(stderr, "Number for size of k is missing.\n");
					return 0;
				} else {

					int k = atoi(argv[i]);
					if (k < 0 || k > 20) {
						fprintf(stderr,
								"%d is not a valid value for k.\nPlease select a number greater than zero and less than 21\n",
								k);
						exit(EXIT_FAILURE);
					}
					config.k = k;
				}
			} else if (strcmp(argv[i], "-q") == 0
					|| strcmp(argv[i], "--quiet") == 0) {
				i++;
				if (i == argc) {
					fprintf(stderr,
							"True/false value for quiet option is missing.\nUsage is \"-q 1\" for suppression OR \"-q 0\" for expansion\n");
					exit(EXIT_FAILURE);
				} else {

					int suppressOutputEnableOption = atoi(argv[i]);
					if (suppressOutputEnableOption == 1
							|| suppressOutputEnableOption == 0) {
						config.suppressOutputEnable =
								suppressOutputEnableOption;
					} else {
						fprintf(stderr,
								"%d is not a valid value for suppress Output Enable Option.\nPlease select either 0 for FALSE or a 1 for TRUE",
								suppressOutputEnableOption);
						exit(EXIT_FAILURE);
					}
				}
			} else if (strcmp(argv[i], "-z") == 0
					|| strcmp(argv[i], "--zthreshold") == 0) {
				i++;
				if (i == argc) {
					fprintf(stderr,
							"Z threshold number is missing\nUsage is \"-z 1000\".\n");
					exit(EXIT_FAILURE);
				} else {
					config.zThresholdEnable = 1;
					config.zThreshold = atoi(argv[i]);
				}
			} else {
				fprintf(stderr, "Ignoring invalid option %s\n", argv[i]);
				if (config.suppressOutputEnable == 0) {
					fprintf(stderr, "Press enter to continue.\n");
					getchar();
				}
			}
			i++;
		}
	}

	return 1;
}
void statistics(unsigned long long * const baseCounter,
		statistics_t * const baseStatistics,
		unsigned long long * const TotalNumSequencesN,
		unsigned long int * const maxNumberOfNodes) {
	const char* nameOfFile = "mer_Base_Stats_Of_";
	const char* outFileExension = ".txt";

	char* stats_out_file_name = (char*) allocate_array(
			strlen("999") + strlen(nameOfFile) + strlen(config.sequence_file)
					+ strlen(outFileExension), sizeof(char));

	sprintf(stats_out_file_name, "%d%s%s%s", config.k, nameOfFile,
			config.sequence_file, outFileExension);

	FILE * stats_out_file_pointer = NULL;
	if ((stats_out_file_pointer = fopen(stats_out_file_name, "w")) == NULL) {
		fprintf(stderr,
				"Out file failed to open\nFile MUST be in current directory.\n");
		exit(EXIT_FAILURE);
	}

	fprintf(stdout,
			"Statistics of occurrences and probability of A, C, G and T respectively: \n");
	for (int i = 0; i < 4; i++) {
		DEBUG_STATISTICS(
				fprintf(stdout, "%u count / %llu baseCounter\n",baseStatistics[i].Count , *baseCounter));

		fprintf(stdout, "%u", baseStatistics[i].Count);

		baseStatistics[i].Probability = (double) baseStatistics[i].Count
				/ *baseCounter;
		if (baseStatistics[i].Probability == 0.0) {
			fprintf(stdout, "Division overflow detected in statistics.\n");
			exit(EXIT_FAILURE);
		}

		fprintf(stdout, ", %Lf\n", baseStatistics[i].Probability);
		fprintf(stats_out_file_pointer, "%Lf\n", baseStatistics[i].Probability);
	}

	fprintf(stdout, "Found %llu valid bases total INSIDE sequences >= k.\n",
			*baseCounter);

	DEBUG(cout << (*TotalNumSequencesN) << " sequences of length k were found."<<endl<<" This is NOT the number of combinations found."
			<< endl

			cout << nodeCounter << " Nodes created " << endl;

			cout << (*maxNumberOfNodes) << " Max possible Nodes expected " << endl;);

	fprintf(stdout, "%0.0f%% tree density.\n",
			((double) (nodeCounter) / (double) (*maxNumberOfNodes)) * 100);

	if (nodeCounter == (*maxNumberOfNodes)) {
		fprintf(stats_out_file_pointer,
				"All possible %dmers combinations were found.\n", config.k);
		fprintf(stdout, "All possible kmer combinations were found.\n");
	} else if (nodeCounter > (*maxNumberOfNodes)) {
		fprintf(stderr,
				"Error! too many nodes were created!\nThere may be a corruption of data!\n");
		fprintf(stats_out_file_pointer,
				"too many nodes were created when looking for %dmers.\n",
				config.k);
	} else {
		fprintf(stdout, "FYI we did not find all possible combinations.\n");
		fprintf(stats_out_file_pointer,
				"did not find all possible %dmers combinations.\n", config.k);
	};

	free(stats_out_file_name);
	fclose(stats_out_file_pointer);

	return;

}
// Gus' Functions
int base2int(char base) {
	int integer = -1;
	if (base == 'A') {
		integer = 0;
	} else if (base == 'T') {
		integer = 3;
	} else if (base == 'C') {
		integer = 1;
	} else if (base == 'G') {
		integer = 2;
	} else if (base == 'N') {
		integer = -2;
	} else if (base == EOF) {
		DEBUG(fprintf(stderr,"End of file \n"));
	} else {
		fprintf(stderr,
				"Unknown character %c processed! File may be corrupted.\n",
				base);
	}

	//close the function
	return integer;
}
char int2base(int integer) {
	char base;
	if (integer == 0) {
		base = 'A';
	} else if (integer == 3) {
		base = 'T';
	} else if (integer == 1) {
		base = 'C';
	} else if (integer == 2) {
		base = 'G';
	} else {
		base = 'e';
	}

	//close the function
	return base;
}
/*
 * Creates a tree node.
 * Brings in the base of the node to create
 *
 */
node_t* node_create(int base) {
	node_t* node = (node_t*) malloc(sizeof(node_t));
	node->base = base;
	node->frequency = 1;
	node->nextNodePtr[0] = NULL;
	node->nextNodePtr[1] = NULL;
	node->nextNodePtr[2] = NULL;
	node->nextNodePtr[3] = NULL;
	nodeCounter++;
	return node;
}
/*
 * Adds a branch if one doesn't already exist.
 * Increments the counter for the node we are going to step into.
 * Returns the pointer to the next node.
 */
node_t* node_branch_enter_and_create(node_t* node, int base) {
	DEBUG_TREE_CREATE(
			fprintf(stdout, "..node_branch_enter_and_create for base %c\n", int2base(base)));
	//if node doe
	if (node->nextNodePtr[base] == NULL) {

		DEBUG_TREE_CREATE(fprintf(stdout, "***Creating Node.\n"));

		node->nextNodePtr[base] = node_create(base);

	} else {

		node->nextNodePtr[base]->frequency++;

		if (node->nextNodePtr[base]->frequency == 0) {
			fprintf(stderr,
					"\n\n!!! COUNTER ROLLOVER DETECTED! \nIncrease the number of bits used for the counter variable if you have the source code, else use a smaller sequence file.\n\n");
			fprintf(stdout,
					"\n\n!!! COUNTER ROLLOVER DETECTED! \nIncrease the number of bits used for the counter variable if you have the source code, else use a smaller sequence file.\n\n");
			exit(EXIT_FAILURE);

		}DEBUG_TREE_CREATE(
				fprintf(stdout, "+++Incrementing counter to %d.\n", node->nextNodePtr[base]->frequency));

	}DEBUG_TREE_CREATE(fprintf(stdout, "..returning next base pointer.\n"));
	return node->nextNodePtr[base];
}
/*
 * Brings in a pointer to head of the tree, an integer array and the size k of the array
 * Checks to see if the array exists, and if the tree already exists.
 * Breaks if array does not exist, but creates the head node if tree does not exist.
 * Traverses the array and creates the tree based on what it finds.
 *
 *
 */
node_t* tree_create(node_t* head, int* const array, int k,
		statistics_t *baseStatistics) {
	if (array != NULL) {
		if (head == NULL) {
			DEBUG_TREE_CREATE(
					fprintf(stdout, "-Creating the head of the tree!\n"));
			head = node_create('H');
		}
		node_t *currentNode = head;

		/*
		 * Traverse the given integer array that is of k size.
		 * For each base within the integer array, create a node and enter the created node.
		 * Also record how often each base (A, C, G or T).
		 */
		for (int i = 0; i < k; i++) {
			DEBUG_TREE_CREATE(
					fprintf(stdout, "-Moving into a branch on depth %d\n", i););
			currentNode = node_branch_enter_and_create(currentNode, array[i]);

		}
	} else {
		fprintf(stderr,
				"kmer array was empty, this is embarrassing. Quitting the program. \n");
		exit(1);
	}
	return head;
}
/*
 * Histogram can be recursive for low numbers of K.
 * If K becomes too high then we may run out of stack/heap memory.
 * It will traverse the tree starting at the head and work its way to the depth provided by k.
 * When it reaches the depth k, it will print out the sequence it saw from the root to this point.
 * It will then work its way back
 * The implementation of freeing the tree is not implemented here.
 */
void histo_recursive(node_t* const head, int * const array, const int depth,
		const int k, unsigned long long * const baseCounter,
		statistics_t * const baseStatistics,
		unsigned long long * const TotalNumSequencesN) {
	DEBUG_HISTO_AND_FREE_RECURSIVE(
			fprintf(stdout, "histo&free @ depth %d of %d has %d\n",depth,k,head != NULL?head->base:-1 ));

	if (head == NULL) {
		DEBUG_HISTO_AND_FREE_RECURSIVE(
				fprintf(stdout, "histo_and_free::Head == NULL. Leaf found.\n"));
	} else {
		if (head->base == 'H') {
			DEBUG_HISTO_AND_FREE_RECURSIVE(fprintf(stdout, "Head found.\n\n"));
		} else {
			array[depth - 1] = head->base;
			DEBUG_HISTO_AND_FREE_RECURSIVE(
					for(int z =0; z < depth; z++) {fprintf(stdout, " ");}fprintf(stdout, "writing %c to array\n", int2base(head->base)));
		}

		//check each branch in this node even if we "think" its the leaf, to be safe
		for (int i = 0; i < 4; i++) {
			DEBUG_HISTO_AND_FREE_RECURSIVE(
					fprintf(stdout, "histo&free @ depth %d of %d has %d checking branch %d\n",depth,k,head != NULL?head->base:-1,i ));
			histo_recursive(head->nextNodePtr[i], array, depth + 1, k,
					baseCounter, baseStatistics, TotalNumSequencesN);
		}

		//once we have exhausted all branches, we check for depth of k.
		if (depth == (k)) {
			statistics_t kmerBaseStatistics[4] = { 0 }; //This will hold data that is only for this single Kmer and not for the entire file.

			DEBUG_STATISTICS(
					for (int i = 0; i < 4; i++) {
						cout << kmerBaseStatistics[i].Count << " = count and "
						<< kmerBaseStatistics[i].Probability
						<< " = probability initially" << endl
						;
					});

			/*
			 * count the number of times each base occurs. GATTACA,
			 * kmerBaseStatistics[base2int('A')].Count = 3,kmerBaseStatistics[base2int('C')].Count = 1,
			 * kmerBaseStatistics[base2int('G')].Count = 1, kmerBaseStatistics[base2int('T')].Count = 2,
			 */
			DEBUG_STATISTICS(cout << "pre traversing kmer" << endl);
			for (int location = 0; location < k; location++) {
				DEBUG_STATISTICS(
						cout << "location == " << location << endl; cout << "array[location] == "
						<< array[location] << endl; cout << "kmerBaseStatistics[array[location]].Count == "
						<< kmerBaseStatistics[array[location]].Count << endl;);

				kmerBaseStatistics[array[location]].Count++; //increment the counter for this letter

				DEBUG(fprintf(stdout, "%c", int2base(array[location])));
			}DEBUG(fprintf(stdout, ", %d\n", head->frequency));

			/*
			 * Calculate Shannon Entropy to determine if a sequence contains information. it could be estimated
			 * # of different bases in the sequence, # of bits estimated
			 * 1, 0
			 * 2, 1
			 * 3, 2
			 * 4, 2
			 * Then multiply by the number of letters in the sequence.
			 * AAAAAAAAA has zero bits of information.
			 * GATTACA has 14 bits of information to encode the entire sequence and still be able to decode it.
			 * H(X) = (over x) Σ P(x) * log2(1/P(x)) in bits
			 * Where P(x) is the probability of the current letter occurring in the current KMER sequence.
			 * TODO create a dynamic shannon entropy limit filter. H > 0 for sure but H > k would be ok, is H > k*2 ok? or k*4
			 */
			//calculate the probability
			for (int i = 0; i < 4; i++) {
				kmerBaseStatistics[i].Probability =
						(double) kmerBaseStatistics[i].Count
								/ (double) config.k;
				DEBUG_STATISTICS(
						cout << kmerBaseStatistics[i].Probability << " = "
						<< kmerBaseStatistics[i].Count << " / "
						<< config.k << endl);
			}

			DEBUG_STATISTICS(
					for (int i = 0; i < 4; i++) {
						cout << kmerBaseStatistics[i].Count << " = count and "
						<< kmerBaseStatistics[i].Probability
						<< " = probability calculate the probability"
						<< endl
						;
					}

			);

			//calculate the number of bits to encode a single symbol
			long double h = 0;
			for (int i = 0; i < 4; i++) {
				if (kmerBaseStatistics[i].Probability != 0)
					h +=
							(double) kmerBaseStatistics[i].Probability
									* log2(
											1
													/ (double) kmerBaseStatistics[i].Probability);
			}

			DEBUG_STATISTICS(cout << h << " = h" << endl
					; );

			//calculate the number of bits to encode the entire sequence.
			long double H = h * config.k;

			DEBUG_STATISTICS(cout << H << " = H" << endl
					; );

			/*
			 * Find the likely hood that this base occurred this many times randomly
			 * based on the proportion of times that we found it in the file.
			 * I.e. we found A 50% of the time, and we found it 3 times,
			 * Run this on all possible bases A,C,G, and T
			 * Thus creating a cumulative probability of finding this kmer based on occurrences of letters in kmer vs letters in entire file
			 */
			double estimatedProportion = 1; //proportion of finding this base in this kmer.
			for (int i = 0; i < 4; i++) {
				estimatedProportion *= pow(
						(double) baseStatistics[i].Probability,
						(double) kmerBaseStatistics[i].Count);

				DEBUG_STATISTICS(
						cout << "baseStatistics[i].Probability == " << baseStatistics[i].Probability << " raised to the " << kmerBaseStatistics[i].Count << "  == kmerBaseStatistics[i]" << endl;

						cout << "estimatedProportion so far == " << estimatedProportion << endl;);

			}

			//Find the Z score which is the normal binomial distribution from previously calculated values.
			unsigned long long n = *TotalNumSequencesN; //total number of bases in the file.
			unsigned long long x = head->frequency; // x = number of successes that I have had given the number of trials (x <= N)
			long double p = estimatedProportion; // probability of success based on occurrences of letters in kmer vs letters in entire file
			long double q = 1 - p; // probability of failure.
			long double standardDev = sqrt(n * p * q);	//standard deviation
			long double mean = n * p; //average (population mean)
			long double z = (x - mean) / standardDev; //calculate the z score

			DEBUG_STATISTICS(
					cout << n << " = n, " << x << " = x, " << p << " = p, " << q
					<< " = q, " << standardDev << " = standardDev, "
					<< mean << " = mean, " << z << " = z" << endl; );

			/*
			 * If there is no z filtering
			 * Or if z filtering is enabled and the z score of this sequence is above it
			 * Then we can print the data to the file
			 *
			 */
			if (config.zThresholdEnable == 0
					|| ((config.zThresholdEnable > 0)
							&& (abs(z) >= config.zThreshold))) {

				//write the information to the file.
				//start a new line.
				fputc('\n', config.out_file_pointer);

				//print out the sequence that we found.
				for (int i = 0; i < k; i++) {
					fputc(int2base(array[i]), config.out_file_pointer);
				}

				// print out the number of bits to encode a single symbol in the sequence.
				fprintf(config.out_file_pointer, ", %LE", h);

				// print out the number of bits to encode the entire sequence.
				fprintf(config.out_file_pointer, ", %LE", H);

				//print out the number of times that we saw the sequence.
				fprintf(config.out_file_pointer, ", %d", head->frequency);

				//There is a test to see if we can do the normal approximation test or not.
				bool canDoNormalApprox = normal_approx_check(n, p, 1 - p);
				if (canDoNormalApprox == true) {
					//http://www.cplusplus.com/reference/cstdio/printf/ was using %Le
					//if the threshold is not enabled OR (if the threshold is enabled and our z score is a minimum the Z threshold.)
					DEBUG_STATISTICS(
							cout << config.zThresholdEnable
							<< " config.zThresholdEnable, " << z << " = z, "
							<< config.zThreshold << " = config.zThreshold"
							<< endl);
					//print out the Z score value if it is greater than or equal to the threshold.
					fprintf(config.out_file_pointer, ", %LE", z);
				}DEBUG_STATISTICS( else {fprintf(stdout,
									"The sequence did not pass the normal approximation test and was not written to the file.\n");});

				// print higher precision, but the length of long double is undefined and in our experiments, we don't have any duplicate Z scores.
				//fprintf(config.out_file_pointer, ", %.10LE", z);
			}

			//OLD STUFF to verify that our procedure is working step by step.
			//			float_n_choose_k(n, x) * pow((double) 1 - p, (double) n - x) * pow((double) p, (double) x)
			DEBUG_STATISTICS(
					cout << float_n_choose_k(n, x) << "   " << pow((double) 1 - p, (double) n - x) << "   " << pow((double) p, (double) x) << endl; cout << float_n_choose_k(n, x) * pow((double) 1 - p, (double) n - x) * pow((double) p, (double) x) << endl;);

			DEBUG_STATISTICS(
					{
						long double answer = float_n_choose_k(
								*TotalNumSequencesN, head->frequency);
						long double binomialDistribution = answer
						* pow((double ) 1 - estimatedProportion,
								(double ) (*TotalNumSequencesN)
								- head->frequency)
						* pow((double ) estimatedProportion,
								(double ) head->frequency);
						fprintf(config.out_file_pointer, ", %Le",
								binomialDistribution)
						;

						cout
						<< "float_n_choose_k(TotalNumSequencesN, head->frequency) == float_n_choose_k( "
						<< (*TotalNumSequencesN) << ", "
						<< head->frequency << endl
						;

						cout
						<< float_n_choose_k((*TotalNumSequencesN),
								head->frequency) << "   "
						<< pow((double ) 1 - estimatedProportion,
								(double ) (*TotalNumSequencesN)
								- head->frequency) << "   "
						<< pow((double ) estimatedProportion,
								(double ) head->frequency) << endl
						;

						cout
						<< float_n_choose_k((*TotalNumSequencesN),
								head->frequency)
						* pow((double ) 1 - estimatedProportion,
								(double ) (*TotalNumSequencesN)
								- head->frequency)
						* pow((double ) estimatedProportion,
								(double ) head->frequency)
						<< endl
						;
					});

		}			//end if for reaching depth of k
	}			//end else if for head == NULL
}			//end histogram function.
/*
 * This function brings in a pointer to an integer array and shifts its contents left.
 * It also brings in an integer to insert into the array at the right most location.
 */
void shift_left_and_insert(int * const array, const int integer_to_insert) {
	//DEBUG_SHIFT_AND_INSERT(printf("shift_left_and_insert:: received: %c. Starting with array : ",int2base(integer_to_insert)); for (int i = 0; i < config.k; i++) {printf("%c", int2base(*(array + i)));}printf("\n"););

	int i = 0;
	for (i = 0; i < config.k - 1; i++) {
		array[i] = array[i + 1];
	}
	array[i] = integer_to_insert;

	DEBUG_SHIFT_AND_INSERT(
			printf("shift_left_and_insert:: ending with array :                "); for (int i = 0; i < config.k; i++) {printf("%c", int2base(*(array + i)));}printf("\n"););
}
/*
 * This function conforms to the description of this program above by reading a text file and creating a histogram of sequences of length k.
 */
node_t * findKmer(node_t * headNode, unsigned long long * const baseCounter,
		statistics_t * const baseStatistics,
		unsigned long long * const TotalNumSequencesN) {

	// Array to hold the kmer of size k. This ensures that we can hold each sequence, but requires the shifting of data in the array.
	int *kmer = (int*) allocate_array(config.k, sizeof(int));
	int i = 0;

	//fill array by inserting a negative one and testing the functionality of the shift function.
	while (i++ < config.k)
		shift_left_and_insert(kmer, -1);

	//stores the character read in from the file.
	char c = 'A';
	//holds the size of the current sequence. NATTAN would have seqSize 4 before the N was encountered to reset it.
	int seqSize = 0;
	//stores the coded value of the character read from the file. 0-3 are valid, else invalid base.
	int codedBase = 0;

	//check for empty sequence file.
	if (fgetc(config.sequence_file_pointer) == EOF) {
		fprintf(stderr, "Sequence File Is Empty, Ending Program");
		exit(EXIT_FAILURE);
	}
	rewind(config.sequence_file_pointer);

	while ((c = fgetc(config.sequence_file_pointer)) != EOF) {

		/* if the below is true then we have found an identifier which has an entire line of non sequence data */
		if (c == '>') {

			//This line is very important! This means that a > in the file will break a sequence and it will treat the line as a comment.
			seqSize = 0;

			if (config.suppressOutputEnable == 0) {
				fprintf(stdout, "Read %llu bases\n%c", *baseCounter, c);

				while ((c = fgetc(config.sequence_file_pointer)) != '\n') {
					fprintf(stdout, "%c", c);
				}
				fprintf(stdout, "\n");
			} else {
				//dump the current line.
				while ((c = fgetc(config.sequence_file_pointer)) != '\n')
					;
			}
		}

		//ignore newlines, but other invalid base data may be necessary to break the sequence, like > or numbers.
		if (c != '\n') {
			codedBase = base2int(c);
			DEBUG(
					cout << "Reading from file. \n"; cout << "The Base we found is " << c; cout << " which when coded is " << codedBase << "\n";);

			/* If the below is true then we have found an invalid base value thus we must break the sequence apart.
			 * else we have found a valid base value
			 */
			if (codedBase < 0) {
				//any character in the file that is not a newline or a > or preceded by a > will break the sequence
				seqSize = 0;
			} else {

				/* Store the coded base into the kmer to be read later. */
				shift_left_and_insert(kmer, codedBase);
				seqSize++;

				/* If the below is true then that means we have found a valid sequence that is either of k size or greater.
				 * Create a tree data structure where each node is a base encountered.
				 * We also keep track of the total number of bases and the number of each base encountered.
				 */
				if (seqSize > config.k) {

					headNode = tree_create(headNode, kmer, config.k,
							baseStatistics);

					(*baseCounter)++;
					baseStatistics[codedBase].Count++;
					(*TotalNumSequencesN)++;

				} else if (seqSize == config.k) {

					//this case will occur less often than seqSize > config.k
					headNode = tree_create(headNode, kmer, config.k,
							baseStatistics);

					for (int i = 0; i < config.k; i++) {
						baseStatistics[kmer[i]].Count++;
						DEBUG_STATISTICS(
								fprintf(stdout,"i == %d, int2base(kmer[i]) == %c, baseStatistics[kmer[i]].Count == %d.\n",i, int2base(kmer[i]),baseStatistics[kmer[i]].Count));

					}DEBUG_STATISTICS(fprintf(stdout,"\n"));
					(*baseCounter) += seqSize;
					(*TotalNumSequencesN)++;
				} //end detection of a kmer of length k or greater.

			} //end end of sequence detection.
		} //end ignore newline character
	} //end while loop to read the file.

	return headNode;
}
void scratch_function() {

	{
		cout << "simple binomial dist " << endl;
		int n = 10;
		int x = 4;
		double p = 1.0 / 5.0;
		double q = 1 - p;
		double standardDev = sqrt(n * p * q);	//standard deviation
		double median = n * p;
		double z = (x - median) / standardDev;

		cout << float_n_choose_k(n, x) << "   "
				<< pow((double) 1 - p, (double) n - x) << "   "
				<< pow((double) p, (double) x) << endl;
		cout
				<< float_n_choose_k(n, x) * pow((double) 1 - p, (double) n - x)
						* pow((double) p, (double) x) << endl;
		bool canDoNormalApprox = normal_approx_check(n, p, 1 - p);
		cout << " normal distribution possible? " << canDoNormalApprox << endl;
		if (canDoNormalApprox == true) {
			cout << z << " == the z score " << endl;
		}
	}
	{
		cout << "simple binomial dist " << endl;
		int n = 7;
		int x = 4;
		double p = 0.9;

		double q = 1 - p;
		double standardDev = sqrt(n * p * q);	//standard deviation
		double median = n * p;
		double z = (x - median) / standardDev;

		cout << float_n_choose_k(n, x) << "   "
				<< pow((double) 1 - p, (double) n - x) << "   "
				<< pow((double) p, (double) x) << endl;
		cout
				<< float_n_choose_k(n, x) * pow((double) 1 - p, (double) n - x)
						* pow((double) p, (double) x) << endl;
		bool canDoNormalApprox = normal_approx_check(n, p, 1 - p);
		cout << " normal distribution possible? " << canDoNormalApprox << endl;
		if (canDoNormalApprox == true) {
			cout << z << " == the z score " << endl;
		}

	}

	{
		cout << "AAA in a 3mer of uptstreams binomial dist calculation "
				<< endl;
		int n = 90148375;
		int x = 2645164;
		double p = 0.0166418;
		double q = 1 - p;
		double standardDev = sqrt(n * p * q);	//standard deviation
		double median = n * p;
		double z = (x - median) / standardDev;

		cout << float_n_choose_k(n, x) << "   "
				<< pow((double) 1 - p, (double) n - x) << "   "
				<< pow((double) p, (double) x) << endl;
		cout
				<< float_n_choose_k(n, x) * pow((double) 1 - p, (double) n - x)
						* pow((double) p, (double) x) << endl;
		bool canDoNormalApprox = normal_approx_check(n, p, 1 - p);
		cout << " normal distribution possible? " << canDoNormalApprox << endl;
		if (canDoNormalApprox == true) {
			cout << z << " == the z score " << endl;
		}
	}

	{
		cout
				<< "3mer normal distribution to approximate a binomial prob. distribution."
				<< endl;
		int n = 90148375;
		int x = 2645164;
		double p = 0.0166418;

		double q = 1 - p;
		double standardDev = sqrt(n * p * q);	//standard deviation
		double median = n * p;
		double z = (x - median) / standardDev;

		cout << float_n_choose_k(n, x) << "   "
				<< pow((double) 1 - p, (double) n - x) << "   "
				<< pow((double) p, (double) x) << endl;
		cout
				<< float_n_choose_k(n, x) * pow((double) 1 - p, (double) n - x)
						* pow((double) p, (double) x) << endl;
		bool canDoNormalApprox = normal_approx_check(n, p, 1 - p);
		cout << " normal distribution possible? " << canDoNormalApprox << endl;
		if (canDoNormalApprox == true) {
			cout << z << " == the z score " << endl;
		}

	}
	{
		cout
				<< "8mer all A's normal distribution to approximate a binomial prob. distribution."
				<< endl;
		int n = 89697072;
		int x = 151071;
		double p = 1.80524e-05;

		double q = 1 - p;
		double standardDev = sqrt(n * p * q);	//standard deviation
		double median = n * p;
		double z = (x - median) / standardDev;

		cout << float_n_choose_k(n, x) << "   "
				<< pow((double) 1 - p, (double) n - x) << "   "
				<< pow((double) p, (double) x) << endl;
		cout
				<< float_n_choose_k(n, x) * pow((double) 1 - p, (double) n - x)
						* pow((double) p, (double) x) << endl;
		bool canDoNormalApprox = normal_approx_check(n, p, 1 - p);
		cout << " normal distribution possible? " << canDoNormalApprox << endl;
		if (canDoNormalApprox == true) {
			cout << z << " == the z score " << endl;
		}

	}

	fprintf(stderr, " ");
	exit(1);
}
/*
 * This is supposed to be a simple recursive destruction of a tree
 * It currently does NOT free the tree and the reason is UNKNOWN
 */
void destroy(node_t *root) {
	DEBUG_FREE(cout << "attempting destroy " << endl);
	// If we have a non-NULL pointer, we need to
	// recursively free its left and right children,
	// and then free it.
	if (root) {
		DEBUG_FREE(cout << "found leaf " << endl);
		for (int i = 0; i < 4; i++) {
			DEBUG_FREE(cout << "branch " << i << endl);

			destroy(root->nextNodePtr[i]);
			root->nextNodePtr[i] = NULL;

		}DEBUG_FREE(
				cout << "attempting free of base " << int2base(root->base)
				<< endl);
		fflush(stdout);
		free(root);

		DEBUG_FREE(cout << "out of free " << endl);
	}

}
unsigned long int estimate_RAM_usage() {

	if (sizeof(int) < 4 || sizeof(long int) < 8 || sizeof(long long int) < 8) {
		cout
				<< "The normal size of bits is less than expected when this program was written."
				<< endl;
		cout
				<< "This may cause rollover or inaccurate data, so please watch for that "
				<< endl;
		cout
				<< "expected : sizeof(  int) = 4, 	sizeof(long int) = 8, 	sizeof(long long int) = 8, 	sizeof(short int) = 2, 	sizeof(unsigned short int) = 2, 	sizeof(char) = 1, 	sizeof(node_t) = 48, 	sizeof(node_t*) = 8"
				<< endl;
		cout << "found:" << endl;
		printf("sizeof(  int) = %lu\n", sizeof(int));
		printf("sizeof( long int) = %lu\n", sizeof(long int));
		printf("sizeof( long long int) = %lu\n", sizeof(long long int));
		printf("sizeof( short int) = %lu\n", sizeof(short int));
		printf("sizeof(unsigned short int) = %lu\n",
				sizeof(unsigned short int));
		printf("sizeof(char) = %lu\n", sizeof(char));
		printf("sizeof(node_t) = %lu\n", sizeof(node_t));
		printf("sizeof(node_t*) = %lu\n", sizeof(node_t*));
		cout << "Press enter to continue or you can abort the program now."
				<< endl;
		if (config.suppressOutputEnable == 0) {
			getchar();
		}
	}

	//add one for the head node. Calculate RAM usage and Harddrive usage.
	unsigned long int maxNumberOfNodes = 1;
	double n = 1;
	while (n <= config.k) {
		maxNumberOfNodes += pow(4.0, n++);
	}
	if (((sizeof(char) * (config.k + 10)) * maxNumberOfNodes)
			>= (1024 * 1024 * 1024)) {
		cout
				<< ((sizeof(char) * (config.k + 10)) * maxNumberOfNodes)
						/ (double) (1024 * 1024 * 1024) << " gibibytes";
	} else {
		cout
				<< ((sizeof(char) * (config.k + 10)) * maxNumberOfNodes)
						/ (double) (1024 * 1024) << " mibibytes";
	}

	cout << " of disk usage and ";

	if (maxNumberOfNodes * sizeof(node_t) >= (1024 * 1024 * 1024)) {
		cout
				<< (maxNumberOfNodes * sizeof(node_t)
						/ (double) (1024 * 1024 * 1024))
				<< " gibibytes of RAM usage likely" << endl;
		;
		cout << "We are stopping here to make sure that is ok with you!"
				<< endl;
		cout << "Hit enter to proceed or else abort the program." << endl;
		if (config.suppressOutputEnable == 0) {
			getchar();
		}
	} else {
		cout << (maxNumberOfNodes * sizeof(node_t) / (double) (1024 * 1024))
				<< " mibibytes of RAM usage likely" << endl;
	}
	return maxNumberOfNodes;
}
int main(int argc, char *argv[]) {

	/* Deal with command line arguments */
	DEBUG(
			int currentArgument =0; while(currentArgument < argc) {fprintf(stdout, "argv[%d]== %s\n",currentArgument, *(argv+currentArgument)); /* %s instead of %c and drop [i]. */
				/* Next arg. */
				currentArgument++;}fprintf(stdout, "\n");)

	init_conf();
	usage();
	while (!parse_arguments(argc, argv))
		usage();
	print_conf(argc);

	unsigned long int maxNumberOfNodes = estimate_RAM_usage(); //Most number of nodes that can be created in memory.

	/* Begin the procedure to extract valid sequences from file */
	fprintf(stdout, "!!!Find The KMER!!!\n");
	fprintf(stdout, "Reading sequence from file\n");
	fprintf(stdout,
			"     2858658142 bases in the reference genome FYI.\nThat is 2,858,658,142 by the way.\n");

	/* variable for the root of the tree */
	node_t * headNode = NULL;
	unsigned long long baseCounter = 0; //number of bases that fit into a kmer in the entire file. GATTACA has baseCounter = 7 if k <= 7
	unsigned long long TotalNumSequencesN = 0; //number of kmers found. if k = 2 then GATA has N=3.
	statistics_t baseStatistics[4] = { 0 };

	headNode = findKmer(headNode, &baseCounter, baseStatistics,
			&TotalNumSequencesN);

	statistics(&baseCounter, baseStatistics, &TotalNumSequencesN,
			&maxNumberOfNodes);

	fprintf(stdout, "Now creating histogram.\n");

	//create a temporary array for the recursive function to keep as scratch memory to hold the sequence.
	//int* histogram_temp = (int*) allocate_array(config.k, sizeof(int));
	int* histogram_temp = (int*) malloc(config.k * sizeof(int));
	if (!histogram_temp) {
		fprintf(stderr, "allocate_array():: memory allocation failed\n");
		exit(EXIT_FAILURE);
	}
	//Initialize memory.
	for (int i = 0; i < config.k; i++) {
		*(histogram_temp + i) = -1;
	}

	/* Output the occurrence of every sequence of length k */

	histo_recursive(headNode, histogram_temp, 0, config.k, &baseCounter,
			baseStatistics, &TotalNumSequencesN);

	//Begin cleanup and closing of files.
	free(histogram_temp);
	histogram_temp = NULL;

	destroy(headNode);

	DEBUG(fprintf(stdout, "\n"));
	fprintf(stdout, "histogram creation finished.\n");

	if (fclose(config.out_file_pointer) == EOF) {
		fprintf(stderr,
				"Out file close error! This is not expected and might mean the data was not written to the file properly before the close.\n");
	}

	fprintf(stdout,
			"Your file can be found in the current directory as: \n    %s\n",
			config.out_file);

	if (fclose(config.sequence_file_pointer) == EOF) {
		fprintf(stderr,
				"Sequence file close error! This is likely ok though.\n");
	}
	//Do not put any code after this point.
	//Since the sequence file name is a constant char pointer, we cannot use free according to the standard since we cannot set it to null later.
	//This will be fixed in the destructor of the config object later.
	free((char *) config.sequence_file);
	free(config.out_file);
	fprintf(stdout, "End of program was reached properly.\n\n");
	fprintf(stderr, " "); //simply to trigger error to notify the eclipse that we are done.
//
//	if (config.suppressOutputEnable == 0) {
//		getchar();
//	}
	return 0;
}
