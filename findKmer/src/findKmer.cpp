/*
 * Copyright (c) 2007 Michela Becchi and Washington University in St. Louis.
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
 *    3. The name of the author or Washington University may not be used
 *       to endorse or promote products derived from this source code
 *       without specific prior written permission.
 *    4. Conditions of any other entities that contributed to this are also
 *       met. If a copyright notice is present from another entity, it must
 *       be maintained in redistributions of the source code.
 *
 * THIS INTELLECTUAL PROPERTY (WHICH MAY INCLUDE BUT IS NOT LIMITED TO SOFTWARE,
 * FIRMWARE, VHDL, etc) IS PROVIDED BY  THE AUTHOR AND WASHINGTON UNIVERSITY
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR WASHINGTON UNIVERSITY
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
// Expects     : The parse file is to be formatted in this way:
//                 '>' designates the beginning of an ID and newlines ('\n') designates the end of an ID
//                 Newlines ('\n') are otherwise ignored completely
//                 Sequences are broken by any letter (other than A, C, T or G)
//                 The letter N is a known letter that occurs in a sequence and is designated as the end of a sequence if found.
//                 All valid bases are in upper case (A, C, T or G) not lower case.
// Bugs        : Known bugs include memory leaks and the tree is created with more depth than intended.
// Compile     : To compile perform: g++ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/findKmer.d" -MT"src/findKmer.d" -o "src/findKmer.o" "../src/findKmer.cpp"
//============================================================================
using namespace std;
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h> //for strcmp(string1,string2) string comparison returns a 0 if they are the same.
#include <stdlib.h> //malloc is in this.
#include <fstream> // basic file operations

//must be a string literal with file extension included.
#define DEFAULT_SEQUENCE_FILE_NAME "homo_sapiensupstream.fas"
//#define DEFAULT_SEQUENCE_FILE_NAME "Full_homo_sapiens.fa"
//#define DEFAULT_SEQUENCE_FILE_NAME "test.txt"
#define DEFAULT_K_VALUE 3

//This option will NOT make the tree in memory and will NOT create a valid histogram. it is only to save memory and count base occurances.
#define COUNT_BASES_ONLY(x) x

//debugging
#define DEBUG(x) //x
#define DEBUG_TREE_CREATE(x) //x
#define DEBUG_HISTO_AND_FREE_RECURSIVE(x) //x
#define DEBUG_SHIFT_AND_INSERT(x) //x
#define DEBUG_STATISTICS(x) //x

//Holds the statistics for a base
struct statistics_t {
	unsigned int Count; //number of said base encountered.
	long double Probability; //probability that this base will be encountered out of all bases.
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
	unsigned char base;
	struct node_t *nextNodePtr[4];
	unsigned int counter;
};

/* structure definition for configuration of file names, pointers, and length of k.*/
static struct conf {
	char *sequence_file; //holds the string representation of the file name.
	FILE *sequence_file_pointer; //holds the FILE pointer to the file itself
	char *out_file;		//holds the string representation of the file name.
	FILE *out_file_pointer;  //holds the FILE pointer to the file itself
	int k; //holds the length of k for the size of the sequence to be recorded.
} config; /* Config is a GLOBAL VARIABLE for configuration of file names, pointers, and length of k.*/

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

/* check that the given file can be read/written */
void check_file(char *filename, char *mode) {
	FILE *file = fopen(filename, mode);
	if (file == NULL) {
		fprintf(stderr, "Unable to open file %s in %s mode\n", filename, mode);
		exit(EXIT_FAILURE);
	} else
		fclose(file);
}

/* initialize the configuration */
void init_conf() {
	config.sequence_file = NULL;
	config.sequence_file_pointer = NULL;
	config.out_file = NULL;
	config.out_file_pointer = NULL;
	config.k = 0;

}

/* This function fills in any gaps in the configuration file.*/
void set_default_conf() {

	if (!config.sequence_file) {
		config.sequence_file = DEFAULT_SEQUENCE_FILE_NAME;
	}

	if (!config.k) {
		config.k = DEFAULT_K_VALUE;
	}

	if (!config.k) {
		fprintf(stdout, "k must be greater than zero. Ending program.\n\n");
		exit(EXIT_FAILURE);
	}

	if (!config.out_file) {
		char* nameOfFile = "mer_Historam_Of_";
		char* outFileExension = ".csv";
		config.out_file = (char*) allocate_array(
				strlen("999") + strlen(nameOfFile)
						+ strlen(config.sequence_file)
						+ strlen(outFileExension), sizeof(char));
		sprintf(config.out_file, "%d%s%s%s", config.k, nameOfFile,
				config.sequence_file, outFileExension);
	}

}

/* print the configuration */
void print_conf() {
	fprintf(stdout, "\nATTEMPTING CONFIGURATION: \n");
	set_default_conf();
	if (config.sequence_file)
		fprintf(stdout, "- sequence_file file: %s\n", config.sequence_file);
	if (config.out_file)
		fprintf(stdout, "- export file: %s\n", config.out_file);
	if (config.k)
		fprintf(stdout, "- k size: %d\n", config.k);
	fprintf(stdout, "\n");

	/* Double check configuration */
	if (config.k < 0) {
		fprintf(stderr,
				"%d is not a valid value for k. Please select a number greater than zero\n",
				config.k);
		exit(EXIT_FAILURE);
	}

	if ((config.sequence_file_pointer = fopen(config.sequence_file, "r"))
			!= NULL) {
		fprintf(stdout, "Sequence file opened properly\n");

	} else {
		fprintf(stderr, "Sequence file failed to open\n\n");
		exit(EXIT_FAILURE);
	}
	if ((config.out_file_pointer = fopen(config.out_file, "w")) != NULL) {
		fprintf(stdout, "Out file opened properly\n");
		fprintf(config.out_file_pointer, "Sequence, Frequency");
	} else {
		fprintf(stderr, "Out file failed to open\n\n");
		exit(EXIT_FAILURE);
	}
}

/* usage */
static void usage() {
	fprintf(stdout, "\n");
	fprintf(stdout, "Usage: findKmer [options]\n");
	fprintf(stdout,
			"             [--parse|-p <sequence_file.txt>] \n               File with DNA sequence data.\n               File must be in current directory.\n                Default is %s.\n\n",
			DEFAULT_SEQUENCE_FILE_NAME);
	fprintf(stdout,
			"             [--export|-e  <out_file.csv>] \n               File to output histogram data to.\n                Default output file name is dynamic.\n\n");
	fprintf(stdout,
			"             [--ksize|-k  <k>] \n               Size of sequence for histogram.\n                Default is %d.\n\n",
			DEFAULT_K_VALUE);
	fprintf(stdout, "\n");

}

int parse_arguments(int argc, char **argv) {
	int i = 1;
	if (argc < 2) {
		return 1;

	} else {

		while (i < argc) {
			if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
				usage();
				return 0;
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
			} else if (strcmp(argv[i], "-k") == 0) {
				i++;
				int k = atoi(argv[i]);
				if (k < 0) {
					fprintf(stderr,
							"%d is not a valid value for k. Please select a number greater than zero",
							k);
				}
				config.k = k;
			} else {
				fprintf(stderr, "Ignoring invalid option %s\n", argv[i]);
			}
			i++;
		}
	}

	return 1;
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
	node_t* node = (node_t*) allocate_array(1, sizeof(node_t));
	node->base = base;
	node->counter = 1;
	node->nextNodePtr[0] = NULL;
	node->nextNodePtr[1] = NULL;
	node->nextNodePtr[2] = NULL;
	node->nextNodePtr[3] = NULL;
	return node;
}
/*
 * Adds a branch if one doesn't already exist.
 * Increments the counter for the node we are going to step into.
 * Returns the pointer to the next node.
 */
node_t* node_branch_enter_and_create(node_t* node, int base) {
	DEBUG_TREE_CREATE(
			fprintf(stdout, "..node_branch_enter_and_create for base %c\n",
					int2base(base)));
	if (node->nextNodePtr[base] == NULL) {
		DEBUG_TREE_CREATE(fprintf(stdout, "***Creating Node.\n"));
		node->nextNodePtr[base] = node_create(base);
	} else {
		node->nextNodePtr[base]->counter++;
		if (node->nextNodePtr[base]->counter == 0) {
			fprintf(stderr,
					"\n\n!!! COUNTER ROLLOVER DETECTED! \nIncrease the number of bits used for the counter variable if you have the source code, else use a smaller sequence file.\n\n");
			fprintf(stdout,
					"\n\n!!! COUNTER ROLLOVER DETECTED! \nIncrease the number of bits used for the counter variable if you have the source code, else use a smaller sequence file.\n\n");
			exit(EXIT_FAILURE);
		}DEBUG_TREE_CREATE(
				fprintf(stdout, "+++Incrementing counter to %d.\n",
						node->nextNodePtr[base]->counter));
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
node_t* tree_create(node_t* head, int* array, int k,
		statistics_t *baseStatistics) {
	if (array != NULL) {
		if (head == NULL) {
			DEBUG_TREE_CREATE(fprintf(stdout, "-Creating the head of the tree!\n"));
			head = node_create('H');
		}
		node_t *currentNode = head;

		/*
		 * Traverse the given integer array that is of k size.
		 * For each base within the integer array, create a node and enter the created node.
		 * Also record how often each base (A, C, G or T).
		 */
		for (int i = 0; i < k; i++) {
			DEBUG_TREE_CREATE(fprintf(stdout, "-Moving into a branch on depth %d\n", i));
			COUNT_BASES_ONLY(
					currentNode = node_branch_enter_and_create(currentNode,
							array[i]));

		}
	}
	return head;
}

/*
 * Histogram and free can be recursive for low numbers of K.
 * If K becomes too high then we may run out of stack/heap memory.
 * It will traverse the tree starting at the head and work its way to the depth provided by k.
 * When it reaches the depth k, it will print out the sequence it saw from the root to this point.
 * It will then free that node and work its way back, hopefully freeing every node as it goes along.
 * The implementation of freeing the tree is not implemented here.
 */
void histo_recursive(node_t* head, int* array, int depth, int k,
		unsigned long long *baseCounter, statistics_t *baseStatistics) {
	DEBUG_HISTO_AND_FREE_RECURSIVE(
			fprintf(stdout, "histo&free @ depth %d of %d has %d\n",depth,*k,head != NULL?head->base:-1 ));

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

		for (int i = 0; i < 4; i++) {
			DEBUG_HISTO_AND_FREE_RECURSIVE(
					fprintf(stdout, "histo&free @ depth %d of %d has %d checking branch %d\n",depth,*k,head != NULL?head->base:-1,i ));
			histo_recursive(head->nextNodePtr[i], array, depth + 1, k,
					baseCounter, baseStatistics);
			//free memory.
			deallocate_array((void**) &head->nextNodePtr[i]);
		}

		if (depth == (k)) {
			unsigned int kmerBaseStatistics[4] = { 0 };
			fputc('\n', config.out_file_pointer);

			cout << "traversing kmer" << endl;
			for (int i = 0; i < k; i++) {
				cout << "i == " << i << endl;
				cout << "array[i] == "<< array[i]<< endl;
				cout << "kmerBaseStatistics[array[i]] == " << kmerBaseStatistics[array[i]] << endl;
				kmerBaseStatistics[array[i]]++;

				DEBUG(fprintf(stdout, "%c", int2base(array[i])));
				fputc(int2base(array[i]), config.out_file_pointer);

			}DEBUG(fprintf(stdout, ", %d\n", head->counter));

			fprintf(config.out_file_pointer, ", %d", head->counter);
			double estimatedProportion = 1;
			for(int i = 0; i < 4; i++){
				estimatedProportion *= pow((double)baseStatistics[i].Probability,(double)kmerBaseStatistics[i]);

				cout << "baseStatistics[i].Probability == "<<baseStatistics[i].Probability << " raised to the " << kmerBaseStatistics[i] << "  == kmerBaseStatistics[i]"<< endl;

				cout << "estimatedProportion so far == " << estimatedProportion<<endl;
			}



				int TotalNumSequencesN = 10000;

				long double answer = float_n_choose_k(TotalNumSequencesN, head->counter);
				cout << float_n_choose_k(TotalNumSequencesN, head->counter) << "   " << pow(1 - estimatedProportion, TotalNumSequencesN-head->counter) <<"   " <<  pow(estimatedProportion, head->counter)<<endl;
				cout << float_n_choose_k(TotalNumSequencesN, head->counter) * pow(1 - estimatedProportion, TotalNumSequencesN-head->counter) * pow(estimatedProportion, head->counter) << endl;
				fprintf(stderr," ");
				exit(1);


		}			//end if for reaching depth of k
	}			//end else if for head == NULL
}			//end histogram function.

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

void shift_left_and_insert(int* array, int integer_to_insert) {
	DEBUG_SHIFT_AND_INSERT(
			printf("shift_left_and_insert:: received: %c. Starting with array : ",int2base(integer_to_insert)); for (int i = 0; i < config.k; i++) {printf("%c", int2base(*(array + i)));}printf("\n"););

	int i = 0;
	for (i = 0; i < config.k - 1; i++) {
		array[i] = array[i + 1];
	}
	array[i] = integer_to_insert;

	DEBUG_SHIFT_AND_INSERT(
			printf("shift_left_and_insert:: ending with array :                "); for (int i = 0; i < config.k; i++) {printf("%c", int2base(*(array + i)));}printf("\n"););
}

node_t * findKmer(node_t * headNode, unsigned long long *baseCounter,
		statistics_t *baseStatistics) {

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

	if (fgetc(config.sequence_file_pointer) == EOF) {
		fprintf(stderr, "Sequence File Is Empty, Ending Program");
		exit(EXIT_FAILURE);
	}
	rewind(config.sequence_file_pointer);

	while (c != EOF) {
		c = fgetc(config.sequence_file_pointer);

		/* if the below is true then we have found an identifier which has an entire line of non sequence data */
		if (c == '>') {
			fprintf(stdout, "Read %llu bases so far.\n%c", *baseCounter, c);

			while ((c = fgetc(config.sequence_file_pointer)) != '\n') {
				fprintf(stdout, "%c", c);
			}
			fprintf(stdout, "\n");
			seqSize = 0;
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
				} else if (seqSize == config.k) {
					//this case will occur less often than seqSize > config.k so put it after in an else statement
					headNode = tree_create(headNode, kmer, config.k,
							baseStatistics);

					DEBUG_STATISTICS(fprintf(stdout,"sequence size equals k, valid sequence found!\n"));
					for (int i = 0; i < config.k; i++) {
						baseStatistics[kmer[i]].Count++;
						DEBUG_STATISTICS(fprintf(stdout,"i == %d, int2base(kmer[i]) == %c, baseStatistics[kmer[i]].Count == %d.\n",i, int2base(kmer[i]),baseStatistics[kmer[i]].Count));
					}DEBUG_STATISTICS(fprintf(stdout,"\n"));
					(*baseCounter) += seqSize;
				}

			} //end end of sequence detection.
		} //end ignore newline character
	} //end while loop to read the file.

	return headNode;
}

void statistics(unsigned long long *baseCounter, statistics_t *baseStatistics) {
	//Statistics section Still needs to be verified and expanded.
	fprintf(stdout,
			"Statistics of occurrences and probability of A, C, G and T respectively: \n");
	for (int i = 0; i < 4; i++) {
		DEBUG_STATISTICS(fprintf(stdout, "%u count / %llu baseCounter\n",baseStatistics[i].Count , *baseCounter));

		fprintf(stdout, "%u", baseStatistics[i].Count);

		baseStatistics[i].Probability = (double) baseStatistics[i].Count
				/ *baseCounter;
		if (baseStatistics[i].Probability == 0.0) {
			fprintf(stdout, "Division overflow detected in statistics.\n");
			exit(EXIT_FAILURE);
		}

		fprintf(stdout, ", %Lf\n", baseStatistics[i].Probability);
	}

	fprintf(stdout, "Found %llu valid bases total INSIDE sequences >= k.\n",
			*baseCounter);
}

int main(int argc, char *argv[]) {


	DEBUG(printf("sizeof(  int) = %lu\n",sizeof( int));
			printf("sizeof( short int) = %lu\n",sizeof( short int));
			printf("sizeof(unsigned short int) = %lu\n",sizeof(unsigned short int));
			printf("sizeof(char) = %lu\n",sizeof(char)));

	/* Deal with command line arguments */
	DEBUG(
			int currentArgument =0; while(currentArgument < argc) {fprintf(stdout, "argv[%d]== %s\n",currentArgument, *(argv+currentArgument)); /* %s instead of %c and drop [i]. */
				/* Next arg. */
				currentArgument++;}fprintf(stdout, "\n");)
	usage();
	init_conf();
	while (!parse_arguments(argc, argv))
		usage();
	print_conf();

	/* Begin the procedure to extract valid sequences from file */
	fprintf(stdout, "!!!Find The KMER!!!\n");
	fprintf(stdout, "Reading sequence from file\n");
	fprintf(stdout,
			"     2858658142 bases in the reference genome FYI.\nThat is 2,858,658,142 by the way.\n");

	/* variables for the root of the tree */
	node_t * headNode = NULL;
	unsigned long long baseCounter = 0;
	statistics_t baseStatistics[4] = { 0 };

	headNode = findKmer(headNode, &baseCounter, baseStatistics);

	statistics(&baseCounter, baseStatistics);

	fprintf(stdout, "Now creating histogram.\n");

	//create a temporary array for the recursive function to keep as scratch memory to hold the sequence.
	int* histogram_temp = (int*) allocate_array(config.k, sizeof(int));
	//Initialize memory.
	for (int i = 0; i < config.k; i++) {
		*(histogram_temp + i) = -1;
	}

	/* print out the occurrence of every sequence of length k */
	histo_recursive(headNode, histogram_temp, 0, config.k, &baseCounter,
			baseStatistics);

	//TODO free the histogram_temp array...This kept throwing errors on me. free(histogram_temp);
	DEBUG(fprintf(stdout, "\n"));
	fprintf(stdout, "histogram creation finished.\n");

	deallocate_array((void**) &headNode);

	if (fclose(config.sequence_file_pointer) == EOF) {
		fprintf(stderr,
				"Sequence file close error! This is likely ok though.\n");
	}
	if (fclose(config.out_file_pointer) == EOF) {
		fprintf(stderr,
				"Out file close error! This is not expected and might mean the data was not written to the file properly before the close.\n");
	}
	fprintf(stdout,
			"Your file can be found in the current directory as: \n    %s\n",
			config.out_file);
	fprintf(stdout, "End of program was reached properly.\n\n");
	fprintf(stderr, " "); //simply to trigger error to notify the user that we are done.
	return 0;
}
