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
// Bugs		   : Known bugs include memory leaks and the tree is created with more depth than intended.
// Compile     : To compile perform: g++ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/findKmer.d" -MT"src/findKmer.d" -o "src/findKmer.o" "../src/findKmer.cpp"
//============================================================================
using namespace std;
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h> //for strcmp(string1,string2) string comparison returns a 0 if they are the same.
#include <stdlib.h> //malloc is in this.
// basic file operations
#include <fstream>

#define MAX_LINE 1001
#define DEBUG(x) //x
#define DEBUG_TREE_CREATE(x) //x
#define DEBUG_HISTO_AND_FREE_RECURSIVE(x) //x
// Data structure for a tree.
struct node_t {
	int base;
	struct node_t *nextNodePtr[4];
	unsigned int counter;
};

/* configuration */
static struct conf {
	char *sequence_file;
	FILE *sequence_file_pointer;
	char *out_file;
	FILE *out_file_pointer;
	int k;
} config;
//enumeration of possible bases in alphabetical order.
enum Base {
	A = 0, C = 1, G = 2, T = 3
};

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
	config.k = 5;

}

/* This function fills in any gaps in the configuration file.*/
void set_default_conf() {

	if (!config.sequence_file) {
		config.sequence_file = "homo_sapiensupstream.fas";
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

	if (!config.k) {
		config.k = 5;
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
}

/* usage */
static void usage() {
	fprintf(stdout, "\n");
	fprintf(stdout, "Usage: findKmer [options]\n");
	fprintf(stdout,
			"             [--parse|-p <sequence_file.txt>] \n               File with DNA sequence data.\n               File must be in current directory.\n                Default is homo_sapiensupstream.fas.\n\n");
	fprintf(stdout,
			"             [--export|-e  <out_file.csv>] \n               File to output histogram data to.\n                Default output file name is dynamic.\n\n");
	fprintf(stdout,
			"             [--ksize|-k  <k>] \n               Size of sequence for histogram.\n                Default is 5.\n\n");
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
		if (node->nextNodePtr[base] == 0) {
			fprintf(stderr,
					"!!! COUNTER ROLLOVER DETECTED! \nIncrease the number of bits used for the counter variable");
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
node_t* tree_create(node_t* head, int* array, int k) {
	if (array != NULL) {
		if (head == NULL) {
			DEBUG_TREE_CREATE(fprintf(stdout, "-Creating the head of the tree!\n"));
			head = node_create('H');
		}
		node_t *currentNode = head;
		for (int i = 0; i < k; i++) {
			DEBUG_TREE_CREATE(fprintf(stdout, "-Moving into a branch on depth %d\n", i));
			currentNode = node_branch_enter_and_create(currentNode, array[i]);
		}
	}
	return head;
}

/*
 * Histogram and free can be recursive for low numbers of K.
 * If K becomes too high then we may run out of stack/heap memory.
 *
 */
void histo_and_free_recursive(node_t* head, int* array, int *k, int depth) {
	DEBUG_HISTO_AND_FREE_RECURSIVE(
			fprintf(stderr, "histo&free @ depth %d of %d has %d\n",depth,*k,head != NULL?head->base:-1 ));
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
			histo_and_free_recursive(head->nextNodePtr[i], array, k, depth + 1);
			//TODO verify that we are freeing memory.
			deallocate_array((void**) &head->nextNodePtr[i]);
		}

		if (depth == (*k)) {
			for (int i = 0; i < *k; i++) {
				DEBUG(fprintf(stdout, "%c", int2base(array[i])));
				fputc(int2base(array[i]), config.out_file_pointer);
			}DEBUG(fprintf(stdout, ", %d\n", head->counter));
			fprintf(config.out_file_pointer, ", %d\n", head->counter);
		}

	}
}

//for testing purposes.
void random_array(int sizeOfArray) {
	//Initializing array to test code. this will come from the pre processed line
	int *array = (int*) allocate_array(sizeOfArray, sizeof(int));
	int* currentArrayPosition = array;
	srand(time(NULL)); //seeding the random function.
	for (int i = 0; i < sizeOfArray; i++) {
		*currentArrayPosition = (rand() % 4);
		//fprintf(stdout, " currentArrayPosition %d has == %d\n", i,	*currentArrayPosition);
		currentArrayPosition++;
	}
	deallocate_array((void**) &array);
}

int main(int argc, char *argv[]) {

	//command line arguments.
	DEBUG(
			int currentArgument =0; while(currentArgument < argc) {fprintf(stdout, "argv[%d]== %s\n",currentArgument, *(argv+currentArgument)); /* %s instead of %c and drop [i]. */
				/* Next arg. */
				currentArgument++;}fprintf(stdout, "\n");)
	usage();
	init_conf();
	while (!parse_arguments(argc, argv))
		usage();
	print_conf();

	fprintf(stdout, "!!!Find The KMER!!!\n");

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
		fprintf(config.out_file_pointer, "Sequence, Frequency\n");
	} else {
		fprintf(stderr, "Out file failed to open\n\n");
		exit(EXIT_FAILURE);
	}

	fprintf(stdout, "Reading sequence from file\n");

	//variables for the nodes
	node_t * headNode = NULL;

	//extract valid sequences from file
	unsigned long long sequenceNumber = 0;
	int integerBuffer[MAX_LINE]; //todo commandline arg?
	char c = 'A';
	int seqSize = 0;
	int codedBase = 0;

	if (fgetc(config.sequence_file_pointer) == EOF) {
		fprintf(stderr, "Sequence File Is Empty, Ending Program");
		exit(EXIT_FAILURE);
	}
	rewind(config.sequence_file_pointer);

	while (c != EOF) {
		c = fgetc(config.sequence_file_pointer);
		codedBase = base2int(c);
		DEBUG(
				cout << "Reading from file. \n"; cout << "The Base we found is " << c; cout << " which when coded is " << codedBase << "\n";);

		//ignore newlines, but other invalid base data may be necessary to break the sequence, like > or numbers.
		if (c != '\n') {

			//If the below is true then that means we have found a valid sequence that is either of k size or greater.
			if (seqSize >= config.k && codedBase < 0) {
				DEBUG(
						cout << "!!!Found " << sequenceNumber << " valid sequences so far!!!" << endl; cout << seqSize << " Length Sequence Found \n";

						for (int n = 0; n < seqSize; n++) cout << int2base(integerBuffer[n]); cout << endl;);

				++sequenceNumber;
				int* currentArrayPosition = integerBuffer;
				//traverse the integer array and create a tree.
				int stopPoint = seqSize - config.k + 1;
				for (int i = 0; i < stopPoint; i++) {

					DEBUG(
							fprintf(stdout, "Extracting sequence: "); int z = 0; while (z < config.k) {fprintf(stdout, "%c ",int2base(*(currentArrayPosition + z)) ); z++;}fprintf(stdout, "\n"););

					headNode = tree_create(headNode, currentArrayPosition++,
							config.k + 1);
				}

			}
			if (codedBase < 0) {
				seqSize = 0;
			} else {
				if (seqSize > MAX_LINE) {
					fflush(stdout);
					fprintf(stderr,
							"\n!!!!Index %d has exceeded array bounds of %d\n",
							seqSize, MAX_LINE);
					fflush(stderr);
				}
				integerBuffer[seqSize] = codedBase;
				seqSize++;
			}
		}

	}

	fprintf(stdout, "Found %d valid sequences.\n", sequenceNumber);
	fprintf(stdout, "Now creating histogram.\n");

	//create a temporary array for the recursive function to keep as scratch memory to hold the sequence.
	int* histogram_temp = (int*) allocate_array(config.k, sizeof(int));
	//Initialize memory.
	for (int i = 0; i < config.k; i++) {
		*(histogram_temp + i) = -1;
	}

	histo_and_free_recursive(headNode, histogram_temp, &config.k, 0);
	deallocate_array((void**) &headNode);
	//TODO free the histogram_temp array...This kept throwing errors on me. free(histogram_temp);
	DEBUG(fprintf(stdout, "\n"));
	cout << "histogram creation finished." << endl;

	if (fclose(config.sequence_file_pointer) == EOF) {
		fprintf(stderr,
				"Sequence file close error! This is likely ok though.\n");
	}
	if (fclose(config.out_file_pointer) == EOF) {
		fprintf(stderr,
				"Out file close error! This is not expected and might mean the data was not written to the file properly before the close.\n");
	}
	cout << "Your file can be found in the current directory as: " << "\n    "
			<< config.out_file << endl;
	fprintf(stdout, "End of program\n");
	return 0;
}
