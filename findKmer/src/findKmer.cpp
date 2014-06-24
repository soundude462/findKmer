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
#define DEBUG(x) x
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

/* check that the given file can be read/written */
void check_file(char *filename, char *mode) {
	FILE *file = fopen(filename, mode);
	if (file == NULL) {
		printf("Unable to open file %s in %s mode\n", filename, mode);
		exit(1);
	} else
		fclose(file);
}
/* initialize the configuration */
void init_conf() {
	config.sequence_file = "homo_sapiensupstream.fas";
	config.sequence_file_pointer = NULL;
	config.out_file = "output.txt";
	config.out_file_pointer = NULL;
	config.k = 5;

}
/* print the configuration */
void print_conf() {
	fprintf(stderr, "\nCONFIGURATION: \n");
	if (config.sequence_file)
		fprintf(stderr, "- sequence_file file: %s\n", config.sequence_file);
	if (config.out_file)
		fprintf(stderr, "- export file: %s\n", config.out_file);
	if (config.k)
		fprintf(stderr, "- k size: %d\n", config.k);
}
//enumeration of possible bases in alphabetical order.
enum Base {
	A = 0, C = 1, G = 2, T = 3
};

/* usage */
static void usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: findKmer [options]\n");
	fprintf(stderr,
			"             [--parse|-p <sequence_file>] file with DNA sequence data with default as homo_sapiensupstream.fas. File must be in current directory.\n");
	fprintf(stderr,
			"             [--export|-e  <out_file>] file to output histogram data to with default as output.txt.\n");
	fprintf(stderr,
			"             [--ksize|-k  <k>] size of sequence for histogram with default as 5.\n");
	fprintf(stderr, "\n");
	exit(0);
}

int parse_arguments(int argc, char **argv) {
	int i = 1;
	if (argc < 7) {
		usage();
		return 0;
	}
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
	return 1;
}

void *allocate_array(int size, size_t element_size) {
	void *mem = malloc(size * element_size);
	if (!mem) {
		printf("allocate_array():: memory allocation failed\n");
		exit(1);
	}
	return mem;
}

void *reallocate_array(void *array, int size, size_t element_size) {
	void *new_array = realloc(array, element_size * size);
	if (!new_array) {
		printf("reallocate_array():: memory reallocation failed\n");
		exit(1);
	}

	return new_array;
}

void deallocate_array(void** array) {
	free(*array);
	*array = NULL;
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
	DEBUG_TREE_CREATE(printf("..node_branch_enter_and_create for base %d\n", base));
	if (node->nextNodePtr[base] == NULL) {
		DEBUG_TREE_CREATE(printf("***Creating Node.\n"));
		node->nextNodePtr[base] = node_create(base);
	} else {
		node->nextNodePtr[base]->counter++;
		DEBUG_TREE_CREATE(
				printf("+++Incrementing counter to %d.\n",
						node->nextNodePtr[base]->counter));
	}DEBUG_TREE_CREATE(printf("..returning next base pointer.\n"));
	return node->nextNodePtr[base];
}

/*
 * Brings in a pointer to head of the tree, an integer array and the size k of the array
 * Checks to see if the array exists, and if the tree already exists.
 * Breaks if array does not exist, but creates the head node if tree does not exist.
 * Traverses the array and creates the tree based on what it finds.
 *
 */
node_t* tree_create(node_t* head, int* array, int k) {
	if (array != NULL) {
		if (head == NULL) {
			DEBUG_TREE_CREATE(printf("-Creating the head of the tree!\n"));
			head = node_create('H');
		}
		node_t *currentNode = head;
		for (int i = 0; i < k; i++) {
			DEBUG_TREE_CREATE(printf("-Moving into a branch on depth %d\n", i));
			currentNode = node_branch_enter_and_create(currentNode, array[i]);
		}
	}
	return head;
}
/*
 * Histogram and free can be recursive for low numbers of K.
 * If K becomes too high then we may run out of stack/heap memory.
 * The histogram has been implemented but the free part of the function has not.
 */
void histo_and_free_recursive(node_t* head, int* array, int *k, int depth) {
	DEBUG_HISTO_AND_FREE_RECURSIVE(printf("histo&free @ depth %d of %d has %d\n",depth,*k,head != NULL?head->base:-1 ));
	if (head == NULL) {
		DEBUG_HISTO_AND_FREE_RECURSIVE(printf("histo_and_free::Head == NULL. Leaf found.\n"));
	} else if (depth == (*k) + 1) {
		for (int i = 0; i < *k; i++) {
			printf("%d", array[i]);
		}
		printf(", %d\n", head->counter);
	} else {
		if (head->base == 'H') {
			printf("Head found.\n\n");
		} else {
			array[depth - 1] = head->base;
			//printf("%d", head->base);
		}

		for (int i = 0; i < 4; i++) {
			DEBUG_HISTO_AND_FREE_RECURSIVE(printf("histo&free @ depth %d of %d has %d checking branch %d\n",depth,*k,head != NULL?head->base:-1,i ));
			histo_and_free_recursive(head->nextNodePtr[i], array, k, depth + 1);
			//TODO verify that we are freeing memory.
			deallocate_array((void**) &head->nextNodePtr[i]);
		}

	}
}

int char2int(char base) {
	int integer = -1;
	if (base == 'A') {
		integer = 0;
	} else if (base == 'C') {
		integer = 1;
	} else if (base == 'G') {
		integer = 2;
	} else if (base == 'T') {
		integer = 3;
	} else if (base == 'N') {
		integer = -2;
	} else if (base == '\n') {
		integer = -3;
	}
	return (integer);
}

int main(int argc, char *argv[]) {

	//command line arguments.
	DEBUG(
			int currentArgument =0; while(currentArgument < argc) { printf("argv[%d]== %s\n",currentArgument, *(argv+currentArgument)); /* %s instead of %c and drop [i]. */
			/* Next arg. */
			currentArgument++; } printf("\n");)

	init_conf();
	while (!parse_arguments(argc, argv))
		usage();
	print_conf();

	char buffer[MAX_LINE]; //todo commandline arg?

	cout << "!!!Find The KMER!!!" << endl;

	if (config.k < 0) {
		cout << config.k
				<< " is not a valid value for k. Please select a number greater than zero"
				<< endl;
		exit(1);
	}

	//file operations.
	char* mode = "r";

	//variables for the nodes
	node_t * headNode = NULL;
	int sizeOfArray = 28;

	//to let the user know of our progress.
	unsigned long long sequenceNumber = 0;

	if ((config.sequence_file_pointer = fopen(config.sequence_file, "r")) != NULL) {
		printf("Sequence file opened properly\n");
	}else{
		printf("Sequence file failed to open");
		exit(1);
	}
	if ((config.out_file_pointer = fopen(config.out_file, "w")) != NULL) {
		printf("Out file opened properly\n");
	}else{
		printf("Out file failed to open");
		exit(1);
	}

	//Initializing array to test code. this will come from the pre processed line
	int *array = (int*) allocate_array(sizeOfArray, sizeof(int));
	int* currentArrayPosition = array;
	srand(time(NULL)); //seeding the random function.
	for (int i = 0; i < sizeOfArray; i++) {
		*currentArrayPosition = (rand() % 4);
		//printf(" currentArrayPosition %d has == %d\n", i,	*currentArrayPosition);
		currentArrayPosition++;
	}
	currentArrayPosition = array;
	//end test code.

	cout << "!!!Found " << ++sequenceNumber << " valid sequences so far!!!" << endl;

	//traverse the integer array and create a tree.
	int stopPoint = sizeOfArray - config.k + 1;
	for (int i = 0; i < stopPoint; i++) {

		DEBUG(
				printf("Extracting sequence: "); int z = 0; while (z < config.k) {printf("%d ",*(currentArrayPosition + z) ); z++;}printf("\n"););

		headNode = tree_create(headNode, currentArrayPosition++, config.k + 1);
	}

	cout << "!!!All Sequences Found, now creating histogram!!!" << endl;

	//create a temporary array for the recursive function to keep as scratch memory to hold the sequence.
	int* histogram_temp = (int*) allocate_array(config.k, sizeof(int));
	//Initialize memory.
	for (int i = 0; i < config.k; i++) {
		*(histogram_temp + i) = -1;
	}

	histo_and_free_recursive(headNode, histogram_temp, &config.k, 0);

	printf("\n");
	cout << "!!!histogram creation finished!!!" << endl;

	deallocate_array((void**) &histogram_temp);
	deallocate_array((void**) &array);

	if (fclose(config.sequence_file_pointer) == EOF) {
		printf("Sequence file close error! This is likely ok though.\n");
	}
	if (fclose(config.out_file_pointer) == EOF) {
		printf(
				"Out file close error! This is not expected and might mean the data was not written to the file properly before the close.\n");
	}

	DEBUG(printf("End of program\n"));
	return 0;
}
