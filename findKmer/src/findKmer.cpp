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
//============================================================================
using namespace std;
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h> //for strcmp(string1,string2) string comparison returns a 0 if they are the same.
#include <stdlib.h> //malloc is in this.
#include "stdinc.h" //from Dr. Becchi's code to help reuse some functions.
// basic file operations
#include <fstream>

// Data structure for a tree.
struct node_t {
	int base;
	struct node_t *nextNodePtr[4];
	unsigned int counter;
};

//enumeration of possible bases in alphabetical order.
enum Base {
	A = 0, C = 1, G = 2, T = 3
};
#define MAX_LINE 1001
#define DEBUG(x) x

void *allocate_array(int size, size_t element_size) {
	void *mem = malloc(size * element_size);
	if (!mem) {
		printf("allocate_array():: memory allocation failed");
		exit(1);
	}
	return mem;
}

void *reallocate_array(void *array, int size, size_t element_size) {
	void *new_array = realloc(array, element_size * size);
	if (!new_array) {
		printf("reallocate_array():: memory reallocation failed");
		exit(1);
	}

	return new_array;
}
/* check that the given file can be read/written */
void check_file(char *filename, char *mode) {
	FILE *file = fopen(filename, mode);
	if (file == NULL) {
		printf("Unable to open file %s in %s mode", filename, mode);
		exit(1);
	} else
		fclose(file);
}

void file_reader(FILE *file, int k) {
	rewind(file);

	//TODO find out how many characters are on the line, then allocate memory for it.

	char *buffer = allocate_char_array(MAX_LINE);
	int i = 0;
	int j = 0;
	unsigned int c = fgetc(file);

	DEBUG(int counter = 0);
	int currentLocation = 0;
	while (NULL != fgets(buffer, MAX_LINE, file)DEBUG(&& counter != 5)) {
		DEBUG(printf("%s", buffer));
		DEBUG(counter++);

		buffer[currentLocation];
	}
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
node_t* node_branch_enter(node_t* node, int base) {
	DEBUG(printf("..node_branch_enter for base %d\n", base));
	if (node->nextNodePtr[base] == NULL) {
		DEBUG(printf("***Creating Node.\n"));
		node->nextNodePtr[base] = node_create(base);
	} else {
		node->nextNodePtr[base]->counter++;
		DEBUG(
				printf("+++Incrementing counter to %d.\n",
						node->nextNodePtr[base]->counter));
	}
	DEBUG(printf("..returning next base pointer.\n"));
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
			DEBUG(printf("-Creating the head of the tree!\n"));
			head = node_create('H');
		}
		node_t *currentNode = head;
		for (int i = 0; i < k; i++) {
			DEBUG(printf("-Moving into a branch on depth %d\n", i));
			currentNode = node_branch_enter(currentNode, array[i]);
		}
	}
	return head;
}

void histo_and_free(node_t* head, int* array, int *k) {

	if (head == NULL) {
		//DEBUG(printf("histo_and_free::Head == NULL\n"));
	} else if (head->nextNodePtr[0] == NULL && head->nextNodePtr[1] == NULL
			&& head->nextNodePtr[2] == NULL && head->nextNodePtr[3] == NULL) {
		printf(", %d\n", head->counter);
	} else {
		if (head->base == 'H') {
			printf("Head found.\n\n");
		} else {
			printf("%d", head->base);
		}
		histo_and_free(head->nextNodePtr[0],array,k);
		histo_and_free(head->nextNodePtr[1],array,k);
		histo_and_free(head->nextNodePtr[2],array,k);
		histo_and_free(head->nextNodePtr[3],array,k);
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

//close the function
	return (integer);
}

int main(int argc, char *argv[]) {
	FILE *fp = NULL;

//command line arguments go here.
	int k = 5;
	char fileName[] = "homo_sapiensupstream.fas";
	char buffer[MAX_LINE]; //todo commandline arg

	char base;
	char* mode = "r";

	cout << "!!!Find The KMER!!!" << endl; // prints !!!Hello World!!!

	if (fileName != NULL)
		check_file(fileName, mode);

	if ((fp = fopen(fileName, "r")) != NULL) {
		DEBUG(printf("File opened properly\n"));
	}

	//variables for the nodes
	node_t * headNode = NULL;
	node_t * currentNode = NULL;
	int sizeOfArray = 28;

	//Initializing array to test code. this will come from the pre processed line
	int *array = (int*) allocate_array(sizeOfArray, sizeof(int));
	int* currentArrayPosition = array;
	srand(time(NULL)); //seeding the random function.
	for (int i = 0; i < sizeOfArray; i++) {
		*currentArrayPosition = (rand()% 4  );
		printf(" currentArrayPosition %d has == %d\n", i,
				*currentArrayPosition);
		currentArrayPosition++;
	}
	currentArrayPosition = array;
	//end test code.

	//traverse the integer array and create a tree.
	int stopPoint = sizeOfArray - k + 1;
	for (int i = 0; i < stopPoint; i++) {
		DEBUG(
				printf("Looking at sequence: "); int z = 0; while (z < k) { printf("%d ",*(currentArrayPosition + z) ); z++;} printf("\n"););
		headNode = tree_create(headNode, currentArrayPosition++, k+1);
	}

	int* histogram_temp = (int*) allocate_array(k,sizeof(int));
	for(int i = 0; i < k; i++){
		*(histogram_temp+i)=0;
	}
	histo_and_free(headNode,histogram_temp,  &k);
	//file_reader(fp, k);

	//char2int(base);

	fclose(fp);
	DEBUG(printf("End of line\n"));
	return 0;
}
