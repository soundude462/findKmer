//============================================================================
// Name        : findKmer.cpp
// Author      : Kalen Brown and Gus Thies
// Version     :
// Copyright   : Do not copy
// Description : This code is being scrapped and replaced by Becchi's code.
//============================================================================

using namespace std;
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h> //for strcmp(string1,string2) string comparison returns a 0 if they are the same.
#include <stdlib.h> //malloc is in this.
#include "stdinc.h"
// basic file operations
#include <fstream>

// Data structure for a tree.
struct node {
	int base;
	struct node *nextBasePtr[4];
	unsigned int counter;
};
//enumeration of possible bases in alphabetical order.
enum Base {
	A = 0, C = 1, G = 2, T = 3
};
#define MAX_LINE 1001
#define DEBUG(x) x

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
 * Function reads in a pointer to an integer array, the size of the array, and an integer.
 * The function shifts the contents of the array to the left and inserts a new element at the tail.
 */

/*
 * Reads in a pointer to an integer array, the size of the array.
 * Manages the tree to hold the string of integers
 *
 */

/*
 * increment the counter for then odes as we step on it.
 */
void enter_the_node(node* node) {
	node->counter++;
}

/*
 * Checks to see if there is a branch in that direction.
 * return true if the branch exists, else return false
 */
bool create_branch(node* node, int base) {
	if(node->nextBasePtr[base] == NULL) {
		return false;
	}
	return true;

}
/*
 * Creates a tree node.
 * Brings in the base of the node to create
 *
 */
node* create_node(int base) {
	node* mem = NULL;
	mem = (node*) malloc(sizeof(node));
	mem->base = base;
	mem->counter = 1;
	mem->nextBasePtr[0] = NULL;
	mem->nextBasePtr[1] = NULL;
	mem->nextBasePtr[2] = NULL;
	mem->nextBasePtr[3] = NULL;
	return mem;
}
/*
 * Adds a branch
 */
void add_branch(node* stumpNode, int base) {
	stumpNode->nextBasePtr[base] = create_node(base);
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
	return integer;
}

int main(int argc, char *argv[]) {
	FILE *fp = NULL;

	//command line arguments go here.
	char fileName[] = "homo_sapiensupstream.fas";
	char buffer[MAX_LINE]; //todo commandline arg
	int k = 10;

	char base;
	char* mode = "r";

	cout << "!!!Find The KMER!!!" << endl; // prints !!!Hello World!!!

	if (fileName != NULL)
		check_file(fileName, mode);

	if ((fp = fopen(fileName, "r")) != NULL) {
		DEBUG(printf("File opened properly\n"));
	}

	file_reader(fp, k);

	char2int(base);

	fclose(fp);
	DEBUG(printf("End of program\n"));
	return 0;
}
