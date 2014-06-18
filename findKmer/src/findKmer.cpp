//============================================================================
// Name        : findKmer.cpp
// Author      : Kalen Brown and Gus Thies
// Version     :
// Copyright   : Do not copy
// Description : This code is being scrapped and replaced by Becchi's code.
//============================================================================

#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h> //for strcmp(string1,string2) string comparison returns a 0 if they are the same.
#include <stdlib.h> //malloc is in this.
#define MAX_LINE 1001
#define DEBUG(x) x

using namespace std;

// Data structure for a tree.
struct node {
	char base;
	struct node *nextBasePtr[4];
	unsigned int counter;
};

int char2int(char base);

int main(int argc, char *argv[]) {
	FILE *fp = NULL;
	char fileName[] = "homo_sapiensupstream.fas";
	char buffer[MAX_LINE]; //todo commandline arg
	char base;
	DEBUG(int counter = 0);
	cout << "!!!Find The KMER!!!" << endl; // prints !!!Hello World!!!

	if ((fp = fopen(fileName, "r")) == NULL) {
		printf("File could not be opened. Ending program.\n");
		return (-1);
	}
	while (fgets(buffer, MAX_LINE, fp) != NULL DEBUG(&& counter != 5)) {
		DEBUG(printf("%s", buffer));
		DEBUG(counter++);
	}


	char2int(base);

	fclose(fp);
	DEBUG(printf("End of program\n"));
	return 0;
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
