/*
 * ============================================================================
 * cmdLineParser.cpp of findKmer project.
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
#include "cmdlineParser.h"

cmdline_parser::cmdline_parser() {
	sequence_file_name = NULL;
	sequence_file_pointer = NULL;
	out_file_name = NULL;
	out_file_pointer = NULL;
	k = 0;
	suppressOutputEnable = -1;
	zThresholdEnable = -1;
	zThreshold = -1;
	memMgr = new memMgt();
}

cmdline_parser::~cmdline_parser() {
	//free the file name strings.
	delete (sequence_file_name);
	delete (out_file_name);
	//free(sequence_file);
	//free(out_file);

	//If we try to use fclose on a null pointer then the program crashes.
	if (sequence_file_pointer) {
		fclose(sequence_file_pointer);
	}

	if (out_file_pointer) {
		fclose(out_file_pointer);
	}

}

void cmdline_parser::usage() {
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
	fprintf(stdout,
			"             [--quiet|-q  < 0 for FALSE | 1 for TRUE >] \n               Suppress file read output and breaks.\n                Default is %s.\n\n",
			DEFAULT_SUPPRESS_OUTPUT_VALUE ? "true" : "false");

	long double tempzThreshold = DEFAULT_Z_THRESHOLD;
	fprintf(stdout,
			"             [--zthreshold|-z  < Threshold_for_Z >] \n               Suppress sequences with Z scores < threshold.\n                Default is %s with a value of %LG.\n\n",
			DEFAULT_Z_THRESHOLD_ENABLE ? "enabled" : "disabled",
			tempzThreshold);
	fprintf(stdout, "\n");
}

int cmdline_parser::parse_arguments(int argc, char** argv) {
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
				this->out_file_name = argv[i];
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

				this->sequence_file_name = argv[i];
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
					this->k = k;
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
						this->suppressOutputEnable = suppressOutputEnableOption;
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
					this->zThresholdEnable = 1;
					this->zThreshold = atoi(argv[i]);
				}
			} else {
				fprintf(stderr, "Ignoring invalid option %s\n", argv[i]);
				if (this->suppressOutputEnable == 0) {
					fprintf(stderr, "Press enter to continue.\n");
					getchar();
				}
			}
			i++;
		}
	}
	//will never reach this line.
	return 1;
}

void cmdline_parser::print_conf(int argc) {
	fprintf(stdout, "\nATTEMPTING CONFIGURATION: \n");
	this->set_default_conf();
	if (this->sequence_file_name)
		fprintf(stdout, "- sequence_file file: %s\n", this->sequence_file_name);
	if (this->out_file_name)
		fprintf(stdout, "- export file: %s\n", this->out_file_name);
	if (this->k)
		fprintf(stdout, "- k size: %d\n", this->k);

	fprintf(stdout, "- %s\n",
			this->suppressOutputEnable > 0 ?
					"Suppressing file read output and breaks." :
					"Showing DNA Sequence identifier and allowing breaks.");

	fprintf(stdout, "- Z score filtering is %s",
			this->zThresholdEnable ? "enabled" : "disabled");

	if (this->zThresholdEnable > 0) {
		fprintf(stdout, "\n    with threshold of %LG", this->zThreshold);
	}
	fprintf(stdout, ".\n");

	//if suppressOutputEnable is false and no command line arguments have been given:
	if (this->suppressOutputEnable == 0 && argc < 2) {
		fprintf(stdout, "Press enter to proceed with this configuration.");
		getchar();
	}

	/* Double check configuration */
	if (this->k < 0 || this->k > 20) {
		fprintf(stderr,
				"%d is not a valid value for k. Please select a number greater than zero\n",
				this->k);
		exit(EXIT_FAILURE);
	}

	if ((this->sequence_file_pointer = fopen(this->sequence_file_name, "r")) != NULL) {
		//fprintf(stdout, "Sequence file opened properly\n");
	} else {
		fprintf(stderr, "Sequence file failed to open\n\n");
		exit(EXIT_FAILURE);
	}

	if ((this->out_file_pointer = fopen(this->out_file_name, "w")) != NULL) {
		//fprintf(stdout, "Out file opened properly\n");
		fprintf(this->out_file_pointer, OUT_FILE_COLUMN_HEADERS);
	} else {
		fprintf(stderr,
				"Out file failed to open\nFile MUST be in current directory.\n");
		exit(EXIT_FAILURE);
	}

	fprintf(stdout, "Sequence file and out file opened properly\n");

	fprintf(stdout, "\n");
}

int cmdline_parser::getThresholdEnable() const {
	return zThresholdEnable;
}

/* check that the given file can be read/written */
void cmdline_parser::check_file(const char* filename, const char* mode) {
	FILE *file = fopen(filename, mode);
	if (file == NULL) {
		fprintf(stderr,
				"Unable to open file %s in %s mode\nFile MUST be in current directory.\n",
				filename, mode);
		exit(EXIT_FAILURE);
	} else
		fclose(file);
}

unsigned long int cmdline_parser::getTotalAllocatedBytes() const {
	return getTotalAllocatedBytes();
}

unsigned long int cmdline_parser::estimate_RAM_usage() {
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
		if (this->getSuppressOutputEnable() == 0) {
			getchar();
		}
	}

	//add one for the head node. Calculate RAM usage and Harddrive usage.
	unsigned long int maxNumberOfNodes = 1;
	double n = 1;
	while (n <= this->getK()) {
		maxNumberOfNodes += pow(4.0, n++);
	}
	if (((sizeof(char) * (this->getK() + 10)) * maxNumberOfNodes)
			>= (1024 * 1024 * 1024)) {
		cout
				<< ((sizeof(char) * (this->getK() + 10)) * maxNumberOfNodes)
						/ (double) (1024 * 1024 * 1024) << " gibibytes";
	} else {
		cout
				<< ((sizeof(char) * (this->getK() + 10)) * maxNumberOfNodes)
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
		if (this->getSuppressOutputEnable() == 0) {
			cout << "Hit enter to proceed or else abort the program." << endl;
			getchar();
		}
	} else {
		cout << (maxNumberOfNodes * sizeof(node_t) / (double) (1024 * 1024))
				<< " mibibytes of RAM usage likely" << endl;
	}
	return maxNumberOfNodes;
}

void cmdline_parser::set_default_conf() {

	if (!this->sequence_file_name) {
		//here we transfer the string literal properly.
		const char* tempFileName = DEFAULT_SEQUENCE_FILE_NAME;
		this->sequence_file_name = tempFileName;
	}

	if (!this->k) {
		this->k = DEFAULT_K_VALUE;
	}

	//double check default and user defined K value.
	if (!this->k) {
		fprintf(stdout, "k must be greater than zero. Ending program.\n\n");
		exit(EXIT_FAILURE);
	}
	//if lessthan zero, then the user did not specify.
	if (this->suppressOutputEnable < 0) {
		this->suppressOutputEnable = DEFAULT_SUPPRESS_OUTPUT_VALUE;
	}

	if (this->zThresholdEnable < 0) {
		this->zThresholdEnable = DEFAULT_Z_THRESHOLD_ENABLE;
		this->zThreshold = DEFAULT_Z_THRESHOLD;
	}

	if (!this->out_file_name) {
		const char* nameOfFile = "mer_Historam_Of_";
		const char* outFileExension = ".csv";
		const char* zScoreFiltered = "zScoreFiltered";
		if (this->zThresholdEnable == 0) {
			this->out_file_name = (char*) this->memMgr->allocate_array(
					strlen("999") + strlen(nameOfFile)
							+ strlen(this->sequence_file_name)
							+ strlen(outFileExension), sizeof(char));
			sprintf(this->out_file_name, "%d%s%s%s", this->k, nameOfFile,
					this->sequence_file_name, outFileExension);
		} else {
			this->out_file_name = (char*) this->memMgr->allocate_array(
					strlen("999") + strlen(nameOfFile)
							+ strlen(this->sequence_file_name)
							+ strlen(outFileExension) + strlen(zScoreFiltered),
					sizeof(char));
			sprintf(this->out_file_name, "%d%s%s%s%s", this->k, nameOfFile,
					this->sequence_file_name, zScoreFiltered, outFileExension);
		}
	}

}

//Getter
int cmdline_parser::getK() const {
	return k;
}

//Getter
const FILE* cmdline_parser::getOutFilePointer() const {
	return out_file_pointer;
}

//Getter
const FILE* cmdline_parser::getSequenceFilePointer() const {
	return sequence_file_pointer;
}

//Getter
int cmdline_parser::getSuppressOutputEnable() const {
	return suppressOutputEnable;
}

//Getter
long double cmdline_parser::getThreshold() const {
	return zThreshold;
}
//Getter
const char* cmdline_parser::getSequenceFileName() const {
	return sequence_file_name;
}
//Getter
char* cmdline_parser::getOutFileName() const {
	return out_file_name;
}
