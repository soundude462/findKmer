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
#include "CmdLineParser.h"

CmdLineParser::CmdLineParser() {
	sequence_file_name_ = NULL;
	sequence_file_pointer_ = NULL;
	out_file_name_ = NULL;
	out_file_pointer_ = NULL;
	k_ = 0;
	suppressOutputEnable_ = -1;
	zThresholdEnable_ = -1;
	zThreshold_ = -1;
	memMgr_ = new memMgt();
	configurationCopied_ = false;
	validConfiguration_ = false;
}

CmdLineParser::~CmdLineParser() {
	if (configurationCopied_ == false) {
		//free the file name strings.
		delete (sequence_file_name_);
		delete (out_file_name_);
		//free(sequence_file);
		//free(out_file);

		//If we try to use fclose on a null pointer then the program crashes.
		if (sequence_file_pointer_) {
			fclose(sequence_file_pointer_);
		}

		if (out_file_pointer_) {
			fclose(out_file_pointer_);
		}
	}
}

void CmdLineParser::usage() {
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

/*
 * Reads command line arguments.
 * Only prints to std out if there is something missing or an error
 * Does not fill in default or unknown values
 * A call to this function should be followed by set default config or a future config checker.
 */
int CmdLineParser::parse_arguments(int argc, char** argv) {
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
				out_file_name_ = argv[i];
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

				sequence_file_name_ = argv[i];
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
					k_ = k;
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
						suppressOutputEnable_ = suppressOutputEnableOption;
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
					zThresholdEnable_ = 1;
					zThreshold_ = atoi(argv[i]);
				}
			} else {
				fprintf(stderr, "Ignoring invalid option %s\n", argv[i]);
				if (suppressOutputEnable_ == 0) {
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

void CmdLineParser::print_conf(int argc) {
	fprintf(stdout, "\nATTEMPTING CONFIGURATION: \n");
	set_default_conf_();
	if (sequence_file_name_)
		fprintf(stdout, "- sequence_file file: %s\n", sequence_file_name_);
	if (out_file_name_)
		fprintf(stdout, "- export file: %s\n", out_file_name_);
	if (k_)
		fprintf(stdout, "- k size: %d\n", k_);

	fprintf(stdout, "- %s\n",
			suppressOutputEnable_ > 0 ?
					"Suppressing file read output and breaks." :
					"Showing DNA Sequence identifier and allowing breaks.");

	fprintf(stdout, "- Z score filtering is %s",
			zThresholdEnable_ ? "enabled" : "disabled");

	if (zThresholdEnable_ > 0) {
		fprintf(stdout, "\n    with threshold of %LG", zThreshold_);
	}
	fprintf(stdout, ".\n");

	//if suppressOutputEnable is false and no command line arguments have been given:
	if (suppressOutputEnable_ == 0 && argc < 2) {
		fprintf(stdout, "Press enter to proceed with this configuration.");
		getchar();
	}

	/* Double check configuration */
	if (k_ < 0 || k_ > 20) {
		fprintf(stderr,
				"%d is not a valid value for k. Please select a number greater than zero\n",
				k_);
		exit(EXIT_FAILURE);
	}

	if ((sequence_file_pointer_ = fopen(sequence_file_name_, "r")) != NULL) {
		//fprintf(stdout, "Sequence file opened properly\n");
	} else {
		fprintf(stderr, "Sequence file failed to open\n\n");
		exit(EXIT_FAILURE);
	}

	if ((out_file_pointer_ = fopen(out_file_name_, "w")) != NULL) {
		//fprintf(stdout, "Out file opened properly\n");
		fprintf(out_file_pointer_, OUT_FILE_COLUMN_HEADERS);
	} else {
		fprintf(stderr,
				"Out file failed to open\nFile MUST be in current directory.\n");
		exit(EXIT_FAILURE);
	}

	fprintf(stdout, "Sequence file and out file opened properly\n");

	fprintf(stdout, "\n");
}

int CmdLineParser::getThresholdEnable() const {
	return zThresholdEnable_;
}

/* check that the given file can be read/written */
void CmdLineParser::check_file(const char* filename, const char* mode) {
	FILE *file = fopen(filename, mode);

	if (file == NULL) {
		fprintf(stderr,
				"Unable to open file %s in %s mode\nFile MUST be in current directory.\n",
				filename, mode);
		exit(EXIT_FAILURE);
	} else
		fclose(file);
}

unsigned long int CmdLineParser::getTotalAllocatedBytes() const {
	return getTotalAllocatedBytes();
}

unsigned long int CmdLineParser::estimate_RAM_usage() {
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
		if (getSuppressOutputEnable() == 0) {
			getchar();
		}
	}

	//add one for the head node. Calculate RAM usage and Harddrive usage.
	unsigned long int maxNumberOfNodes = 1;
	double n = 1;
	while (n <= getK()) {
		maxNumberOfNodes += pow(4.0, n++);
	}
	if (((sizeof(char) * (getK() + 10)) * maxNumberOfNodes)
			>= (1024 * 1024 * 1024)) {
		cout
				<< ((sizeof(char) * (getK() + 10)) * maxNumberOfNodes)
						/ (double) (1024 * 1024 * 1024) << " gibibytes";
	} else {
		cout
				<< ((sizeof(char) * (getK() + 10)) * maxNumberOfNodes)
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
		if (getSuppressOutputEnable() == 0) {
			cout << "Hit enter to proceed or else abort the program." << endl;
			getchar();
		}
	} else {
		cout << (maxNumberOfNodes * sizeof(node_t) / (double) (1024 * 1024))
				<< " mibibytes of RAM usage likely" << endl;
	}
	return maxNumberOfNodes;
}

/*
 * Once the configuration is copied, the destructor will no longer free memory or close the files.
 * Configuration can only be copied after set_default_conf has been called (may be via print_conf.)
 * Returns EXIT_SUCCESS on success
 */
int CmdLineParser::getConfiguration(configuration *config) {
	if (validConfiguration_ == true) {
		config->k = k_;
		config->out_file_name = out_file_name_;
		config->sequence_file_name = sequence_file_name_;
		config->sequence_file_pointer = sequence_file_pointer_;
		config->zThreshold = zThreshold_;
		config->zThresholdEnable = zThresholdEnable_;
		configurationCopied_ = true;
		return EXIT_SUCCESS;
	} else {
		return EXIT_FAILURE;
	}
}

bool CmdLineParser::isValidConfiguration() const {
	return validConfiguration_;
}

bool CmdLineParser::isConfigurationCopied() const {
	return configurationCopied_;
}



/*
 * Sets anything that has not been set yet using the macros in the header file.
 * This sets the configuration as valid.
 */
void CmdLineParser::set_default_conf_() {

	if (!sequence_file_name_) {
		//here we transfer the string literal properly.
		const char* tempFileName = DEFAULT_SEQUENCE_FILE_NAME;
		sequence_file_name_ = tempFileName;
	}

	if (!k_) {
		k_ = DEFAULT_K_VALUE;
	}

	//double check default and user defined K value.
	if (!k_) {
		fprintf(stdout, "k must be greater than zero. Ending program.\n\n");
		exit(EXIT_FAILURE);
	}
	//if lessthan zero, then the user did not specify.
	if (suppressOutputEnable_ < 0) {
		suppressOutputEnable_ = DEFAULT_SUPPRESS_OUTPUT_VALUE;
	}

	if (zThresholdEnable_ < 0) {
		zThresholdEnable_ = DEFAULT_Z_THRESHOLD_ENABLE;
		zThreshold_ = DEFAULT_Z_THRESHOLD;
	}

	/*
	 * strlen is used here incase of a change in the variables, increases flexibility.
	 */
	if (!out_file_name_) {
		const char* nameOfFile = "mer_Historam_Of_";
		const char* outFileExension = ".csv";
		const char* zScoreFiltered = "zScoreFiltered";
		int stringLength = strlen("999") + strlen(nameOfFile)
				+ strlen(sequence_file_name_) + strlen(outFileExension);

		if (zThresholdEnable_ == 0) {

			memMgr_->allocateArray(sizeof(char), stringLength, out_file_name_);

			sprintf(out_file_name_, "%d%s%s%s", k_, nameOfFile,
					sequence_file_name_, outFileExension);
		} else {

			memMgr_->allocateArray(sizeof(char),
					stringLength + strlen(zScoreFiltered), out_file_name_);

			sprintf(out_file_name_, "%d%s%s%s%s", k_, nameOfFile,
					sequence_file_name_, zScoreFiltered, outFileExension);
		}
	}

	//This verified that the configuration is finally valid.
	validConfiguration_ = true;
}

//Getter
int CmdLineParser::getK() const {
	return k_;
}

//Getter
const FILE* CmdLineParser::getOutFilePointer() const {
	return out_file_pointer_;
}

//Getter
const FILE* CmdLineParser::getSequenceFilePointer() const {
	return sequence_file_pointer_;
}

//Getter
int CmdLineParser::getSuppressOutputEnable() const {
	return suppressOutputEnable_;
}

//Getter
long double CmdLineParser::getThreshold() const {
	return zThreshold_;
}
//Getter
const char* CmdLineParser::getSequenceFileName() const {
	return sequence_file_name_;
}
//Getter
char* CmdLineParser::getOutFileName() const {
	return out_file_name_;
}
