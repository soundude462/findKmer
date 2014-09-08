/*
 * ============================================================================
 * memMgt.cpp of findKmer project.
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

#include "MemMgt.h"

memMgt::memMgt() {
	TotalAllocatedBytes = 0;
}

memMgt::~memMgt() {
	if(TotalAllocatedBytes != 0){
		fprintf(stderr, "Memory Leak detected!");
	}
}

unsigned long int memMgt::getTotalAllocatedBytes() const {
	return TotalAllocatedBytes;
}

/*
 * returns EXIT_SUCCESS if successful. ... THIS code will probably not be able to see the number of bytes from a void pointer.
 * need to have overloaded methods that accept different pointer types to know, or have an input argument for the sizeof thing...
 */
int memMgt::deallocateArray(void** array,int size, size_t element_size) {
	int Return=EXIT_SUCCESS;
	if ((*array)!=NULL) {
		//TODO I am unsure of this, math with size_t may not be stable? unsure of size_t type/use...
		this->TotalAllocatedBytes -= (size* element_size);
		DEBUG(printf("Freeing %d bytes of memory",(size* element_size)));
		free(*array);
		*array = NULL;
	}else{
		Return = 1;
	}
	return Return;
}

/*
 * returns EXIT_SUCCESS if successful
 */
int memMgt::allocateArray( size_t elementSize, int numberOfElements, void* destinationPointer) {
	int Return = EXIT_SUCCESS;
	destinationPointer = malloc(numberOfElements * elementSize);
	if (!destinationPointer) {
		fprintf(stderr, "allocate_array():: memory allocation failed\n");
		exit(EXIT_FAILURE);
	}
	TotalAllocatedBytes += numberOfElements * elementSize;
	return Return;
}
