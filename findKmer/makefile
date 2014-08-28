
RM := rm -rf

# All Target
all: findKmer
findKmer:
# Tool invocations
	@echo 'Building target: $@'
	@echo 'Invoking: Cross G++'
	g++ ./src/findKmer.cpp -o findKmer -O3 -w
	@echo 'Finished building target: $@'
	@echo 'Goodbye from handwritten make file for findKmer:'

# Other Targets
clean:
	-$(RM) findKmer
	-@echo ' '

