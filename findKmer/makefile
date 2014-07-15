
RM := rm -rf

# All Target
all: findKmer
findKmer:
# Tool invocations
	@echo 'Building target: $@'
	@echo 'Invoking: Cross G++'
	g++ ./src/findKmer.cpp -o findKmer -O3 -w
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) findKmer
	-@echo ' '

