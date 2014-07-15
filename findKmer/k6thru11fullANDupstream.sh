# Written on Ubuntu 14.04 LTS bash scripting by Kalen Brown using MANY google searches!

# This will launch a process for each k in the for loop above 
# Running these as background may cause harddisk thrashing as many different location reads of files will occur.
# >& /dev/null drops the findKmer output. 
# Wrapping in () is necessary otherwise the & is consumed by the command line arguments of findKmer. 
# Nice is applied before the ./ to keep the CPU priority LOW
# Nice is required to launch all 14 "tasks" or it would shutdown the system.
# Performing this script uses less than a gig of ram. possibly half a gig. 
# We use -q argument to suppress any user input requirements and minimize printing to the user. 
for count in {6..11}; do echo "starting background run for k = $count for homo_sapiensupsream.fas"; 
((nice ./Debug/findKmer -q 1 -k $count -z 100 -p homo_sapiensupstream.fas >& /dev/null)&);
done

for count in {6..11}; do echo "starting background  run for k = $count for Full_homosapiens.fa"; 
((nice ./Debug/findKmer -q 1 -k $count -z 1000 -p Full_homo_sapiens.fa >& /dev/null)&); 
done

#pkill findKmer
echo "type \"pkill -f findKmer\" to end running processes";


command watch -n 5 -t top -b -n 1 -p$(pgrep findKmer | head -20 | tr "\\n" "," | sed 's/,$//')
