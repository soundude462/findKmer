# Written on Ubuntu 14.04 LTS bash scripting by Kalen Brown using MANY google searches!fa

command watch -n 5 -t top -b -n 1 -p pgrep findKmer | head -20 | tr "\\n" "," | sed 's/,$//'
