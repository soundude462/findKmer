for count in {6..11}; do echo "starting run for k = $count for homo_sapiensupsream.fas"; 
time ./Debug/findKmer -q -k $count -p homo_sapiensupstream.fas; 
done
for count in {6..11}; do echo "starting run for k = $count for Full_homosapiens.fa"; 
time ./Debug/findKmer -k $count -p Full_homo_sapiens.fa; 
done

