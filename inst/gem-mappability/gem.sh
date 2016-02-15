#!/usr/bin/sh

#from: http://wiki.bits.vib.be/index.php/Create_a_mappability_track

reference="ce10.fa"
idxpref="ce10_index"
thr=12; # use 8 cores

gem-indexer -T ${thr} -c dna -i ${reference} -o ${idxpref}

pref="ce10_mappability"
kmer=50

# compute mappability data
gem-mappability -T ${thr} -I ${idxpref}.gem -l ${kmer} -o ${pref}_${kmer}
 
# convert results to wig and bigwig
gem-2-wig -I ${idxpref}.gem -i ${pref}_${kmer}.mappability -o ${pref}_${kmer}
./seqLengthsFromFa ${reference} > ${pref}.sizes
wigToBigWig ${pref}_${kmer}.wig ${pref}.sizes ${pref}_${kmer}.bw


# In loop
for kmer in $(seq 50 50 250); do
 
  # compute mappability data
  gem-mappability -T ${thr} -I ${idxpref}.gem -l ${kmer} -o ${pref}_${kmer}
 
  # convert results to wig and bigwig
  gem-2-wig -I ${idxpref}.gem -i ${pref}_${kmer}.mappability -o ${pref}_${kmer}
  wigToBigWig ${pref}_${kmer}.wig ${pref}.sizes ${pref}_${kmer}.bw
 
done
