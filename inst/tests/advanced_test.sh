C=GSM733780_hg19_wgEncodeBroadHistoneK562ControlStdAlnRep1

REF=hg19.2bit
MAP=MappabilityTrackHg19.gff

beads sampleGenomeGC $REF -length 200 -step 50 > $REF.plusGC.gff

#CONTROL PREP
samtools view -q 10 $C.bam > $C.sam
beads sam2gff $C.sam > $C.gff
beads extend -threePrime 164 $C.gff > $C.ext.gff
beads getGC $REF $C.ext.gff > $C.ext.plusGC.gff
beads gcHist $C.ext.plusGC.gff > $C.bg.gcHist.dat
beads gcHist $REF.plusGC.gff > genome.bg.gcHist.dat

beads gcWeigh $C.ext.plusGC.gff $C.bg.gcHist.dat genome.bg.gcHist.dat > $C.gcw.gff
beads tagCount -base 50 $C.gcw.gff > $C.gcw.binned.50bp.gff
cat $C.gcw.binned.50bp.gff | sed 's/^/chr/g' > $C.gcw.binned.50bp.gff.fix
beads mapCorr $MAP $C.gcw.binned.50bp.gff.fix -maxMap 400 > $C.gcw-map.binned.50bp.gff

#SAMPLE PREP

I=GSM733680_hg19_wgEncodeBroadHistoneK562H3k4me3StdAlnRep2
beads bed2gff $I.bed > $I.ER.gff
ER=GSM733680_hg19_wgEncodeBroadHistoneK562H3k4me3StdAlnRep2.ER.gff
cat $ER | sed 's/*/\./g' > $ER.fix

samtools view -q 10 $I.bam > $I.sam
beads sam2gff $I.sam > $I.gff
beads extend -threePrime 164 $I.gff > $I.ext.gff
beads getGC $REF $I.ext.gff > $I.ext.plusGC.gff

beads mask $ER.fix $I.ext.plusGC.gff > $I.bg.ext.plusGC.gff
beads gcHist $I.bg.ext.plusGC.gff > $I.bg.gcHist.dat
beads mask $ER.fix $REF.plusGC.gff > $I.genome.bg.plusGC.gff
beads gcHist $I.genome.bg.plusGC.gff > $I.genome.bg.gcHist.dat

beads gcWeigh $I.ext.plusGC.gff $I.bg.gcHist.dat $I.genome.bg.gcHist.dat > $I.gcw.gff
beads tagCount -base 50 $I.gcw.gff > $I.gcw.binned.50bp.gff
cat $I.gcw.binned.50bp.gff | sed 's/^/chr/g' > $I.gcw.binned.50bp.gff.fix
beads mapCorr $MAP $I.gcw.binned.50bp.gff.fix -maxMap 400 > $I.gcw-map.binned.50bp.gff

#DIV
beads divide $I.gcw-map.binned.50bp.gff $C.gcw-map.binned.50bp.gff > H3K4me3_rep1.gcw-map-div.binned.50bp.gff