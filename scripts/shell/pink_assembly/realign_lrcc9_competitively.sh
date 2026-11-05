# Inputs
BAM=m84185_240406_063605_s4.hifi_reads.filt.aligned.sort.lrrc9.bam
HAP1=Pink23E_4.bp.hap1.p_ctg.fa
HAP2=Pink23E_4.bp.hap2.p_ctg.fa
REF=h1tg000509l_h2tg001052l.fa
OUT=lrrc9.vs_h1tg000509l_h2tg001052l

# 1. Build the two-contig reference
samtools faidx $HAP1 h1tg000509l  > $REF
samtools faidx $HAP2 h2tg001052l >> $REF
samtools faidx $REF

# 2. BAM -> FASTQ -> minimap2 -> sorted BAM
#    -F 0x900 drops secondary/supplementary records so each read is only emitted once
samtools fastq -@ 4 -F 0x900 $BAM \
  | minimap2 -ax map-hifi -t 16 \
      -R '@RG\tID:m84185_s4\tSM:Pink23E_4\tPL:PACBIO' \
      $REF - \
  | samtools sort -@ 8 -o ${OUT}.bam -

# 3. Index and check
samtools index ${OUT}.bam
samtools flagstat ${OUT}.bam
samtools idxstats ${OUT}.bam

ALN=lrrc9.vs_h1tg000509l_h2tg001052l.bam   # new minimap2 alignment from before
BAM=m84185_240406_063605_s4.hifi_reads.filt.aligned.sort.lrrc9.bam

# 1. Read names whose primary alignment is on each contig
#    -F 0x904 drops unmapped, secondary, and supplementary records
samtools view -F 0x904 $ALN h1tg000509l | cut -f1 | sort -u > h1.names.txt
samtools view -F 0x904 $ALN h2tg001052l | cut -f1 | sort -u > h2.names.txt

# 2. Pull those reads out of the original BAM
samtools view -@ 8 -b -N h1.names.txt -o ${BAM%.bam}.h1tg000509l.bam $BAM
samtools view -@ 8 -b -N h2.names.txt -o ${BAM%.bam}.h2tg001052l.bam $BAM

# 3. Index and sanity check
for f in ${BAM%.bam}.h1tg000509l.bam ${BAM%.bam}.h2tg001052l.bam; do
  samtools index $f
  samtools view -c $f
done
wc -l h1.names.txt h2.names.txt
comm -12 h1.names.txt h2.names.txt | wc -l
