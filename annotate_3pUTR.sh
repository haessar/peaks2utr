if [ ! -f "$1.forward.bam" ]; then
    samtools view -b -F 20 "$1.bam" > "$1.forward.bam"
fi

if [ ! -f "$1.reverse.bam" ]; then
    samtools view -b -f 16 "$1.bam" > "$1.reverse.bam"
fi

macs2 callpeak -t "$1.reverse.bam" -n reverse  --nomodel --extsize 100 --broad
macs2 callpeak -t "$1.forward.bam" -n forward  --nomodel --extsize 100 --broad

./peaks_to_UTR.py

cat Orig.PRFA01000011.gff three_prime_UTRs.gff > full.gff
gt gff3 -sort -retainids -tidy full.gff > full.sorted.gff
