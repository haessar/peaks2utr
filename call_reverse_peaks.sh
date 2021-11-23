if [ ! -f "$1.reverse.bam" ]; then
  samtools view --threads 6 -b -f 16 "$1.bam" > "$1.reverse.bam"
fi

if [ ! -f reverse_peaks.broadPeak ]; then
  macs2 callpeak -t "$1.reverse.bam" -n reverse --nomodel --extsize 100 --broad
fi

