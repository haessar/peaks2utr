if [ ! -f "$1.forward.bam" ]; then
  samtools view --threads 6 -b -F 20 "$1.bam" > "$1.forward.bam"
fi

if [ ! -f forward_peaks.broadPeak ]; then
  macs2 callpeak -t "$1.forward.bam" -n forward --nomodel --extsize 100 --broad
fi
