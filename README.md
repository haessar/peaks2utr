# peaks2utr: a robust, parallelized Python CLI for annotating 3' UTR
<img width="400" src="https://user-images.githubusercontent.com/11962461/172829916-c2fa81e6-7ae5-4c9a-a758-c3ba4c4198cb.png">
peaks2utr is a Python command-line tool that annotates 3' untranslated regions (UTR) for a given set of aligned sequencing reads in BAM format, and canonical annotation in GFF or GTF format. peaks2utr uses MACS3 (https://pypi.org/project/MACS3/) to call broad "peaks" of significant read coverage in the BAM file, and uses those peaks that pass a set of criteria as a basis to annotate novel 3' UTRs. This favours BAM files from the likes of 10x Chromium runs, where signal is inherently concentrated at the distal ends of the 3' or 5' UTRs. Reads containing soft-clipped bases and polyA-tails of a given length are detected, and their end bases tallied as "truncation points". When piled up, each co-occurring truncation point is used to determine the precise end base of a given UTR. peaks2utr can be tuned to extend, override or ignore any pre-existing 3' UTR annotations in the input GFF file.

## Installation
Install latest release with:
```
pip install peaks2utr
```
Alternatively, to install from source:
```
git clone https://github.com/haessar/peaks2utr.git
cd peaks2utr
python -m build
pip install dist/*.tar.gz
```
## Quick start
Download demo reference annotations <a href="https://github.com/haessar/peaks2utr/raw/master/demo/Tb927_01_v5.1.gff" target="_blank" >Tb927_01_v5.1.gff</a> and bam file <a href="https://github.com/haessar/peaks2utr/raw/master/demo/Tb927_01_v5.1.slice.bam" target="_blank" >Tb927_01_v5.1.slice.bam</a> and run
```
peaks2utr Tb927_01_v5.1.gff Tb927_01_v5.1.slice.bam
```
When complete, you should see a file `Tb927_01_v5.1.new.gff` which contains original annotations as well as 3' UTRs with source "peaks2utr".
