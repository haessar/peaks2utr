# peaks2utr: a robust, parallelized Python CLI for annotating 3' UTR
<!-- ## Quick Start -->
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
