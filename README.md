# peaks2utr: a robust, parallelized Python CLI for annotating 3' UTR
![CI](https://github.com/haessar/peaks2utr/actions/workflows/ci.yml/badge.svg?branch=master)
[![PYPI - Version](https://img.shields.io/pypi/v/peaks2utr.svg)](https://pypi.org/project/peaks2utr/)
[![PYPI - Python Version](https://img.shields.io/pypi/pyversions/peaks2utr.svg)](https://pypi.org/project/peaks2utr/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/428231055.svg)](https://zenodo.org/badge/latestdoi/428231055) 

<img width="400" src="https://user-images.githubusercontent.com/11962461/172829916-c2fa81e6-7ae5-4c9a-a758-c3ba4c4198cb.png">

peaks2utr is a Python command-line tool that annotates 3' untranslated regions (UTR) for a given set of aligned sequencing reads in BAM format, and canonical annotation in GFF or GTF format. peaks2utr uses MACS (https://pypi.org/project/MACS2/) to call broad "peaks" of significant read coverage in the BAM file, and uses those peaks that pass a set of criteria as a basis to annotate novel 3' UTRs. This favours BAM files from the likes of 10x Chromium runs, where signal is inherently concentrated at the distal ends of the 3' or 5' UTRs. Reads containing soft-clipped bases and polyA-tails of a given length are detected, and their end bases tallied as "truncation points". When piled up, each co-occurring truncation point is used to determine the precise end base of a given UTR. peaks2utr can be tuned to extend, override or ignore any pre-existing 3' UTR annotations in the input GFF file.

## Installation
Install latest release with:
```
pip install peaks2utr
```
Alternatively, to install from source:
```
git clone https://github.com/haessar/peaks2utr.git
cd peaks2utr
python3 -m build
python3 -m pip install dist/*.tar.gz
```
### Dependencies
Installation instructions assume a Debian / Ubuntu system with root privileges. Follow the links for instructions for other systems.
#### Required
[bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
```
apt-get install bedtools
```
#### Optional
[GenomeTools](https://github.com/genometools/genometools#building-and-installation) (for post-processing of output gff3)
```
apt-get install genometools
```
## Quick start
To check that peaks2utr has installed correctly, simply run the following in your terminal to initiate a short run with default parameters
```
peaks2utr-demo
```
This uses a small demo set of input files contained in the repository: <a href="https://github.com/haessar/peaks2utr/blob/master/peaks2utr/demo/Tb927_01_v5.1.gff" target="_blank" >Tb927_01_v5.1.gff</a> & <a href="https://github.com/haessar/peaks2utr/blob/master/peaks2utr/demo/Tb927_01_v5.1.slice.bam" target="_blank" >Tb927_01_v5.1.slice.bam</a>. When complete, you should see a file `Tb927_01_v5.1.new.gff` which contains original annotations as well as 3' UTRs with source "peaks2utr".
