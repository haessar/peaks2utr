# 3pUTR_annotation
## Quick Start
### 1. Install dependencies

Execute these commands as root
```
apt update
apt -y upgrade
apt install -y --no-install-recommends build-essential libssl-dev libffi-dev python3.8 python3-dev python3-pip git samtools genometools wget file
```

Install docopts
```
git clone https://github.com/docopt/docopts.git
cd docopts
./get_docopts.sh && \
cp docopts docopts.sh /usr/local/bin
cd ..
```

Install python packages
```
pip3 install build gffutils MACS2
```

### 2. Install 3pUTR_annotation code
```
git clone git@github.com:haessar/3pUTR_annotation.git
cd 3pUTR_annotation
python3 -m build
pip install dist/3pUTR_annotation-<version>.tar.gz
```
