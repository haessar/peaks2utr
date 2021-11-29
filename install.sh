sudo -s
apt update
apt -y upgrade
apt install -y --no-install-recommends build-essential libssl-dev libffi-dev python3.8 python3-dev python3.8-venv git samtools genometools wget file

wget https://bootstrap.pypa.io/get-pip.py
python3 get-pip.py

git clone https://github.com/docopt/docopts.git
cd docopts || exit
./get_docopts.sh
cp docopts docopts.sh /usr/local/bin
cd ..

pip install build gffutils MACS2

git clone https://github.com/haessar/3pUTR_annotation.git
cd 3pUTR_annotation || exit
python3 -m build
pip install dist/3pUTR_annotation-*.tar.gz
