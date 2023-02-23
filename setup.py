from setuptools import setup

setup(
    name="peaks2utr",
    packages=["peaks2utr"],
    install_requires=[
        "pybedtools",
        "gffutils==0.10.1",
        "pysam",
        "numpy>=1.21.4",
        "tqdm",
        "asgiref",
        "psutil",
        "six"
    ]
)
