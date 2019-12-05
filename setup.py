
import setuptools


def readme():
    with open('README.md', 'r') as fh:
        return fh.read()

def license():
    with open('LICENSE', 'r') as f:
        return f.read()    

setuptools.setup(
    name="hiseq", # Replace with your own username
    version="0.0.1",
    author="Ming Wang",
    author_email="wangm08@hotmail.com",
    description="A collection of tools for high-throughput sequencing data analysis",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/bakerwm/hiseq",
    license=license(),
    keywords='hiseq RNAseq ATACseq ChIPseq smRNAseq',
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': ['hiseq=hiseq.hiseq:main'],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    include_package_data=True,
    zip_safe=False,
)
