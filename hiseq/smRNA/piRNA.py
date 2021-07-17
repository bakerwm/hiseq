#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
piRNA analysis pipeline

Analysis in details for piRNA

Example:
1. prepare data: trim 3'adapter, remove 5'UMI,3'barcode, remove PCR dup, (...stat...)

1. remove structural RNAs (uniqe + multi)
2. remove miRNAs (unique + multi) 
3. remove reads not in [23-29] nt 
4. collapse reads ?
5. split into 4 groups: (1U+/-,10A+/-;)
6. map to TE consensus, (unique + multi), only 1U_not_10A piwi?!
7. map to genome (unique + multi) 
8. Direct map to genome

Functions:
rename fastq reads: piR0000001-0000001
piR (piRNA), (piRNA number) {reads number}
mapping, unique + multiple/unique

"""


import os
import pathlib
import argparse
from hiseq.atac.atac_rd import AtacRd
from hiseq.atac.atac_rx import AtacRx
from hiseq.utils.utils import log, update_obj, Config, get_date
from hiseq.utils.file import file_abspath




















"""
pass

"""