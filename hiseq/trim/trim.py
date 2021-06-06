#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Trim adapters, low quality bases, using cutadapt: main:=> TrimRn

1. Trim adapter,
2. remove dup
3. trim n-bases from read

3.1 TruSeq (NSR)  
    - cut 7, -7; from both ends of read (subseq)    

3.2 TruSeq (iCLIP) 
    - cut 9; from read1    

3.3 TruSeq (eCLIP)  
    - cut 10; -7 from read1  
    - cut 7, -10 from read2

optional
1. --rm-untrim, --save-too-short, ...
"""

from hiseq.utils.utils import update_obj
from hiseq.trim.trim_rn import TrimRn, get_args

class Trim(object):
    # a copy of TrimRn
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.rn = TrimRn(**self.__dict__)
        self = update_obj(self, self.rn.__dict__, force=True)


    def run(self):
#         TrimRn(**self.__dict__).run()
        self.rn.run()


def main():
    args = vars(get_args().parse_args())
    Trim(**args).run()


if __name__ == '__main__':
    main()

#