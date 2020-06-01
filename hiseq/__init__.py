import os
import sys
import subprocess
import tempfile
import logging
import shutil
from xopen import xopen
from .hiseq import Hiseq
from .trim import trimmer
from .qc import fastqc
from .align import alignment
from .utils import helper
from .utils import argsParser
from .peak import call_peak