import subprocess
import tempfile
import logging
from shutil import which
from xopen import xopen
from .main.hiseq import Hiseq
from .trim import trimmer
from .qc import fastqc
from .align import align
from .utils import utils
from .utils import argsParser
