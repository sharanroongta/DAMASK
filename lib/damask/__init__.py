# $Id$
import sys

from .environment import Environment      # only one class
from .asciitable  import ASCIItable       # only one class
from .config      import Material         # will be extended to debug and numerics
#from .block       import Block            # only one class
from .result      import Result           # one class with subclasses
from .geometry    import Geometry         # one class with subclasses
from .solver      import Solver           # one class with subclasses
from .test        import Test
from util import *

try:
  from .          import core
except ImportError:
  sys.stderr.write('\nWARNING: Core module (Fortran code) not available, try to run setup_processing.py\nError Message when importing core.so: \n\n')
  core = None # from http://www.python.org/dev/peps/pep-0008/
  
  
