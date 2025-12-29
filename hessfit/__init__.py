#from .main import main
#
## Call main function when package is imported
#main()

#__all__ = ['main']
#
#from .main import main

from os.path import dirname, basename, isfile, join
import glob
from .hessfit import main

modules = glob.glob(join(dirname(__file__), "*.py"))
__all__ = [ basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]
