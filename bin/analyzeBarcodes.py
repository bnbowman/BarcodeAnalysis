#!/usr/bin/env python

import sys
import cProfile
from BarcodeAnalysis.main import main

if __name__ == '__main__':
    sys.exit(main())
    #cProfile.run('main()')
