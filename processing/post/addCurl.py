#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import numpy as np

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Add column(s) containing curl of requested column(s).
Operates on periodic ordered three-dimensional data sets of vector and tensor fields.
""", version = scriptID)

parser.add_option('-p','--pos','--periodiccellcenter',
                  dest = 'pos',
                  type = 'string', metavar = 'string',
                  help = 'label of coordinates [%default]')
parser.add_option('-l','--label',
                  dest = 'labels',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'label(s) of field values')

parser.set_defaults(pos = 'pos',
                   )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

if options.labels is None: parser.error('no data column specified.')

for name in filenames:
    damask.util.report(scriptName,name)

    table = damask.Table.load(StringIO(''.join(sys.stdin.read())) if name is None else name)
    grid,size,origin = damask.grid_filters.cell_coord0_gridSizeOrigin(table.get(options.pos))

    for label in options.labels:
        field = table.get(label)
        shape = (3,) if np.prod(field.shape)//np.prod(grid) == 3 else (3,3)                         # vector or tensor
        field = field.reshape(tuple(grid)+(-1,),order='F').reshape(tuple(grid)+shape)
        curl  = damask.grid_filters.curl(size,field)
        table = table.add('curlFFT({})'.format(label),
                          curl.reshape(tuple(grid)+(-1,)).reshape(-1,np.prod(shape),order='F'),
                          scriptID+' '+' '.join(sys.argv[1:]))

    table.save((sys.stdout if name is None else name), legacy=True)
