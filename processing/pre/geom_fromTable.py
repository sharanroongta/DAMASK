#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,re,string,math
import scipy.spatial, numpy as np
from optparse import OptionParser
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]



# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Generate geometry description and material configuration from position, phase, and orientation (or microstructure) data.

""", version = scriptID)

parser.add_option('--coordinates',
                  dest = 'coordinates',
                  type = 'string', metavar = 'string',
                  help = 'coordinates label')
parser.add_option('--phase',
                  dest = 'phase',
                  type = 'string', metavar = 'string',
                  help = 'phase label')
parser.add_option('--microstructure',
                  dest = 'microstructure',
                  type = 'string', metavar = 'string',
                  help = 'microstructure label')
parser.add_option('-t', '--tolerance',
                  dest = 'tolerance',
                  type = 'float', metavar = 'float',
                  help = 'angular tolerance for orientation squashing [%default]')
parser.add_option('-e', '--eulers',
                  dest = 'eulers',
                  type = 'string', metavar = 'string',
                  help = 'Euler angles label')
parser.add_option('-d', '--degrees',
                  dest = 'degrees',
                  action = 'store_true',
                  help = 'angles are given in degrees [%default]')
parser.add_option('-m', '--matrix',
                  dest = 'matrix',
                  type = 'string', metavar = 'string',
                  help = 'orientation matrix label')
parser.add_option('-a',
                  dest='a',
                  type = 'string', metavar = 'string',
                  help = 'crystal frame a vector label')
parser.add_option('-b',
                  dest='b',
                  type = 'string', metavar = 'string',
                  help = 'crystal frame b vector label')
parser.add_option('-c',
                  dest = 'c',
                  type =  'string', metavar='string',
                  help = 'crystal frame c vector label')
parser.add_option('-q', '--quaternion',
                  dest = 'quaternion',
                  type = 'string', metavar='string',
                  help = 'quaternion label')
parser.add_option('--axes',
                  dest = 'axes',
                  type = 'string', nargs = 3, metavar = ' '.join(['string']*3),
                  help = 'orientation coordinate frame in terms of position coordinate frame [same]')
parser.add_option('-s', '--symmetry',
                  dest = 'symmetry',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'crystal symmetry %default {{{}}} '.format(', '.join(damask.Symmetry.lattices[1:])))
parser.add_option('--homogenization',
                  dest = 'homogenization',
                  type = 'int', metavar = 'int',
                  help = 'homogenization index to be used [%default]')
parser.add_option('--crystallite',
                  dest = 'crystallite',
                  type = 'int', metavar = 'int',
                  help = 'crystallite index to be used [%default]')

parser.set_defaults(symmetry       = [damask.Symmetry.lattices[-1]],
                    tolerance      = 0.0,
                    degrees        = False,
                    homogenization = 1,
                    crystallite    = 1,
                   )

(options,filenames) = parser.parse_args()

input = [options.eulers     != None,
         options.a          != None and \
         options.b          != None and \
         options.c          != None,
         options.matrix     != None,
         options.quaternion != None,
         options.microstructure != None,
        ]

if np.sum(input) != 1:
  parser.error('need either microstructure label or exactly one orientation input format.')
if options.axes != None and not set(options.axes).issubset(set(['x','+x','-x','y','+y','-y','z','+z','-z'])):
  parser.error('invalid axes {} {} {}.'.format(*options.axes))

(label,dim,inputtype) = [(options.eulers,3,'eulers'),
                         ([options.a,options.b,options.c],[3,3,3],'frame'),
                         (options.matrix,9,'matrix'),
                         (options.quaternion,4,'quaternion'),
                         (options.microstructure,1,'microstructure'),
                        ][np.where(input)[0][0]]                                                    # select input label that was requested
toRadians = math.pi/180.0 if options.degrees else 1.0                                               # rescale degrees to radians
options.tolerance *= toRadians                                                                      # ensure angular tolerance in radians

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:
    table = damask.ASCIItable(name = name,
                              outname = os.path.splitext(name)[0]+'.geom' if name else name,
                              buffered = False)
  except: continue
  table.croak('\033[1m'+scriptName+'\033[0m'+(': '+name if name else ''))

# ------------------------------------------ read head ---------------------------------------  

  table.head_read()                                                                                 # read ASCII header info

# ------------------------------------------ sanity checks ---------------------------------------  

  coordDim = table.label_dimension(options.coordinates)

  errors = []
  if not 3 >= coordDim >= 2:
    errors.append('coordinates {} need to have two or three dimensions.'.format(options.coordinates))
  if not np.all(table.label_dimension(label) == dim):
    errors.append('input {} needs to have dimension {}.'.format(label,dim))
  if options.phase != None and table.label_dimension(options.phase) != 1:
    errors.append('phase column {} is not scalar.'.format(options.phase))
  
  if errors  != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue

  table.data_readArray([options.coordinates,label]+([] if options.phase == None else [options.phase]))
  
  if coordDim == 2:
    table.data = np.insert(table.data,2,np.zeros(len(table.data)),axis=1)                           # add zero z coordinate for two-dimensional input
  if options.phase == None:
    table.data = np.column_stack((table.data,np.ones(len(table.data))))                             # add single phase if no phase column given

# --------------- figure out size and grid ---------------------------------------------------------

  coords = [np.unique(table.data[:,i]) for i in xrange(3)]
  mincorner = np.array(map(min,coords))
  maxcorner = np.array(map(max,coords))
  grid   = np.array(map(len,coords),'i')
  size   = grid/np.maximum(np.ones(3,'d'), grid-1.0) * (maxcorner-mincorner)                        # size from edge to edge = dim * n/(n-1) 
  size   = np.where(grid > 1, size, min(size[grid > 1]/grid[grid > 1]))                             # spacing for grid==1 equal to smallest among other spacings
  delta  = size/np.maximum(np.ones(3,'d'), grid)
  origin = mincorner - 0.5*delta                                                                    # shift from cell center to corner

  N = grid.prod()

  if  N != len(table.data):
    errors.append('data count {} does not match grid {}.'.format(len(table.data),' x '.join(map(repr,grid))))
  if   np.any(np.abs(np.log10((coords[0][1:]-coords[0][:-1])/delta[0])) > 0.01) \
    or np.any(np.abs(np.log10((coords[1][1:]-coords[1][:-1])/delta[1])) > 0.01) \
    or np.any(np.abs(np.log10((coords[2][1:]-coords[2][:-1])/delta[2])) > 0.01):
    errors.append('regular grid spacing {} violated.'.format(' x '.join(map(repr,delta))))

  if errors  != []:
    table.croak(errors)
    table.close(dismiss = True)
    continue
  
# ------------------------------------------ process data ------------------------------------------

  colOri = table.label_index(label)+(3-coordDim)                                                      # column(s) of orientation data (following 3 or 2 coordinates that were expanded to 3!)

  if inputtype == 'microstructure':
    microstructure = table.data[:,colOri]
    nGrains = len(np.unique(microstructure))
  else:
    colPhase = colOri + np.sum(dim)                                                                   # column of phase data comes after orientation
    index = np.lexsort((table.data[:,0],table.data[:,1],table.data[:,2]))                             # index of rank when sorting x fast, z slow
    rank  = np.argsort(index)                                                                         # rank of index
    KDTree = scipy.spatial.KDTree((table.data[:,:3]-mincorner) / delta)                               # build KDTree with dX = dY = dZ = 1 and origin 0,0,0
  
    microstructure = np.zeros(N,dtype = 'uint32')                                                     # initialize empty microstructure
    symQuats = []                                                                                     # empty list of sym equiv orientations
    phases   = []                                                                                     # empty list of phase info
    nGrains = 0                                                                                       # counter for detected grains
    myRank  = 0                                                                                       # rank of current grid point
    for z in xrange(grid[2]):
      for y in xrange(grid[1]):
        for x in xrange(grid[0]):
          if (myRank+1)%(N/100.) < 1: table.croak('.',False)
          myData = table.data[index[myRank]]
          mySym = options.symmetry[min(int(myData[colPhase]),len(options.symmetry))-1]                # select symmetry from option (take last specified option for all with higher index)
          if inputtype == 'eulers':
            o = damask.Orientation(Eulers = myData[colOri:colOri+3]*toRadians,
                                   symmetry = mySym).reduced()
          elif inputtype == 'matrix':
            o = damask.Orientation(matrix = myData[colOri:colOri+9].reshape(3,3).transpose(),
                                   symmetry = mySym).reduced()
          elif inputtype == 'frame':
            o = damask.Orientation(matrix = np.hstack((myData[colOri[0]:colOri[0]+3],
                                                       myData[colOri[1]:colOri[1]+3],
                                                       myData[colOri[2]:colOri[2]+3],
                                                      )).reshape(3,3),
                                   symmetry = mySym).reduced()
          elif inputtype == 'quaternion':
            o = damask.Orientation(quaternion = myData[colOri:colOri+4],
                                   symmetry = mySym).reduced()

          oInv = o.quaternion.conjugated()
          neighbors = KDTree.query_ball_point([x,y,z], 3)                                             # search points within radius
          breaker = False

          for n in neighbors:                                                                         # check each neighbor
            if myRank <= rank[n] or table.data[n,colPhase] != myData[colPhase]: continue              # skip myself, anyone further ahead (cannot yet have a grain ID), and other phases
            for q in symQuats[microstructure[rank[n]]-1]:
              if abs((q*oInv).asAngleAxis()[0]) <= options.tolerance:                                 # found existing orientation resembling me
                microstructure[myRank] = microstructure[rank[n]]
                breaker = True; break
            if breaker: break

          if microstructure[myRank] == 0:                                                             # no other orientation resembled me
            nGrains += 1                                                                              # make new grain ...
            microstructure[myRank] = nGrains                                                          # ... and assign to me
            symQuats.append(o.equivalentQuaternions())                                                # store all symmetrically equivalent orientations for future comparison
            phases.append(myData[colPhase])                                                           # store phase info for future reporting

          myRank += 1

    table.croak('')

# --- generate header ----------------------------------------------------------------------------

  info = {
          'grid':    grid,
          'size':    size,
          'origin':  origin,
          'microstructures': nGrains,
          'homogenization':  options.homogenization,
         }

  table.croak(['grid     a b c:  %s'%(' x '.join(map(str,info['grid']))),
               'size     x y z:  %s'%(' x '.join(map(str,info['size']))),
               'origin   x y z:  %s'%(' : '.join(map(str,info['origin']))),
               'homogenization:  %i'%info['homogenization'],
               'microstructures: %i'%info['microstructures'],
              ])
    
# --- write header ---------------------------------------------------------------------------------

  formatwidth = 1+int(math.log10(info['microstructures']))

  if inputtype == 'microstructure':
    config_header = []
  else:
    config_header = ['<microstructure>']
    for i,phase in enumerate(phases):
      config_header += ['[Grain%s]'%(str(i+1).zfill(formatwidth)),
                        'crystallite %i'%options.crystallite,
                        '(constituent)\tphase %i\ttexture %s\tfraction 1.0'%(phase,str(i+1).rjust(formatwidth)),
                       ]
  
    config_header += ['<texture>']
    for i,quats in enumerate(symQuats):
      config_header += ['[Grain%s]'%(str(i+1).zfill(formatwidth)),
                        'axes\t%s %s %s'%tuple(options.axes) if options.axes != None else '',
                        '(gauss)\tphi1 %g\tPhi %g\tphi2 %g\tscatter 0.0\tfraction 1.0'%tuple(np.degrees(quats[0].asEulers())),
                       ]
  
  table.labels_clear()
  table.info_clear()
  table.info_append([
    scriptID + ' ' + ' '.join(sys.argv[1:]),
    "grid\ta {}\tb {}\tc {}".format(*info['grid']),
    "size\tx {}\ty {}\tz {}".format(*info['size']),
    "origin\tx {}\ty {}\tz {}".format(*info['origin']),
    "homogenization\t{}".format(info['homogenization']),
    "microstructures\t{}".format(info['microstructures']),
    config_header,
    ])
  table.head_write()
  
# --- write microstructure information ------------------------------------------------------------

  table.data = microstructure.reshape(info['grid'][1]*info['grid'][2],info['grid'][0])
  table.data_writeArray('%%%ii'%(formatwidth),delimiter=' ')
  
#--- output finalization --------------------------------------------------------------------------

  table.close()
