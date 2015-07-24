#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import threading,time,os,subprocess,shlex,string
import os,re,sys,math,string
import numpy as np
from optparse import OptionParser
import damask
from collections import defaultdict


scriptID   = string.replace('$Id: geom_directionalNonEquiaxedGrain.py 4290 2015-07-24 08:41:08Z hm.zhang $','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

def execute(cmd,streamIn=None,wd='./'):
  '''
    executes a command in given directory and returns stdout and stderr for optional stdin
  '''
  initialPath=os.getcwd()
  os.chdir(wd)
  process = subprocess.Popen(shlex.split(cmd),stdout=subprocess.PIPE,stderr = subprocess.PIPE,stdin=subprocess.PIPE)
  if streamIn != None:
    out,error = process.communicate(streamIn.read())
  else:
    out,error = process.communicate()
  os.chdir(initialPath)
  return out,error

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------
identifiers = {
        'grid':   ['a','b','c'],
        'size':   ['x','y','z'],
        'origin': ['x','y','z'],
          }
mappings = {
        'grid':            lambda x: int(x),
        'size':            lambda x: float(x),
        'origin':          lambda x: float(x),
        'homogenization':  lambda x: int(x),
        'microstructures': lambda x: int(x),
          }

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Generate the geometry description of a directional non-equiaxed grain structure, e.g., RVE cutted from a cold-rolling sheet.
The initial equiaxed grain structure is generated by standard Voronoi tessellation, '--reduct' specifies the thickness
reduction after rolling, and '-n' specified the number of samples cutted, e.g., n=5, then five samples will be cutted from the
sheet along \\theta = 0 (the rolling direction), \\theta = 22.5, \\theta = 45, \\theta = 67.5, and \\theta = 90 (the
transversal direction ).

""", version = scriptID)

parser.add_option('-N', dest='N', type='int', metavar='int', 
                  help='number of seed points to distribute [%default]')
parser.add_option('-r', '--rnd', dest='randomSeed', type='int', nargs = 2, metavar=' '.join(['int']*2), 
                  help='seed of random number generator [%default]')
parser.add_option('-m', '--microstructure', dest='microstructure', type='int', metavar='int',
                  help='first microstructure index [%default]')
parser.add_option('-g', '--grid', dest='grid', type='int', nargs = 3, metavar=' '.join(['int']*3),
                  help='a,b,c grid of hexahedral box [from seeds file]')
parser.add_option('-s', '--size', dest='size', type='float', nargs = 3, metavar=' '.join(['float']*3),
                  help='x,y,z size of hexahedral box [1.0 along largest grid point number]')
parser.add_option('--phase', dest='phase', type='int', metavar = 'int',
                  help='phase index to be used [%default]')
parser.add_option('--crystallite', dest='crystallite', type='int', metavar = 'int',
                  help='crystallite index to be used [%default]')
parser.add_option('-c', '--configuration', dest='config', action='store_true',
                  help='output material configuration [%default]')
parser.add_option('--secondphase', type='float', dest='secondphase', metavar= 'float',
                  help='volume fraction of randomly distribute second phase [%default]')
parser.add_option('-l', '--laguerre', dest='laguerre', action='store_true',
                  help='use Laguerre (weighted Voronoi) tessellation [%default]')
parser.add_option('-n', dest='number', type='int', metavar='int', 
                  help='the angle(degree) between the longitudinal direction of RVE and the rolling direction [%default]')
parser.add_option('--reduct', dest='reduction', type='float', metavar='float', 
                  help='thickness reduction of rolling [%default]')
parser.set_defaults(
                    N = 500,
                    grid   = (200,100,50),
                    size   = (2.0,1.0,0.5),
                    phase          = 1,
                    crystallite    = 1,
                    secondphase    = 0.0,
                    microstructure = 1,
                    laguerre = False,
                    randomSeed = (None,None),
                    config = False,
                    number = 5,
                    reduction = 0.6
                  )
(options,filenames) = parser.parse_args()
options.grid = np.array(options.grid)

sizeX = sizeY = max(options.size[0], options.size[1])
gridX, gridY = int(np.ceil(sizeX/options.size[0]*options.grid[0]))+1,int(np.ceil(sizeX/options.size[1]*options.grid[1]))+1
gridx, gridy, gridz = options.grid; sizex, sizey, sizez = options.size
nGrids = gridx*gridy*gridz; avgGrids = nGrids/options.N
dx, dy = options.size[0]/options.grid[0], options.size[1]/options.grid[1]

Ngrains = int(np.ceil(sizeX*sizeY/options.size[0]/options.size[1]*options.N))
filename = 'grains'+str(Ngrains)+'_'+str(gridX)+str(gridY)+str(options.grid[2])
thickness = 1.0-options.reduction

print 'run seeds_fromRandom'
execute('seeds_fromRandom -N %i -g %i %i %i %s.seeds'%(Ngrains,int(gridX*thickness)+1, gridY, int(gridz/thickness)+1, filename))
print 'run geom_fromVoronoiTessellation'
execute('geom_fromVoronoiTessellation -s %s %s %s < %s.seeds'%(sizeX*thickness, sizeY, sizez/thickness,filename))
print 'run geom_rescale'
execute('geom_rescale -g %i %i %i -s %s %s %s < %s.geom > %s_scale.geom'%(gridX, gridY, gridz,
        sizeX, sizeY, sizez,filename,filename) )
print ('the size of the cutted RVE is %sX%sX%s'%(sizex, sizey, sizez) )
print ('the thickness reduction is %s'%('{:.1%}'.format(options.reduction)))

# --- loop over input files -------------------------------------------------------------------------
filenames = [filename+'.geom']

for theta in np.linspace(0, np.pi/2, options.number):
  postfix = str(int(np.round(theta*180.0/np.pi)))
  c_t, s_t = np.cos(theta), np.sin(theta)
  offsetX, offsetY = 0.5*( sizeX - (sizex*c_t-sizey*s_t) ), 0.5*( sizeY - (sizex*s_t+sizey*c_t) )
  for name in filenames:
    if name == 'STDIN':
      file = {'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout, 'croak':sys.stderr}
      file['croak'].write('\033[1m'+scriptName+'\033[0m\n')
    else:
      if not os.path.exists(name): continue
      file = {'name':name, 'input':open(name), 'output':open(name+postfix+'_tmp','w'), 'croak':sys.stderr}
      file['croak'].write('\033[1m'+scriptName+'\033[0m: '+file['name']+'\n')

    table = damask.ASCIItable(file['input'],file['output'],buffered=False)                            # make unbuffered ASCII_table
    table.head_read()                                                                                 # read ASCII header info
#--- interpret header ----------------------------------------------------------------------------
    info = {
          'grid':    np.zeros(3,'i'),
          'size':    np.array((0.0,0.0,0.0)),
          'origin':  np.zeros(3,'d'),
          'microstructures':  0,
          'homogenization': 0,
           }
    newInfo = {
          'microstructures': 0,
           }
    extra_header = []

    for header in table.info:
      headitems = map(str.lower,header.split())
      if len(headitems) == 0: continue
      if headitems[0] in mappings.keys():
        if headitems[0] in identifiers.keys():
          for i in xrange(len(identifiers[headitems[0]])):
            info[headitems[0]][i] = \
              mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])
        else:
          info[headitems[0]] = mappings[headitems[0]](headitems[1])
      else:
        extra_header.append(header)
    newInfo['microstructures'] = info['microstructures']
    if 0 not in options.grid:                                                                         # user-specified grid
      info['grid'] = np.array(options.grid)

    for i in xrange(3):
      if info['size'][i] <= 0.0:                                                                      # any invalid size?
        info['size'][i] = float(info['grid'][i])/max(info['grid'])
        file['croak'].write('rescaling size %s...\n'%{0:'x',1:'y',2:'z'}[i])
    if np.any(info['grid'] < 1):
      file['croak'].write('invalid grid a b c.\n')
      continue
    if np.any(info['size'] <= 0.0):
      file['croak'].write('invalid size x y z.\n')
      continue

    # read the topological data from maternal RVE file
    print ( 'cut the RVE along the direction of theta = %i'%(np.round(theta*180/np.pi)) )
    GrainNo = np.chararray((gridX, gridY, gridz), itemsize=6)
    for i in xrange(gridY*gridz):
      content = file['input'].readline().split()
      for j in xrange(gridX): GrainNo[j, np.mod(i,gridY), i/gridY] = content[j]

    # cut a sub-RVE cooresponding to the specified direction from the maternal RVE
    subGrainNo = np.chararray((gridx, gridz, gridy), itemsize=6)
    for i in xrange(gridx):
      for j in xrange(gridy):
        I = int(np.floor( (dx*(i+0.5)*c_t - dy*(j+0.5)*s_t + offsetX)/dx))
        J = int(np.floor( (dx*(i+0.5)*s_t + dy*(j+0.5)*c_t + offsetY)/dy))
        I = min(I, gridX-1); J = min(J, gridY-1)
        for k in xrange(gridz): subGrainNo[i,k,j] = GrainNo[I,J,k]
    subGrainNoVec = subGrainNo.reshape(nGrids)

    # count the number of grains in the sub-RVE
    index = defaultdict(list)
    for i in xrange(nGrids): index[subGrainNoVec[i]].append(i)
    ngrains = len(index)

    # count the broken (scattered) grains due to the cutting, and merge them.
    if ngrains > options.N*1.1:
      N1, N2 = 0, nGrids
      for key in index:
        if len(index[key])>=0.4*avgGrids: 
          N1+=1; N2-=len(index[key])  # N1: valid grains; N2: number of grids needed to be re-assigned orientation
      ngrid2 = min(int(0.8*avgGrids), N2/(options.N-N1)+1) # grid in each grain

      newGrains = N2/ngrid2+1
      a = [ngrid2]*(newGrains)
      for key in index:
        if len(index[key])<0.4*avgGrids:
          for i in xrange(newGrains):
            if a[i] >= len(index[key]):
              for j in index[key]: subGrainNoVec[j] = -i-1
              a[i] -= len(index[key])
              break

    index = defaultdict(list)
    for i in xrange(nGrids): index[subGrainNoVec[i]].append(i)
    ngrains = len(index)

    # assign orientations
    grainsNo = np.arange(ngrains); np.random.shuffle(grainsNo)
    for i,key in enumerate(index):
      no = str(grainsNo[i]+1)
      for j in index[key]: subGrainNoVec[int(j)] = no

    newInfo['microstructures'] = ngrains
    if newInfo['microstructures'] == 0:
      file['croak'].write('no grain info found.\n')
      continue

#--- write header ---------------------------------------------------------------------------------
    table.labels_clear()
    table.info_clear()
    table.info_append(extra_header+[
      scriptID + ' ' + ' '.join(sys.argv[1:]),
      "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2],),
      "size\tx %f\ty %f\tz %f"%(options.size[0],options.size[1],options.size[2],),
      "origin\tx %f\ty %f\tz %f"%(info['origin'][0],info['origin'][1],info['origin'][2],),
      "homogenization\t%i"%info['homogenization'],
      "microstructures\t%i"%(newInfo['microstructures']),
      ])
    table.head_write()
# --- write microstructure information ------------------------------------------------------------
    formatwidth = 1+int(math.log10(newInfo['microstructures']))
    table.data = np.array(map(int, subGrainNoVec)).reshape(gridx, gridy*gridz).T
    table.data_writeArray('%%%ii'%(formatwidth),delimiter=' ')

#--- output finalization -------------------------------------------------------------------------- 
    if file['name'] != 'STDIN':
      prefix = os.path.splitext(file['name'])[0].replace('grains'+str(Ngrains),'grains'+str(ngrains))
      os.rename(name+postfix+'_tmp', 
         prefix+'%s'%('_material.config' if options.config else '_'+postfix+'.geom'))
  table.close()
