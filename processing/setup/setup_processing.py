#!/usr/bin/env python

# Makes postprocessing routines acessible from everywhere.

import os,sys,glob,string
from damask import Environment
from optparse import OptionParser, Option


# -----------------------------
class extendableOption(Option):
# -----------------------------
# used for definition of new option parser action 'extend', which enables to take multiple option arguments
# taken from online tutorial http://docs.python.org/library/optparse.html
  
  ACTIONS = Option.ACTIONS + ("extend",)
  STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
  TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
  ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

  def take_action(self, action, dest, opt, value, values, parser):
    if action == "extend":
      lvalue = value.split(",")
      values.ensure_value(dest, []).extend(lvalue)
    else:
      Option.take_action(self, action, dest, opt, value, values, parser)

      
      
########################################################
# MAIN
########################################################

parser = OptionParser(option_class=extendableOption, usage='%prog options', description = """
Sets up the pre and post processing tools of DAMASK

""" + string.replace('$Id$','\n','\\n')
)

compilers = ['intel','ifort','intel32','gfortran','gnu95']

parser.add_option('--F90', '--f90',     dest='compiler', type='string', \
                                        help='name of F90 compiler')
                                        
parser.set_defaults(compiler  = 'ifort')  
(options,filenames) = parser.parse_args()

if options.compiler not in compilers:
  parser.error('compiler switch "--F90" has to be one out of: %s'%(', '.join(compilers)))

f2py_compiler = {
                  'gfortran': 'gnu95   --f90flags="-fno-range-check -xf95-cpp-input -std=f2008 -fall-intrinsics -DSpectral -fdefault-real-8 -fdefault-double-8 -DFLOAT=8 -DINT=4 -I${DAMASK_ROOT}/lib"',
                  'gnu95':    'gnu95   --f90flags="-fno-range-check -xf95-cpp-input -std=f2008 -fall-intrinsics -DSpectral -fdefault-real-8 -fdefault-double-8 -DFLOAT=8 -DINT=4 -I${DAMASK_ROOT}/lib"',
                  'intel32':  'intel   --f90flags="-fpp -stand f03 -diag-disable 5268 -assume byterecl -DSpectral -real-size 64 -integer-size 32 -DFLOAT=8 -DINT=4 -I${DAMASK_ROOT}/lib"',
                  'intel':    'intelem --f90flags="-fpp -stand f03 -diag-disable 5268 -assume byterecl -DSpectral -real-size 64 -integer-size 32 -DFLOAT=8 -DINT=4 -I${DAMASK_ROOT}/lib"',
                  'ifort':    'intelem --f90flags="-fpp -stand f03 -diag-disable 5268 -assume byterecl -DSpectral -real-size 64 -integer-size 32 -DFLOAT=8 -DINT=4 -I${DAMASK_ROOT}/lib"',
                }[options.compiler]
compiler = {
                  'gfortran': 'gfortran',
                  'gnu95':    'gfortran',
                  'intel32':  'ifort',
                  'intel':    'ifort',
                  'ifort':    'ifort',
                }[options.compiler]

damaskEnv = Environment()
baseDir = damaskEnv.relPath('processing/')
codeDir = damaskEnv.relPath('code/')

if   'ikml' in damaskEnv.pathInfo and damaskEnv.pathInfo['ikml'] != '':
  lib_lapack = ''  # TODO!!
elif 'acml' in damaskEnv.pathInfo and damaskEnv.pathInfo['acml'] != '':
  lib_lapack = '-L%s/%s64/lib -lacml'%(os.path.join(damaskEnv.pathInfo['acml']),compiler)                                    # can we use linker flag?
elif 'lapack' in damaskEnv.pathInfo and damaskEnv.pathInfo['lapack'] != '':
  lib_lapack = '-L%s -llapack'%(damaskEnv.pathInfo['lapack'])                                     # see http://cens.ioc.ee/pipermail/f2py-users/2003-December/000621.html

#define ToDo list
bin_link = { \
        'pre' : [
                'marc_addUserOutput.py',
                'mentat_pbcOnBoxMesh.py',
                'mentat_spectralBox.py',
                'OIMang_hex2cub.py',
                'patchFromReconstructedBoundaries.py',
                'spectral_geomCanvas.py',
                'spectral_geomCheck.py',
                'spectral_geomFromAng.py',
                'spectral_geomFromMinimalSurface.py',
                'spectral_geomFromVoronoiTessellation.py',
                'spectral_geomMultiply.py',
                'spectral_geomPack.py',
                'spectral_geomTranslate.py',
                'spectral_geomUnpack.py',
                'spectral_geomVicinityOffset.py',
                'spectral_randomSeeding.py',
                ],
        'post' : [
                '3Dvisualize.py',
                'addCauchy.py',
                'addCalculation.py',
                'addDeterminant.py',
                'addDivergence.py',
                'addCurl.py',
                'addMises.py',
                'addNorm.py',
                'addPK2.py',
                'addStrainTensors.py',
                'addCompatibilityMismatch.py',
                'addDeformedConfiguration.py',
                'averageDown.py',
                'binXY.py',
                'deleteColumn.py',
                'deleteInfo.py',
                'filterTable.py',
                'mentat_colorMap.py',
                'postResults.py',
                'spectral_iterationCount.py',
                'spectral_parseLog.py',
                'nodesFromCentroids.py',
                'tagLabel.py',
                'table2ang',
                ],
            }

compile = { \
        'pre' : [
                ],
        'post' : [
                ],
            }
            

execute = { \
          'coreModule' : [ 
                        'make tidy',
                        # The following command is used to compile the fortran files and make the functions defined
                        # in damask.core.pyf available for python in the module core.so
                        # It uses the fortran wrapper f2py that is included in the numpy package to construct the
                        # module core.so out of the fortran code in the f90 files
                        # For the generation of the pyf file use the following lines:
                        ###########################################################################
                        #'f2py -h damask.core.pyf' +\
                        #' --overwrite-signature --no-lower prec.f90 DAMASK_spectral_interface.f90 math.f90 mesh.f90',
                        ###########################################################################
                        'f2py damask.core.pyf' +\
                        ' --build-dir ./' +\
                        ' -c --no-lower --fcompiler=%s'%(f2py_compiler) +\
                        ' %s'%'prec.f90'+\
                        ' %s'%'DAMASK_spectral_interface.f90'+\
                        ' %s'%'IO.f90'+\
                        ' %s'%'numerics.f90'+\
                        ' %s'%'debug.f90'+\
                        ' %s'%'math.f90'+\
                        ' %s'%'FEsolving.f90'+\
                        ' %s'%'mesh.f90'+\
                        ' %s'%'spectral_quit.f90'+\
                        ' -L%s/lib -lfftw3'%(damaskEnv.pathInfo['fftw'])+\
                        ' %s'%lib_lapack,
                        'mv %s `readlink -f %s`' %(os.path.join(codeDir,'core.so'),os.path.join(damaskEnv.relPath('lib/damask'),'core.so')),
                        ]
            }


for dir in compile:
  for file in compile[dir]:
    src = os.path.abspath(os.path.join(baseDir,dir,file))
    if os.path.isfile(src):
      try:
        os.system('rm %s.exe'%(os.path.splitext(src)[0]))
        print 'removing %s.exe '%(os.path.splitext(src)[0])
      except:
        pass
      print 'compiling ',src,'using',compiler
      os.system('%s -O2 -o %s.exe %s'%(compiler,os.path.splitext(src)[0],src))
 
os.chdir(codeDir)                   # needed for compilation with gfortran and f2py
for tasks in execute:
  for cmd in execute[tasks]:
    try:
      print 'executing...:',cmd
      os.system(cmd)
    except:
      print 'failed..!'
      pass

os.chdir(damaskEnv.relPath('processing/setup/'))
modules = glob.glob('*.mod')
for module in modules:
  print 'removing', module
  os.remove(module)
  
for dir in bin_link:
  for file in bin_link[dir]:
    src = os.path.abspath(os.path.join(baseDir,dir,file))
    if (file == ''):
      sym_link = os.path.abspath(os.path.join(damaskEnv.binDir(),dir))
    else:
      sym_link = os.path.abspath(os.path.join(damaskEnv.binDir(),os.path.splitext(file)[0]))
    print sym_link,'-->',src
    if os.path.lexists(sym_link):
      os.remove(sym_link)    
    os.symlink(src,sym_link)
    
