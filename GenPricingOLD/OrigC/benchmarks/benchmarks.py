#!/usr/bin/python
#
# CA: this program executes a number of timing measures on a number of contracts,
#       recording GPU and CPU results and runtimes
# USAGE: run the program, which will produce a log file with the measurements
#

import subprocess as S
import numpy as N
import re
import sys
#-------------------------------------------------------------------------
# global VARS
#-------------------------------------------------------------------------
_TIMING_RECORD_ON=True
#_TEST_RUNS=20 # will be overwritten from command line
_TIMING_RUNS_MULTIPLIER=[1,1,1] #[2,1,1] # on GTX680, the first contract can be so fast that it requires more samples to reduce the timing variance
_DEBUG=5

contracts={
  1:"./examples",
  2:"./barrier_rev_convert_contract",
  3:"./worst_off_contract"
  }

_OPTIMIZATIONS_FLAGS=[
  # first graph
  {'floats': True,  'cost_vect_priv':2, 'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  {'floats': True,  'cost_vect_priv':2, 'sobol_strength_red_recurr': False, 'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  {'floats': True,  'cost_vect_priv':1, 'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  {'floats': True,  'cost_vect_priv':1, 'sobol_strength_red_recurr': False, 'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  #
  {'floats': False, 'cost_vect_priv':2, 'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  {'floats': False, 'cost_vect_priv':2, 'sobol_strength_red_recurr': False, 'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  {'floats': False, 'cost_vect_priv':1, 'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  {'floats': False, 'cost_vect_priv':1, 'sobol_strength_red_recurr': False, 'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  # second graph
  #{'floats': True,  'cost_vect_priv':2,  'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  {'floats': True,  'cost_vect_priv':2, 'sobol_strength_red_recurr': True,  'branch_divergence': False, 'memory_coalescence': True,  'tile_size': 16},
  #{'floats': True,  'cost_vect_priv':1, 'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  {'floats': True,  'cost_vect_priv':1, 'sobol_strength_red_recurr': True,  'branch_divergence': False, 'memory_coalescence': True,  'tile_size': 16},
  #
  #{'floats': False, 'cost_vect_priv':2,  'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  {'floats': False, 'cost_vect_priv':2, 'sobol_strength_red_recurr': True,  'branch_divergence': False, 'memory_coalescence': True,  'tile_size': 16},
  #{'floats': False, 'cost_vect_priv':1, 'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  {'floats': False, 'cost_vect_priv':1, 'sobol_strength_red_recurr': True,  'branch_divergence': False, 'memory_coalescence': True,  'tile_size': 16},
  # third graph
  #{'floats': True,  'cost_vect_priv':2,  'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  {'floats': True,  'cost_vect_priv':2, 'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': False, 'tile_size': 16},
  #{'floats': True,  'cost_vect_priv':1, 'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  {'floats': True,  'cost_vect_priv':1, 'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': False, 'tile_size': 16},
  #
  #{'floats': False, 'cost_vect_priv':2, 'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  {'floats': False, 'cost_vect_priv':2, 'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': False, 'tile_size': 16},
  #{'floats': False, 'cost_vect_priv':1, 'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
  {'floats': False, 'cost_vect_priv':1, 'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': False, 'tile_size': 16},
  # tiling
  #{'cost_vect_priv':True,  'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size':   8}, # segmentation fault in contract 1
  #{'cost_vect_priv':True,  'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size':  32},
  #{'cost_vect_priv':True,  'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size':  64},
  #{'cost_vect_priv':True,  'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 128}, # segmentation fault in contract 2
  #
  # leaving with full optimization, with cost model depending on hardware and contract type
  {'floats': True,  'cost_vect_priv':0, 'sobol_strength_red_recurr': True,  'branch_divergence': True,  'memory_coalescence': True,  'tile_size': 16},
]

_RE_FLAGS=re.IGNORECASE|re.DOTALL|re.MULTILINE
_OPTIMIZATIONS_FILENAME="./Optimizations.h"
_NULL_HANDLE=open("/dev/null","w")
# log file
filename_log=sys.argv[0].replace(".py",".log")
logH=open(filename_log,"w")
#-------------------------------------------------------------------------
def log(msg="",copy_on_stdout=True,eol="\n",fH=logH):
  fH.write("%s%s" % (msg,eol) )
  fH.flush()
  if copy_on_stdout:
    print msg
    sys.stdout.flush()
#-------------------------------------------------------------------------
def run_cmd(cmd):
  p=S.Popen(cmd,stdout=S.PIPE,stderr=S.PIPE)
  (out,err)=p.communicate()
  retcode=p.returncode
  if retcode<0:
    log(err)
    sys.exit(retcode)
  else:
    return (out,err)
#-------------------------------------------------------------------------
def optimizations_write(
    optimization_flags,
    filename=_OPTIMIZATIONS_FILENAME
    ):
  fH=open(filename,"w")
  # floats/double
  fH.write("#define _OPTIMIZATION_USE_FLOATS ")
  if optimization_flags['floats']: fH.write("1\n")
  else: fH.write("0\n")
  # cost model: 0, vectorized: 1, privatized: 2
  fH.write("#define _OPTIMIZATION_COST_MODEL_OR_FORCE_VECT_OR_FORCE_PRIV %d " % optimization_flags['cost_vect_priv'])
  fH.write("// cost model: 0, vectorized: 1, privatized: 2\n");
  # Sobol: strength reduction
  if not optimization_flags['sobol_strength_red_recurr']: fH.write("// ")
  fH.write("#define _OPTIMIZATION_SOBOL_STRENGTH_RED_RECURR\n")
  # branch divergence
  if not optimization_flags['branch_divergence']: fH.write("// ")
  fH.write("#define _OPTIMIZATION_BRANCH_OPT\n")
  # memory coalescence
  if not optimization_flags['memory_coalescence']: fH.write("// ")
  fH.write("#define _OPTIMIZATION_MEM_COALES_ON\n")
  # tiling size
  if not optimization_flags['tile_size']: fH.write("// ")
  fH.write("#define _OPTIMIZATION_TILE %s\n" % optimization_flags['tile_size'])
  #
  fH.close()
#-------------------------------------------------------------------------
def make_vars_prepare(contract_num,pricer_dir):
  return [
    'PRICER_DIR="%s"' % pricer_dir,
    'CONTR_NUM=%d' % contract_num
    ]
#-------------------------------------------------------------------------
def make_clean_compile(contract_num,pricer_dir):
  if _DEBUG>=20:
    log( "# DEBUG: make clean && make compile: contract %d" % contract_num )
  #
  cmd=['make','clean']; cmd.extend( make_vars_prepare(contract_num,pricer_dir) )
  if _DEBUG>=50: print "# DEBUG: cmd: '%s'" % ' '.join(cmd)
  run_cmd(cmd)
  #
  cmd=['make','compile']; cmd.extend( make_vars_prepare(contract_num,pricer_dir) )
  if _DEBUG>=50: print "# DEBUG: cmd: '%s'" % ' '.join(cmd)
  run_cmd(cmd)
#-------------------------------------------------------------------------
def result_gpu_cpu_sanity_check(run_output,sanity_show_flag=False):
  result_ok_flag=False
  try:
    (result_gpu,result_cpu)=map(
      float,
      re.search("FINAL RESULT GPU.+?IS:\s+([\d\.]+).+?FINAL RESULT CPU.+?IS:\s+([\d\.]+)",run_output,_RE_FLAGS).groups()
      )
    if abs(result_gpu/result_cpu-1)<1e-5:
      result_ok_flag=True
      if sanity_show_flag==True or _DEBUG>=10:
        log( "# DEBUG: result of contract %d: OK: consistent between GPU and CPU (%.3f,%.3f)" % (contract_num,result_gpu,result_cpu) )
    else:
      log( "# WARNING: result of contract %d: NOT consistent between GPU and CPU (%.3f,%.3f)" % (contract_num,result_gpu,result_cpu) )
  except:
    log( "# WARNING: run of contract %d: NOT readable" % (contract_num,) )
  #
  return result_ok_flag
#-------------------------------------------------------------------------
def test_run_contract(contract_num,pricer_dir,run_id="",sanity_show_flag=False):
  if _DEBUG>=20:
    log( "# DEBUG: make run_pricer: contract %d" % contract_num )
  cmd=['make','run_pricer']; cmd.extend( make_vars_prepare(contract_num,pricer_dir) )
  if _DEBUG>=50: print "# DEBUG: cmd: '%s'" % ' '.join(cmd)
  (run_output,err)=run_cmd(cmd)
  # sanity check for output coherence
  result_ok_flag=result_gpu_cpu_sanity_check(run_output,sanity_show_flag=sanity_show_flag)
  (timing_gpu,timing_cpu)=(0,0)
  if result_ok_flag:
    # record timing in microseconds
    (timing_gpu,timing_cpu)=map(
      float,
      re.search("GPU Run Time:\s+([\d\.]+).+?CPU Run Time:\s+([\d\.]+)",run_output,_RE_FLAGS).groups()
      )
    if _DEBUG>=5:
      log( "# DEBUG: timing %s of contract %d: %d,%d (%.3fx)" % (run_id,contract_num,timing_gpu,timing_cpu,timing_cpu/timing_gpu) )
  #
  return (timing_gpu,timing_cpu,result_ok_flag)
#-------------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------------
if len(sys.argv)<2:
  print "# USAGE: %s [number_of_tests_per_configuration]" % sys.argv[0]
  sys.exit(1)
# read how many tests to run
_TEST_RUNS=int(sys.argv[1])

# log file
print "# DEBUG: saving output to:",filename_log
log("# DEBUG: timing executed %d times" % _TEST_RUNS)

# set optimization configuration
tests_db={}
for (opt_j,optimization_set) in enumerate(_OPTIMIZATIONS_FLAGS):
  optimizations_write(optimization_set,_OPTIMIZATIONS_FILENAME)
  if _DEBUG>=5:
    log("# DEBUG: optimization set %d: %s" % (opt_j, optimization_set))
  # run tests over all the contracts
  tests_db[opt_j]={}
  for contract_j,contract_num in enumerate(sorted(contracts.keys())):
    pricer_dir=contracts[contract_num]
    # clean and compile
    make_clean_compile(contract_num=contract_num,pricer_dir=pricer_dir)
    # run timing tests per contract
    if _TIMING_RECORD_ON:
      tests_db[opt_j][contract_num]=[]
      sanity_show_flag=True
      # how many runs?
      runs_number=_TEST_RUNS*_TIMING_RUNS_MULTIPLIER[contract_j]
      # perform the timing
      for runs_j in xrange(runs_number):
        (timing_gpu,timing_cpu,result_ok_flag)=test_run_contract(
          run_id="(%2d/%2d)" % (runs_j+1,runs_number),
          contract_num=contract_num,pricer_dir=pricer_dir,
          sanity_show_flag=sanity_show_flag
          )
        sanity_show_flag=False # don't show results after the first run
        if result_ok_flag:
          tests_db[opt_j][contract_num].append( (timing_gpu,timing_cpu) )
        else:
          break # error in one of these runs? probably errors in all the others as well. Exiting
      # statistics : point estimate of the Cauchy distribution given as ratio of two Gaussians
      #    http://mathworld.wolfram.com/NormalRatioDistribution.html
      a=N.array(tests_db[opt_j][contract_num])
      if len(a):
        a=a[:,1]/a[:,0]
        log( "# Speedup for optimization set %d, contract %d: %.3f +/- %.3f" % (
            opt_j,contract_num,a.mean(),a.std()
            )
          )

# output: save timings
log(tests_db)

# log: close
logH.close()
