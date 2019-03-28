#!/usr/bin/env python3
import argparse
import logging
import os
import pathlib
import subprocess
import sys
import time
assert sys.version_info.major >= 3, 'Python 3 required'

SIM_PARAMS = {
  '--pcr-error': 0.001,
  '--seq-error': 0.001,
  '--n-frags': 100000,
  '--cycles': 25,
  '--read-len': 251,
  '--frag-len': 600,
  '--efficiency-decline': 1.05,
}
PROCESSES = 31
USAGE = "%(prog)s [options] ref out/dir -- [sim_arg1 [sim_arg2 [..]]"

def make_argparser():
  parser = argparse.ArgumentParser(usage=USAGE)
  parser.add_argument('ref', type=pathlib.Path)
  parser.add_argument('outdir', type=pathlib.Path)
  parser.add_argument('sim_args', nargs='+',
    help='Arguments for sim.py.')
  parser.add_argument('--no-execute', '-n', dest='execute', action='store_false', default=True)
  parser.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = parser.add_mutually_exclusive_group()
  volume.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser

def main(argv):
  parser = make_argparser()
  args = parser.parse_args(argv[1:])
  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  this_script_dir = pathlib.Path(__file__).resolve().parent
  commands = []

  if args.execute:
    with (args.outdir/'started.txt').open('w') as started_file:
      started_file.write('{:0.0f}\n'.format(time.time()))

  # Run simulator.
  sim_outdir = args.outdir / 'sim'
  sim_paths = {
    'ref': args.ref,
    'frag': sim_outdir / 'frags.fq',
    'muts': sim_outdir / 'muts.tsv',
    'bars': sim_outdir / 'bars.tsv',
    'fq1': sim_outdir / 'reads_1.fq',
    'fq2': sim_outdir / 'reads_2.fq',
  }
  sim_script = pathlib.Path('~/code/dunovo-clean/utils/sim.py').expanduser()
  param_args = get_sim_param_args(SIM_PARAMS, args.sim_args)
  sim_cmd = [
    str(sim_script), '--out-format', 'fastq', '--seed', '1337'] + param_args + [args.ref,
    '--frag-file', sim_paths['frag'], '--mutations', sim_paths['muts'],
    '--barcodes', sim_paths['bars'], '-1', sim_paths['fq1'], '-2', sim_paths['fq2']
  ]
  commands.append({'cmd':sim_cmd, 'stderr':None})
  if not sim_script.is_file():
    fail('Error: Could not find {}'.format(sim_script))
  if args.execute and not sim_outdir.is_dir():
    os.mkdir(sim_outdir)

  # Run Du Novo on simulated data.
  dunovo_outdir = args.outdir / 'dunovo'
  dunovo_paths = {
    'tmp': dunovo_outdir,
    'logs': dunovo_outdir / 'logs',
    'fq1': sim_paths['fq1'],
    'fq2': sim_paths['fq2'],
    'out': dunovo_outdir,
    'stderr': dunovo_outdir / 'logs/dunovo.err.log'
  }
  dunovo_script = pathlib.Path('~/code/dunovo-clean/dunovo.py').expanduser()
  dunovo_cmd = [
    dunovo_script, '--verbose', '--tempdir', dunovo_paths['tmp'], '--log-dir', dunovo_paths['logs'],
    '--processes', PROCESSES, '--threads', 4, '--dist', 3, '--cons-thres', 0.7, '--min-reads', 3,
    '--qual', 25, dunovo_paths['fq1'], dunovo_paths['fq2'], '-o', dunovo_paths['out']
  ]
  commands.append({'cmd':dunovo_cmd, 'stderr':dunovo_paths['stderr']})
  if args.execute:
    if not dunovo_outdir.is_dir():
      os.mkdir(dunovo_outdir)
    if not dunovo_paths['logs'].is_dir():
      os.mkdir(dunovo_paths['logs'])

  # Run output through errstats.py.
  err_outdir = args.outdir / 'errstats'
  err_paths = {
    'fams': args.outdir / 'dunovo/families.msa.tsv',
    'sscs': args.outdir / 'dunovo/sscs',
    'out':  err_outdir,
  }
  mkerr_script = this_script_dir / '../make-errstats.sh'
  err_cmd = [
    'bash', mkerr_script, err_paths['fams'], err_paths['sscs'], args.ref, err_paths['out']
  ]
  commands.append({'cmd':err_cmd, 'stderr':None})
  if not mkerr_script.is_file():
    fail('Error: Missing {}'.format(mkerr_script))
  if args.execute and not err_outdir.is_dir():
    os.mkdir(err_outdir)

  # Run commands.
  for cmd in commands:
    command = [str(arg) for arg in cmd['cmd']]
    print('$ '+' '.join(command))
    sys.stdout.flush()
    if not args.execute:
      continue
    kwargs = {}
    stderr_file = None
    try:
      if cmd['stderr'] is not None:
         stderr_file = cmd['stderr'].open('w')
         kwargs['stderr'] = stderr_file
      subprocess.check_call(command, **kwargs)
    finally:
      if stderr_file is not None:
        stderr_file.close()

def get_sim_param_args(static_args, user_args):
  args = []
  final_args = static_args.copy()
  last_arg = None
  for user_arg in user_args:
    if last_arg in static_args:
      final_args[last_arg] = user_arg
    last_arg = user_arg
  for option, value in final_args.items():
    args.append(option)
    args.append(str(value))
  return args

def fail(message):
  logging.critical(message)
  sys.exit(1)

sys.exit(main(sys.argv))
