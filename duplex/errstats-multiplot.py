#!/usr/bin/env python3
import argparse
import collections
import itertools
import math
import matplotlib.pyplot
import logging
import os
import pathlib
import sys
assert sys.version_info.major >= 3, 'Python 3 required'

DISPLAY_NAMES = {
  'exo7ml':        'chrM exonuclease',
  'chrM':          'chrM LR-PCR',
  'pRBIR1a':       'fragile site A',
  'pRBAT1':        'fragile site B',
  'GT7B':          'GT7 microsat',
  'Br-2-G131':     'chrM brain',
  'Oo-b11sh-lg':   'chrM linear oocyte',
  'G131-Br-sh-lg': 'chrM linear Br',
  'G131-M-sh-lg':  'chrM linear M',
  'G131-Ht-sh-lg': 'chrM linear Ht',
  'simulated':     'simulated',
}


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.add_argument('errstats', nargs='*', type=expanded_path,
    help='The path to a supercollated errstats tsv. The sample name will be inferred from the path. '
         'This works for normal sample directories under duplex/multi or duplex/sim.')
  parser.add_argument('-s', '--sample', dest='sample_strs', action='append',
    help='Samples and their supercollated errstats tsv. Each is a comma-delimited string. '
         'If 2 values: the sample name and the path to the errstats tsv. '
         'If 3 values: the sample name, its display name, and the path to the errstats tsv.')
  parser.add_argument('-S', '--sim-errstats', type=expanded_path,
    help='The supercollated output of errstats.py for the simulated data (if any).')
  parser.add_argument('-o', '--outdir', type=expanded_path, required=True,
    help='Where to save the plots.')
  parser.add_argument('-n', '--display-name', dest='display_names', action='append',
    help='The display name for a sample. Give two comma-delimited values: The sample id, and the '
         'display name.')
  parser.add_argument('-m', '--max-famsize', type=int,
    help='The maximum family size to plot. Default: The largest family size that appears in all '
         'the samples.')
  parser.add_argument('-X', '--x-label', default='Proportion of total errors',
    help='X axis label. Default: %(default)s')
  parser.add_argument('-Y', '--y-label', default='# of reads error occurs in',
    help='Y axis label. Default: %(default)s')
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

  errstats_paths, display_names = parse_sample_strings(args.sample_strs)
  for errstats_path in args.errstats:
    sample_id = get_sample_id_from_path(errstats_path)
    errstats_paths[sample_id] = errstats_path
  if args.sim_errstats:
    sample_info['simulated'] = {'errstats':args.sim_errstats, 'name':'simulated'}

  for display_name_str in args.display_names:
    sample_id, display_name = parse_display_name(display_name_str)
    display_names[sample_id] = display_name
  for sample_id, display_name in DISPLAY_NAMES.items():
    if sample_id not in display_names:
      display_names[sample_id] = display_name

  all_errstats = {}
  for sample, errstats_path in sample_info.items():
    all_errstats[sample] = read_errstats(errstats_path)

  # Get the min/max family sizes in the samples.
  min_famsize, max_famsize = get_famsize_extents(all_errstats)
  if args.max_famsize:
    max_famsize = args.max_famsize

  # Make individual plots.
  for famsize in range(min_famsize, max_famsize+1):
    figure = matplotlib.pyplot.figure(dpi=180, figsize=(8, 6))
    axes = figure.add_subplot(1, 1, 1)
    handles = plot_famsize(axes, famsize, all_errstats)
    set_ticks(axes, famsize, -5.5, 0.5)
    axes.legend(handles=handles, labels=[sample_names[sample] for sample in all_errstats])
    axes.set_xlabel(args.x_label)
    axes.set_ylabel(args.y_label)
    matplotlib.pyplot.savefig(args.outdir / 'famsize{}.png'.format(famsize))

  # Make multiplot.
  figure = matplotlib.pyplot.figure(dpi=180, figsize=(12, 9))
  figure.subplots(3, 3, sharex=True, sharey=True, squeeze=True)
  for i, axes in enumerate(figure.axes):
    famsize = (i+1)*2+3
    handles = plot_famsize(axes, famsize, all_errstats)
  set_ticks(axes, famsize, -5.5, 0.5)
  figure.legend(handles=handles, labels=[sample_names[sample] for sample in all_errstats])
  figure.suptitle('Repeated errors', size=20)
  figure.text(0.055, 0.64, args.x_label, fontsize=16, horizontalalignment='center', rotation='vertical')
  figure.text(0.51, 0.06, args.y_label, fontsize=16, horizontalalignment='center')
  matplotlib.pyplot.savefig(args.outdir / 'combined.png')


def read_errstats(errstats_path):
  errstats = []
  with errstats_path.open() as errstats_file:
    for line in errstats_file:
      fields = line.rstrip('\r\n').split('\t')
      if len(fields) != 3:
        continue
      # Fields: family size, repeat count, count
      line = [int(value) for value in fields]
      errstats.append(line)
  return errstats


def get_famsize_extents(all_errstats):
  max_min_famsize = 0
  min_max_famsize = None
  logging.info('sample       min   max')
  for sample, errstats in all_errstats.items():
    min_famsize = min([line[0] for line in errstats])
    max_famsize = max([line[0] for line in errstats])
    logging.info('{:13s}  {}  {:4d}'.format(sample, min_famsize, max_famsize))
    if min_famsize > max_min_famsize:
      max_min_famsize = min_famsize
    if min_max_famsize is None or max_famsize < min_max_famsize:
      min_max_famsize = max_famsize
  return max_min_famsize, min_max_famsize


def plot_famsize(axes, famsize, all_errstats, ymin=-5.5):
  ymax = 0.5
  handles = []
  for sample, errstats in all_errstats.items():
    total_errors = sum([line[2] for line in errstats if line[0] == famsize])
    repeats = [line[1] for line in errstats if line[0] == famsize]
    counts = [math.log10(line[2]/total_errors) for line in errstats if line[0] == famsize]
    retval = axes.plot(repeats, counts)
    line = retval[0]
    handles.append(line)
    if sample == 'simulated':
      line.set_color('black')
  axes.set_title('Families with {} reads'.format(famsize))
  return handles


def set_ticks(axes, largest_famsize, ymin=-5.5, ymax=0.5):
  xmax = 8.5 + 5*int((largest_famsize-8)/10)
  axes.set_xlim(0.5, xmax)
  axes.set_ylim(ymin, ymax)
  xstep = int(xmax/9)+1
  xticks = list(range(xstep, int(xmax+1), xstep))
  axes.set_xticks(xticks)
  yticks = list(range(math.ceil(ymin), math.ceil(ymax)))
  axes.set_yticks(yticks)
  yticklabels = []
  for tick in yticks:
    pct = 100*10**tick
    if pct >= 1.0:
      pct = int(pct)
    yticklabels.append(str(pct)+'%')
  axes.set_yticklabels(yticklabels)


def expanded_path(path_str):
  return pathlib.Path(path_str).expanduser()


def parse_sample_strings(sample_strings):
  errstats_paths = collections.OrderedDict()
  display_names = collections.OrderedDict()
  if sample_strings is None:
    return errstats_paths, display_names
  for sample_string in sample_strings:
    fields = sample_string.split(',')
    if len(fields) == 2:
      sample_id, path_str = fields
      disp_name = sample_id
    elif len(fields) == 3:
      sample_id, disp_name, path_str = fields
    else:
      fail('Error: Encountered sample spec with too many or few comma-delimited fields: {!r}'
           .format(sample_string))
    errstats_paths[sample_id] = expanded_path(path_str)
    display_names[sample_id] = disp_name
  return errstats_paths, display_names


def get_sample_id_from_path(path):
  components = path.split(os.sep)
  for i in range(len(components)-1, -1, -1):
    if i-1 >= 0 and components[i-1] == 'multi':
      return components[i]
    elif i-2 >= 0 and components[i-2] == 'sim':
      return components[i]
  return None


def parse_display_name(display_name_str):
  fields = display_name_str.split(',')
  assert len(fields) == 2, display_name_str
  return fields


def fail(message):
  logging.critical(message)
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception('Unrecoverable error')


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except BrokenPipeError:
    pass
