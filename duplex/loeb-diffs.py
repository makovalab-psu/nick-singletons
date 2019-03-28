#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import sys
path1 = sys.argv[1]
path2 = sys.argv[2]
header = ['sample', 'chr', 'pos', 'a', 'c', 'g', 't', 'cvrg', 'alleles', 'major', 'minor', 'maf', 'bias']
def read_vars(infile_path):
  sites = {}
  with open(infile_path) as infile:
    for line in infile:
      fields = line.strip().split()
      site = {}
      for label, value in zip(header, fields):
        site[label] = value
      try:
        pos = int(site['pos'])
      except ValueError:
        continue
      sites[pos] = site
  return sites
sites1 = read_vars(path1)
sites2 = read_vars(path2)
positions = set(sites1.keys())
positions.update(sites2.keys())
print('pos', 'bias1', 'bias2', 'maf1', 'maf2', 'diff', 'pct', sep='\t')
for pos in sorted(positions):
  if pos not in sites1 or pos not in sites2:
    continue
  site1 = sites1[pos]
  site2 = sites2[pos]
  maf1 = float(sites1[pos]['maf'])
  bias1 = sites1[pos]['bias']
  maf2 = float(sites2[pos]['maf'])
  bias2 = sites2[pos]['bias']
  diff = maf1 - maf2
  if maf2 > 0:
    ratio = maf1/maf2
  else:
    ratio = 0
  pct = 100*abs(1 - ratio)
  print(pos, bias1, bias2, maf1, maf2, diff, pct, sep='\t')
