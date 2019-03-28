import sys
low = 0
high = 0
total = 0
line_num = 0
for line in sys.stdin:
  line_num += 1
  if line_num % 4 == 0:
    for char in line:
      score = ord(char)-33
      if char == '\n':
        pass
      elif score > 40:
        high += 1
        char = 'I'
      elif score < 0:
        low += 1
        char = '!'
      sys.stdout.write(char)
    total += len(line)-1
  else:
    sys.stdout.write(line)
sys.stderr.write('low: {}\nhigh: {}\ntotal: {}\n'.format(low, high, total))
