BEGIN {
  OFS = "\t"
  trans["A"] = "T"; trans["C"] = "G"; trans["G"] = "C"; trans["T"] = "A"
}
function revcomp(seq) {
  rc = ""
  for (i = length(seq); i > 0; i--) {
    rc = rc trans[substr(seq, i, 1)]
  }
  return rc
}
function min(a, b) {
  if (a <= b) {
    return a
  } else {
    return b
  }
}
NR % 4 == 2 {
  tot++
  seq = $1
  rc = revcomp($2)
  seq_len = length(seq)
  rc_len = length(rc)
  for (off = -3; off <= 3; off++) {
    if (off < 0) {
      min_len = min(rc_len+off, seq_len)
      if (substr(seq, 1, min_len) == substr(rc, 1-off, min_len)) {
        counts[off]++
      }
    } else {
      min_len = min(rc_len, seq_len-off)
      if (substr(seq, 1+off, min_len) == substr(rc, 1, min_len)) {
        counts[off]++
      }
    }
  }
}
END {
  printf("%s\t", tot)
  for (off in counts) {
    printf("%0.2f\t", 100*counts[off]/tot)
  }
  printf("\n")
}
