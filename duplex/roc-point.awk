BEGIN {
  if (! thres) {
    thres = 0
  }
  if (! min_covg) {
    min_covg = 0
  }
  true = 0
  false = 0
  FS = "\t";
  OFS = "\t";
}

NF == 7 {
  muts[$5] = $7
}

NF == 13 && $12 >= thres && $8 >= min_covg {
  if (muts[$3] && muts[$3] == $11) {
    true++
  } else {
    false++
  }
}

END {
  if (nfalse || ntrue) {
    print false/nfalse, true/ntrue
  } else {
    print false, true
  }
}
