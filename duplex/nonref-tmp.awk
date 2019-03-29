BEGIN {
  FS = "\t"; OFS = "\t";
}
$2 != last_bar {
  print NR, "new bar";
  if ($1 == last_code) {
    print NR, "same code";
    diffs1 = diffs;
    reads1 = reads;
    bar1 = last_bar;
  } else {
    print NR, "new code";
    if (reads1) {
      print NR, "reads1:", reads1;
      diffs2 = diffs;
      reads2 = reads;
      bar2 = last_bar;
      avg1 = diffs1/reads1;
      avg2 = diffs2/reads2;
      diff_avg = avg1 - avg2;
      if (diff_avg < 0) diff_avg = -diff_avg;
      print last_code, sprintf("%0.1f", diff_avg), bar1, sprintf("%0.1f", diffs1/reads1), bar2, sprintf("%0.1f", diffs2/reads2);
    }
    diffs1 = 0;
    reads1 = 0;
    diffs2 = 0;
    reads2 = 0;
  }
  diffs = 0;
  reads = 0;
}
{
  last_code = $1;
  last_bar = $2;
  diffs += $4;
  reads++;
}

