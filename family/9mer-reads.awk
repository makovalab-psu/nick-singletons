BEGIN {
  FS = "\t"
  margin = 5
  mer2 = "CCCCCTCTACCCCCTCTA"
  mer2_len = length(mer2)
  mer3 = "CCCCCTCTACCCCCTCTACCCCCTCTA"
  mer3_len = length(mer3)
}
{
  i1 = index($3, mer2)
  i2 = index($4, mer2)
  from_start1 = i1 - 1
  from_end1 = length($3) - i1 - mer2_len + 1
  from_start2 = i2 - 1
  from_end2 = length($4) - i2 - mer2_len + 1
  if ((from_start1 > margin && from_end1 > margin) ||
      (from_start2 > margin && from_end2 > margin)) {
    i1 = index($3, mer3)
    i2 = index($4, mer3)
    from_start1 = i1 - 1
    from_end1 = length($3) - i1 - mer3_len + 1
    from_start2 = i2 - 1
    from_end2 = length($4) - i2 - mer3_len + 1
    if ((from_start1 > margin && from_end1 > margin) ||
        (from_start2 > margin && from_end2 > margin)) {
      printf("%s\n%s\n+\n%s\n", $1, $3, $7) >> "SC8C1-bl.3mer_1.fq"
      printf("%s\n%s\n+\n%s\n", $2, $4, $8) >> "SC8C1-bl.3mer_2.fq"
    } else {
      printf("%s\n%s\n+\n%s\n", $1, $3, $7) >> "SC8C1-bl.2mer_1.fq"
      printf("%s\n%s\n+\n%s\n", $2, $4, $8) >> "SC8C1-bl.2mer_2.fq"
    }
  }
}

