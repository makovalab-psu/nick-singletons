$1 != last {
  count++
  bar = $1
  print ">" count
  print bar
  print ">" count ":rev"
  print swap_halves(bar)
}
{
  last = $1
}
function swap_halves(str) {
  half = length(str)/2
  alpha = substr(str, 1, half)
  beta = substr(str, half+1)
  return beta alpha
}
