BEGIN { OFS = "" }
{ print "$1 == \"",$1,"\" && $2 >= ",$2," && $2 <= ",$3," { print $0 }" }