# Evalue < 1e-5
awk -F'\t' '$1 !~ /^#/ && $3 < 1e-5 {count++} END{print count}' E745_emapper.emapper.annotations
# lowest, highest, average score
awk -F'\t' '
$1 !~ /^#/ && $4 != "-" && $4 != "" {
  if(n==0 || $4 < min) min=$4;
  if(n==0 || $4 > max) max=$4;
  sum += $4;
  n++;
}
END{
  print "min score:", min;
  print "max score:", max;
  print "mean score:", sum/n;
}' E745_emapper.emapper.annotations
# SCORE
awk -F'\t' '
$1 !~ /^#/ && $1!="query" && $4!="score" && $4!="-" && $4!="" {
  total++;
  if($4 < 50) s50++;
  if($4 < 100) s100++;
  if($4 < 200) s200++;
  if($4 < 500) s500++;
}
END{
  print "Total annotations with score:", total;
  print "score < 50:", s50;
  print "score < 100:", s100;
  print "score < 200:", s200;
  print "score < 500:", s500;
}' E745_emapper.emapper.annotations



