{
  echo -e "ID\tGene\tERR1797969\tERR1797970\tERR1797971\tERR1797972\tERR1797973\tERR1797974"
  paste \
    <(awk '{print $1}' ERR1797969.counts.txt) \
    <(awk '{if (NF==2) print $1; else print $2}' ERR1797969.counts.txt) \
    <(awk '{print $NF}' ERR1797969.counts.txt) \
    <(awk '{print $NF}' ERR1797970.counts.txt) \
    <(awk '{print $NF}' ERR1797971.counts.txt) \
    <(awk '{print $NF}' ERR1797972.counts.txt) \
    <(awk '{print $NF}' ERR1797973.counts.txt) \
    <(awk '{print $NF}' ERR1797974.counts.txt)
} > counts_E745.txt
