### Convert a multi-line fasta to single-line
awk '{if(match($0,/^>/)){if(length(SEQ)>0){print HEADER"\n"SEQ};HEADER=$0;SEQ=""}else{SEQ=SEQ$0}}' input.fasta > output.fasta

### calculate GC (this is  about twice as fast as using grep and tee, and about 10x faster than Perl)
awk '{if(!gsub(/">"/,"")) {$0=toupper($0);COUNT=COUNT+gsub(/G|C/,"")}} END {print COUNT}' input.fasta
awk '{if(!gsub(/">"/,"")) {$0=toupper($0);COUNT=COUNT+gsub(/A|T/,"")}} END {print COUNT}' input.fasta
awk '{if(!gsub(/">"/,"")) {$0=toupper($0);COUNT=COUNT+gsub(/[^ATCG]/,"")}} END {print COUNT}' input.fasta
# here is the single command version (if you want to count ambigueos bases, write something yourself (actually I might do this as it would be fun))
awk '{if(!gsub(/">"/,"")) {$0=toupper($0);COUNT_GC=COUNT_GC+gsub(/G|C/,"");COUNT_AT=COUNT_AT+gsub(/A|T/,"");COUNT_ELSE=COUNT_ELSE+gsub(/[^ATCG]/,"")}} END {print "GC:\t"COUNT_GC"\nAT:\t"COUNT_AT"\nOthers:\t"COUNT_ELSE}' input.fasta

### extract sequences from a fasta file based on a list of unformated headers (use the above script to convert the fasta file to single line)
# This is a bit of a hack which will search through the fasta file for every header in your list - if either the fasta file or list are long use a different method
# -P(perl) -o(only matching) grep for a partial header in a text file. Second part returns matches to the header +1 line
grep -Po "FGSG_\d*" a_list.txt|xargs -I HEADER grep -A 1 HEADER singled.fasta > a_list.fa

### some USEARCH stuff
# create a udb database from a transcriptome
usearch -makeudb_usearch QUORN_transcriptome.fa -output QUORN_transcriptome.udb
# local search example
usearch -usearch_local a_list.fa -db QUORN_transcriptome.udb -id 0.7 -evalue 1e-25 -strand plus -blast6out results.out 
# global seach example
usearch -usearch_global a_list.fa -db QUORN_transcriptome.udb -id 0.7 -evalue 1e-25 -strand plus -blast6out results.out 
