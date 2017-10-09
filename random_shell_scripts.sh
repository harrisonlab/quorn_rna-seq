### Convert a multi-line fasta to single-line
awk '{if(match($0,/^>/)){if(length(SEQ)>0){print HEADER"\n"SEQ};HEADER=$0;SEQ=""}else{SEQ=SEQ$0}}' input.fasta > output.fasta

### calculate GC (this is  about twice as fast as using grep and tee, and about 10x faster than Perl)
awk '{if(!gsub(/">"/,"")) {$0=toupper($0);COUNT=COUNT+gsub(/G|C/,"")}} END {print COUNT}' input.fasta
awk '{if(!gsub(/">"/,"")) {$0=toupper($0);COUNT=COUNT+gsub(/A|T/,"")}} END {print COUNT}' input.fasta
awk '{if(!gsub(/">"/,"")) {$0=toupper($0);COUNT=COUNT+gsub(/[^ATCG]/,"")}} END {print COUNT}' input.fasta

# here is the single command version (if you want to count ambigueos bases, write something yourself (actually I might do this as it would be fun))
awk '{if(!gsub(/">"/,"")) {$0=toupper($0);COUNT_GC=COUNT_GC+gsub(/G|C/,"");COUNT_AT=COUNT_AT+gsub(/A|T/,"");COUNT_ELSE=COUNT_ELSE+gsub(/[^ATCG]/,"")}} END {print "GC:\t"COUNT_GC"\nAT:\t"COUNT_AT"\nOthers:\t"COUNT_ELSE}' input.fasta
