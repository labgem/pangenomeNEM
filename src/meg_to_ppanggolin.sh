#!/usr/bin/bash

#$1 is pangenome file
#$2 is gff directory
#$3 is modified pangenome file
#$4 is modified gff directory

less $1 | awk 'NF>2{for(i=2;i<=NF;i++){print $1, $i}} NF==2{print $0}' | sed "s/\.[ib]/\./g" > $3
mkdir -p $4
for i in `ls $2`; do less $2/$i | sed "s/^\([^ \t]*\)\(.*\)ID=[^_]*_\(.*;\).*/\1\2ID=\1_\3;/g" > $4/$i; done
#for i in `ls $4`; do less $4/$i | grep -o "ID=[^;]*" | awk -F '=' '{print $2}' | sort > /tmp/pangenome_all; done
#less $3 | awk '{print $2}' | sort > /tmp/pangenome_families
#pangenome_families=$(wc -l $1 | cut -d" " -f1)
#comm -23 /tmp/pangenome_all /tmp/pangenome_families | awk 'BEGIN{id='"$pangenome_families"'} {print id, $1; id++;}' >> $3
