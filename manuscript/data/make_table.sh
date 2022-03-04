
printf "tool\tsrc\tdest\tcons_src\tcons_dest\tlost_src\tlost_dest\n" > table.tsv
for src in "NCBI" "ensembl" "UCSC" "gencode"
do
    for dest in "NCBI" "ensembl" "UCSC" "gencode"
    do
        if [ "$src" = "$dest" ]; then
            continue
        else
            R_CONS_S=$(wc -l recontig_conserved_source/${src}_2_${dest}.conserved.txt | cut -f1 -d' ')
            R_CONS_D=$(wc -l recontig_conserved_dest/${src}_2_${dest}.conserved.txt | cut -f1 -d' ')
            R_LOST_S=$(wc -l recontig_lost_from_source/${src}_2_${dest}.lost.txt | cut -f1 -d' ')
            R_LOST_D=$(wc -l recontig_lost_from_dest/${src}_2_${dest}.lost.txt | cut -f1 -d' ')

            G_CONS_S=$(wc -l gnu_conserved_source/${src}_2_${dest}.conserved.txt | cut -f1 -d' ')
            G_CONS_D=$(wc -l gnu_conserved_dest/${src}_2_${dest}.conserved.txt | cut -f1 -d' ')
            G_LOST_S=$(wc -l gnu_lost_from_source/${src}_2_${dest}.lost.txt | cut -f1 -d' ')
            G_LOST_D=$(wc -l gnu_lost_from_dest/${src}_2_${dest}.lost.txt | cut -f1 -d' ')

            printf "recontig\t$src\t$dest\t$R_CONS_S\t$R_CONS_D\t$R_LOST_S\t$R_LOST_D\n" >> table.tsv
            printf "gnu\t$src\t$dest\t$G_CONS_S\t$G_CONS_D\t$G_LOST_S\t$G_LOST_D\n" >> table.tsv
        fi
    done
done
