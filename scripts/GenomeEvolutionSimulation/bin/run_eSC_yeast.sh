
source treeyeast.sh

r1='root100'

df='fewcycles'
dn='nocycles'


function runsim {
    echo $1 $2 $tname$3
    simCtrl_runSim.py --inputNewick $tree --outDir ../simresults_$1$tname$3'_'$2/ --rootDir ../$1/ --rootName root --paramsDir ../params_$2/ --jobTree ../jobTree_$1$tname$3'_'$2 --stepLength $step --noMEs --maxThreads 4
    cd ../simresults_$1$tname$3'_'$2/
    simCtrl_postSimMafExtractor.py --simDir . --maxBlkWidth=99999
    cd root/
    mafFilter --maf root.maf --includeSeq $tchr > chr20ao.maf
    cat root.maf > chr20ao.maf
    for chr in "${atchr[@]}" 
        do
        mafFilter --maf chr20ao.maf --excludeSeq $chr > chr20ao.maf.new
        mv -T -f chr20ao.maf.new chr20ao.maf
        done
    mafTransitiveClosure --maf chr20ao.maf > chr20ao_tc.maf
    mafDuplicateFilter --maf chr20ao_tc.maf > chr20ao_tc_df.maf
    mafStrander --maf chr20ao_tc_df.maf --seq hg18.chr20 > chr20ao_tc_df_p.maf
    cat chr20ao_tc_df_p.maf | awk '{print $1" "$2" "$3" "$4" "$5}' > chr20ao_tc_df_p.maf.cut
    mafFilter --maf chr20ao_tc_df_p.maf --includeSeq $lchr > chr20ao_tc_df_p_leafs.maf
    simCtrl_postSimFastaExtractor.py --simDir ../ --allCycles
    cd ../../bin
    ./selectfasta.py all $tchr ../simresults_$1$tname$3'_'$2/
    ./selectfasta.py leafs $lchr ../simresults_$1$tname$3'_'$2/
    echo $1 $2 $tname$3

    }

for r in $r1
    do
    for d in $df $dn
        do 
        for n in 0 1 2 3 4 5 6 7 8 9
            do 
            time runsim $r $d $n &
            done
        done
    done


