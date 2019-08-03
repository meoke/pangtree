cd ..
THEDATE=`date +%d%m_%H%M%S`
OUTPUTDIR="data/Simulated1/output/"analysis_$THEDATE
mkdir $OUTPUTDIR
STATSFILEPATH=$OUTPUTDIR"/stats.csv"
touch $STATSFILEPATH
echo "id,sample,p,runtime,memorypick" >> $STATSFILEPATH

source venv/bin/activate
i=0
for sample in {0..2}; do
  values=(0.25 0.5 1 2 4)
  labels=("025" "05" "1" "2" "4")
  for p in {0..4}; do
    p_value=${values[$p]}
    i=$((i+1))
    ITERDIR=$OUTPUTDIR"/"$sample"_"${labels[$p]}"_output"
    mkdir $ITERDIR
    v=$( { /usr/bin/time -f "%e,%M" python3 -m pangtreebuild --output_dir $ITERDIR --multialignment "data/Simulated1/input/multialignment"$sample".maf" --consensus tree --p $p_value --stop 0.99 --metadata "data/Simulated1/input/metadata.csv" --q; } 2>&1 )
    echo $i","$sample","$p_value","$v >> $STATSFILEPATH
    done
done

python3 scripts/make_charts.py $STATSFILEPATH "Simulated1" $OUTPUTDIR"/time_chart.png" $OUTPUTDIR"/memory_chart.png"

deactivate
