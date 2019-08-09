cd ..
THEDATE=`date +%d%m_%H%M%S`
#dataset_name="Ebola"
dataset_name="Simulated3_nocycles"
#dataset_name="Simulated1"
OUTPUTDIR="data/"$dataset_name"/output/"analysis_$THEDATE
mkdir $OUTPUTDIR
STATSFILEPATH=$OUTPUTDIR"/stats.csv"
touch $STATSFILEPATH
echo "id,sample,p,runtime,memorypick" >> $STATSFILEPATH

source venv/bin/activate
i=0
for sample in {0..1}; do
  values=(0.25 0.5 1 2 4)
  labels=("025" "05" "1" "2" "4")
  for p in {0..4}; do
      p_value=${values[$p]}
      i=$((i+1))
      ITERDIR=$OUTPUTDIR"/"$sample"_"${labels[$p]}"_"$t"_output"
      mkdir $ITERDIR
      v=$( { /usr/bin/time -f "%e,%M" python3 -m pangtreebuild --output_dir $ITERDIR --multialignment "data/"$dataset_name"/input/yeast"$sample"_nocycles.maf" --consensus tree --p $p_value --stop 0.9999 --q; } 2>&1 )
      echo $i","$sample","$p_value","$v >> $STATSFILEPATH
    done
done

python3 scripts/make_charts.py $STATSFILEPATH "$dataset_name" $OUTPUTDIR"/chart.png" $OUTPUTDIR"/data.csv"

deactivate
