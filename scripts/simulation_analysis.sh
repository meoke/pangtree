cd ..
THEDATE=`date +%d%m_%H%M%S`
dataset_name="Ebola"
OUTPUTDIR="data/"$dataset_name"/output/"analysis_$THEDATE
mkdir $OUTPUTDIR
STATSFILEPATH=$OUTPUTDIR"/stats.csv"
touch $STATSFILEPATH
echo "id,sample,p,runtime,memorypick" >> $STATSFILEPATH

source venv/bin/activate
i=0
for sample in {0..9}; do
  values=(0.98 0.99 0.995)
  labels=("098" "099" "0995")
  for hbmin in {0..2}; do
      hbmin_value=${values[$hbmin]}
      i=$((i+1))
      ITERDIR=$OUTPUTDIR"/"$sample"_"${labels[$hbmin]}"_"$t"_output"
      mkdir $ITERDIR
      v=$( { /usr/bin/time -f "%e,%M" python3 -m pangtreebuild --output_dir $ITERDIR --multialignment "data/"$dataset_name"/input/multialignment.maf" --metadata "data/"$dataset_name"/input/metadata.csv"" --fasta_provider ncbi --cache --consensus poa --hbmin $hbmin_value --q; } 2>&1 )
#      v=$( { /usr/bin/time -f "%e,%M" python3 -m pangtreebuild --output_dir $ITERDIR --multialignment "data/"$dataset_name"/input/yeast"$sample"_fewcycles.maf" --consensus tree --p $p_value --stop 0.9999 --q --v ; } 2>&1 )
      echo $i","$sample","$hbmin_value","$v >> $STATSFILEPATH
    done
done

python3 scripts/make_charts.py $STATSFILEPATH "$dataset_name" $OUTPUTDIR"/chart.png" $OUTPUTDIR"/data.csv"

deactivate
