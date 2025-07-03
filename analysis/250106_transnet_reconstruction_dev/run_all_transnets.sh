#for S in 3 5 7 9; do
for S in 9; do
#    for M in 1e-4 1e-5 1e-6; do
    for M in 1e-6; do
#    	Rscript run_transnet_juniper.R -i input_data         -m ${M} -M -s ${S} -o rsv_m${M}fix_s${S}_incvcfs
#    	Rscript run_transnet_juniper.R -i input_data         -m ${M}    -s ${S} -o rsv_m${M}var_s${S}_incvcfs
#    	Rscript run_transnet_juniper.R -i input_data_noivars -m ${M} -M -s ${S} -o rsv_m${M}fix_s${S}_novcfs
    	Rscript run_transnet_juniper.R -i input_data_noivars -m ${M}    -s ${S} -o rsv_m${M}var_s${S}_novcfs
    done
done

