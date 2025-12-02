#!/bin/bash

declare -n _parms_=$1;
#--------------------------------------------- INPUT
      project="${_parms_['project']}"      ### project name
  scripts_dir="${_parms_['scripts_dir']}"  ### folder with scripts
    input_dir="${_parms_['input_dir']}"    ### transfer folder
   output_dir="${_parms_['output_dir']}"   ### folder with results
         cpus="${_parms_['cpus']}"         ### CPUs to use
          log="${_parms_['log']}"          ### logfile
          zip="${_parms_['zip']}"          ### zipped files?
          run="${_parms_['run']}"          ### tool to run
       dryrun="${_parms_['dryrun']}"       ### dryrun?
      options="${_parms_['options']}"      ### options for the scripts
#---------------------------------------------

source ${scripts_dir}libBash.source

local  in=${input_dir}/                    ### folder with the sequencing data
local  Rs=$(jq -r ".Rs" <<< $options);     ### folder with data (R1/2 files)
local out=${output_dir}${project}/${run}/  ### folder for output

if [[ ! -d $out ]]; then mkdir $out; fi    ### create output folder if not exists

logr $log 'Filtering reads...'

IFS_="$IFS"; IFS=$' '
   local         _metadata=$(jq  -c "._metadata?"             <<< $options   )
   local   trimmomatic_run=($(jq -r ".filt.trimmomatic.run?"  <<< $_metadata))
   local   trimmomatic_opt=($(jq -r ".filt.trimmomatic.opt?"  <<< $_metadata))
   local  trimmomatic_type=($(jq -r ".filt.trimmomatic.type?" <<< $_metadata))
   local        FastQC_run=($(jq -r ".filt.FastQC.run?"       <<< $_metadata))
   local        FastQC_opt=($(jq -r ".filt.FastQC.opt?"       <<< $_metadata))
   local            method=($(jq -r ".filt.method?"           <<< $_metadata))
   local      nanofilt_run=($(jq -r ".filt.nanofilt.run?"     <<< $_metadata))
   local      nanofilt_opt=($(jq -r ".filt.nanofilt.opt?"     <<< $_metadata))
   local     minion_qc_run=($(jq -r ".filt.minion_qc.run?"    <<< $_metadata))
   local     minion_qc_opt=($(jq -r ".filt.minion_qc.opt?"    <<< $_metadata))
   local      porechop_run=($(jq -r ".filt.porechop.run?"     <<< $_metadata))
   local      porechop_opt=($(jq -r ".filt.porechop.opt?"     <<< $_metadata))
   local      filtlong_run=($(jq -r ".filt.filtlong.run?"     <<< $_metadata))
   local      filtlong_opt=($(jq -r ".filt.filtlong.opt?"     <<< $_metadata))
IFS="$IFS_"

###################################### DEFAULTS ##########################################
if [[ $method           == 'null' ]]; then           method='trimmomatic'; fi
if [[ $trimmomatic_type == 'null' ]]; then trimmomatic_type=('PE'); fi
if [[ $trimmomatic_opt  == 'null' ]]; then  trimmomatic_opt=('ILLUMINACLIP:/exports/mm-run/Bioinformatics/tools/NANNY/adapters/GenomeScan.fa:2:30:7:5:true' 'LEADING:3' 'TRAILING:3' 'SLIDINGWINDOW:4:15' 'MINLEN:40'); fi
##########################################################################################

local suffix=.gz

IFS_="$IFS"; IFS=$'\n'
for _sample_ in $(jq -c ".samples[]" <<< $options); do
   local sample_id=$(jq -r ".sample_id" <<< $_sample_);

   local R1s=($(get_string ".R1s[]" "$_sample_"));
   local R2s=($(get_string ".R2s[]" "$_sample_"));

   local log=${out}${sample_id}.log 
   
   if [[ -f $log ]]; then log_ '' > $log; fi ### initialize log file
    
   logr $log $(_ver_ $run)
   logr $log '----------------'
   logr $log 'Sample         :' $sample_id
   logr $log 'Input          :' $Rs
   logr $log 'Output         :' $out
   logr $log 'FastQC run     :' $FastQC_run
   logr $log 'FastQC         :' ${FastQC_opt[@]}
   logr $log 'Method         :' $method
   logr $log 'Trimmomatic run:' $trimmomatic_run
   logr $log '            opt:' ${trimmomatic_opt[@]}
   logr $log 'Minion_qc run  :' $minion_qc_run
   logr $log '          opt  :' ${minion_qc_opt[@]}
   logr $log 'Porechop  run  :' $porechop_run
   logr $log '          opt  :' ${porechop_opt[@]}
   logr $log 'Filtlong  run  :' $filtlong_run
   logr $log '          opt  :' ${filtlong_opt[@]}
   logr $log 'Nanofilt  run  :' $nanofilt_run
   logr $log '          opt  :' ${nanofilt_opt[@]}

   for i in ${!R1s[*]}; do

      local R1=${R1s[$i]};
      local R2=${R2s[$i]};

      local R1_=${R1%"$suffix"}
      local R2_=${R2%"$suffix"}

      logr $log 'R1             :' $R1
      logr $log 'R2             :' $R2

      local trimmomatic_type=PE
      local            reads=()
      local           reads_=()

      if [[   $R1 ]] && [[   $R2 ]]; then reads=(${Rs}${R1} ${Rs}${R2}); reads_=(${out}${R1_} ${out}${R1_}_ ${out}${R2_} ${out}${R2_}_); trimmomatic_type=PE; fi
      if [[   $R1 ]] && [[ ! $R2 ]]; then reads=(${Rs}${R1}           ); reads_=(${out}${R1_}                            );              trimmomatic_type=SE; fi
      if [[ ! $R1 ]] && [[   $R2 ]]; then reads=(           ${Rs}${R2}); reads_=(                           ${out}${R2_} );              trimmomatic_type=SE; fi

      logr $log 'Mode           :' ${trimmomatic_type}

      if [[ $dryrun > 0 ]]; then continue; fi

      if [[ $method == 'trimmomatic' ]]
      then
         ${FastQC_run} -t ${cpus} -o $out ${FastQC_opt[@]} ${reads[@]}
         java -jar ${trimmomatic_run} ${trimmomatic_type} -threads ${cpus} ${reads[@]} ${reads_[@]} ${trimmomatic_opt[@]}

         for r in ${out}${R1_} ${out}${R2_}
         do 
            $zip -f ${r}
            _rm_ ${r}_
         done
      fi

      if [[ $method == 'MinKNOW' ]]
      then
         if [[ ${filtlong_run} != 'null' ]]; then ${filtlong_run} ${filtlong_opt[@]} ${Rs}${R1} | $zip > ${out}${R1_}.gz; fi
         if [[ ${porechop_run} != 'null' ]]; then ${porechop_run} -i ${Rs}${R1} -o ${out}${R1_} -t ${cpus} ${porechop_opt[@]}; $zip -f ${out}${R1_}; fi
      fi

   done
   logr $log 'Done: ' $sample_id
done

IFS="$IFS_"
