#!/bin/bash

declare -n _parms_=$1;
#------------------------------------------- INPUT
    project="${_parms_['project']}"      ### project name
scripts_dir="${_parms_['scripts_dir']}"  ### folder with scripts
  input_dir="${_parms_['input_dir']}"    ### transfer folder
 output_dir="${_parms_['output_dir']}"   ### folder with results
        log="${_parms_['log']}"          ### logfile
       cpus="${_parms_['cpus']}"         ### cpus to use
        run="${_parms_['run']}"          ### dryrun?
     dryrun="${_parms_['dryrun']}"       ### dryrun?
    options="${_parms_['options']}"      ### options for the scripts
#-------------------------------------------

source ${scripts_dir}libBash.source

local tools=/exports/mm-run/tools/           ### folder with installed tools
local   zip=gzip                             ### zipping software name: gzip, pigz, ...
local    Rs=$(jq -r ".Rs" <<< $options)      ### folder with data (R1/2 files)
local   out=${output_dir}${project}/${run}/  ### folder for output

if [[ ! -d $out ]]; then mkdir $out; fi      ### create output folder if not exists

logr $log 'Excluding reads mapped to...'

IFS_="$IFS"; IFS=$' '
   local         _metadata=$(jq -c "._metadata?"               <<< $options   )
   local            method=($(jq -r ".excl.method?"            <<< $_metadata))
   local          database=($(jq -r ".excl.database?"          <<< $_metadata))
   local            genome=($(jq -r ".excl.genome?"            <<< $_metadata))
   local       nohuman_run=($(jq -r ".excl.nohuman.run?"       <<< $_metadata))
   local       nohuman_opt=($(jq -r ".excl.nohuman.opt?"       <<< $_metadata))
   local       bowtie2_run=($(jq -r ".excl.bowtie2.run?"       <<< $_metadata))
   local       bowtie2_opt=($(jq -r ".excl.bowtie2.opt?"       <<< $_metadata))
   local      samtools_run=($(jq -r ".excl.samtools.run?"      <<< $_metadata))
   local      samtools_opt=($(jq -r ".excl.samtools.opt?"      <<< $_metadata))
   local bowtie2_build_run=($(jq -r ".excl.bowtie2_build.run?" <<< $_metadata))
   local bowtie2_build_opt=($(jq -r ".excl.bowtie2_build.opt?" <<< $_metadata))
IFS="$IFS_"

local suffix=_R1.fastq.gz

if [[ ! -f ${genome}.1.bt2 ]]; then python3 $bowtie2_build_run ${genome} ${genome}; fi

IFS_="$IFS"; IFS=$'\n'
for _sample_ in $(jq -c ".samples[]" <<< $options); do

   local sample_id=$(jq -r ".sample_id" <<< $_sample_)

   local       R1s=($(get_string ".R1s[]" "$_sample_"));
   local       R2s=($(get_string ".R2s[]" "$_sample_"));

   local       log=${out}${sample_id}.log
   
   if [[ -f $log ]]; then log_ '' > $log; fi ### initialize log file

   logr $log $(_ver_ $run)
   logr $log '-----------------'
   logr $log 'Sample          :' $sample_id
   logr $log 'Input           :' $Rs
   logr $log 'Output          :' $out
   logr $log 'Method          :' $method
   logr $log 'Database        :' $database
   logr $log 'Genome          :' $genome
   logr $log 'nohuman run     :' $nohuman_run
   logr $log 'nohuman opt     :' ${nohuman_opt[@]}
   logr $log 'bowtie2 run     :' $bowtie2_run
   logr $log 'bowtie2 opt     :' ${bowtie2_opt[@]}
   logr $log 'samtools run    :' $samtools_run
   logr $log 'bowtie2build run:' $bowtie2_build_run

   for i in ${!R1s[*]}
   do
      local R1=${R1s[$i]}
      local R2=${R2s[$i]}

      local  mode=SE
      local reads=()
      
      if [[   $R1 ]] && [[   $R2 ]]; then  reads=('-1' ${Rs}${R1}   '-2' ${Rs}${R2})  ; mode=PE; fi
      if [[   $R1 ]] && [[ ! $R2 ]]; then  reads=('-U' ${Rs}${R1})                    ; mode=SE; fi
      if [[ ! $R1 ]] && [[   $R2 ]]; then  reads=('-U' ${Rs}${R2})                    ; mode=SE; fi

      if [[   $R1 ]] && [[   $R2 ]]; then reads_=('-1' ${out}${R1}_ '-2' ${out}${R2}_); fi
      if [[   $R1 ]] && [[ ! $R2 ]]; then reads_=('-o' ${out}${R1}_)                  ; fi
      if [[ ! $R1 ]] && [[   $R2 ]]; then reads_=('-o' ${out}${R2}_)                  ; fi

      logr $log '              R1:' $R1
      logr $log '              R2:' $R2
      logr $log 'Mode            :' $mode

      if [[ ${dryrun} > 0 ]]; then continue; fi
      
      if [[ ${method} == 'nohuman' ]]
      then
         if [[ ! -f ${database} ]]; then ${nonhuman_run} --db ${database} --download; fi

         PATH=$PATH:/exports/mm-run/tools/kraken2-2.1.1
         ${nohuman_run} --threads ${cpus} --db ${database} ${nohuman_opt[@]} \
                        -o ${out}${R1} -O ${out}${R2} \
                           ${Rs}${R1}     ${Rs}${R2}
         continue
      fi

      local name=${R1%"$suffix"}_${gname}
      local  SAM=${out}${name}.sam
      local  BAM=${out}${name}.bam
      local BAMS=${out}${name}.bams

      if [[ ! -f ${genome}.1.bt2 ]]; then ${bowtie2_build_run} $genome $genome; fi

      $bowtie2_run  ${bowtie2_opt[@]} -p $cpus -x $genome ${reads[@]} -S $SAM   ### mapping
      $samtools_run sort              -@ $cpus    $SAM -o $BAMS                 ### ${SAM} -> ${BAMS} sorted BAM

      _rm_ $SAM

      ${samtools_run} index  -@ $cpus  $BAMS
      ${samtools_run} bam2fq -f 12 -@ $cpus $BAMS ${reads_[@]}; fi              ### unmapped -f 12 (clean); mapped -F 12 (human)

      _rm_ $BAMS ${BAMS}.bai

      for r in $R1 $R2
      do
         $zip -c ${out}${r}_ > ${out}${r} 
         chmod 644             ${out}${r}
         _rm_ ${out}${r}_
      done
   done
   logr $log 'Done: ' $sample_id
done

IFS="$IFS_"
