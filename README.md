# gridss

performs somatic genomic rearrangement detection and classification

## Overview

## Dependencies

* [GRIDSS](https://github.com/PapenfussLab/gridss)
* [hmftools](https://github.com/hartwigmedical/hmftools)


## Usage

### Cromwell
```
java -jar cromwell.jar run gridss.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`tumorName`|String|Name of tumor input, will be used in the output vcf header
`normalName`|String|Name of the normal input, will be used in the output vcf header
`tumorBam`|File|Input tumor file (bam)
`tumorBai`|File|Input tumor file index (bai)
`normBam`|File|Input normal file (bam)
`normBai`|File|Input normal file index (bai)


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`assemblyChunks`|Int|4|How many chunks to use for assembly job, may be calculated based on input size
`outputFileNamePrefix`|String|basename("~{tumorBam}",".filter.deduped.realigned.recalibrated.bam")|Output file prefix
`genomeVersion`|String|"38"|Genome identification, will be used in the output vcf header


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`prepTumor.junctions`|File?|None|Optional exisiting junctions file
`prepTumor.blocklist`|String|"$HMFTOOLS_DATA_ROOT/sv/gridss_blacklist.38.bed.gz"|bed file with regions to ignore
`prepTumor.refFasta`|String|"$HG38_ROOT/hg38_random.fa"|Reference FASTA file
`prepTumor.refFastaVersion`|String|"38"|version of fasta file, initially it is 38 only
`prepTumor.svprepScript`|String|"java  -Xmx80G -jar $HMFTOOLS_ROOT/svprep.jar"|path to the pre-processing script
`prepTumor.ensembldata`|String|"$HMFTOOLS_DATA_ROOT/ensembl_data"|path to ENSEMBL data (hmftools resource)
`prepTumor.knownfusion`|String|"$HMFTOOLS_DATA_ROOT/sv/known_fusions.38.bedpe"|path to known fusions (hmftools resource)
`prepTumor.workingDir`|String|"~{basename(inputBam)}.gridss.working"|Working directory
`prepTumor.memory`|Int|86|Memory allocated for this job (GB)
`prepTumor.timeout`|Int|30|Hours before task timeout
`prepTumor.threads`|Int|4|Requested CPU threads
`prepTumor.partition`|Int|10000|Partition size
`prepNormal.blocklist`|String|"$HMFTOOLS_DATA_ROOT/sv/gridss_blacklist.38.bed.gz"|bed file with regions to ignore
`prepNormal.refFasta`|String|"$HG38_ROOT/hg38_random.fa"|Reference FASTA file
`prepNormal.refFastaVersion`|String|"38"|version of fasta file, initially it is 38 only
`prepNormal.svprepScript`|String|"java  -Xmx80G -jar $HMFTOOLS_ROOT/svprep.jar"|path to the pre-processing script
`prepNormal.ensembldata`|String|"$HMFTOOLS_DATA_ROOT/ensembl_data"|path to ENSEMBL data (hmftools resource)
`prepNormal.knownfusion`|String|"$HMFTOOLS_DATA_ROOT/sv/known_fusions.38.bedpe"|path to known fusions (hmftools resource)
`prepNormal.workingDir`|String|"~{basename(inputBam)}.gridss.working"|Working directory
`prepNormal.memory`|Int|86|Memory allocated for this job (GB)
`prepNormal.timeout`|Int|30|Hours before task timeout
`prepNormal.threads`|Int|4|Requested CPU threads
`prepNormal.partition`|Int|10000|Partition size
`preprocessNormal.gridssScript`|String|"$GRIDSS_ROOT/gridss --jar $GRIDSS_ROOT/gridss-2.13.2-gridss-jar-with-dependencies.jar"|Script to run GRIDSS
`preprocessNormal.workingDir`|String|"~{basename(inputBam)}.gridss.working"|Working dir for processing, automatically generated by gridds
`preprocessNormal.memory`|Int|16|Memory allocated for this job (GB)
`preprocessNormal.timeout`|Int|12|Hours before task timeout
`preprocessNormal.threads`|Int|4|Requested CPU threads
`preprocessTumor.gridssScript`|String|"$GRIDSS_ROOT/gridss --jar $GRIDSS_ROOT/gridss-2.13.2-gridss-jar-with-dependencies.jar"|Script to run GRIDSS
`preprocessTumor.workingDir`|String|"~{basename(inputBam)}.gridss.working"|Working dir for processing, automatically generated by gridds
`preprocessTumor.memory`|Int|16|Memory allocated for this job (GB)
`preprocessTumor.timeout`|Int|12|Hours before task timeout
`preprocessTumor.threads`|Int|4|Requested CPU threads
`assembleBam.gridssScript`|String|"$GRIDSS_ROOT/gridss --jar $GRIDSS_ROOT/gridss-2.13.2-gridss-jar-with-dependencies.jar"|Script to run GRIDSS
`assembleBam.workingDirNorm`|String|"~{basename(normBam)}.gridss.working"|Working dir for normal bam, automatically generated by gridds
`assembleBam.workingDirTumr`|String|"~{basename(tumorBam)}.gridss.working"|Working dir for tumor bam, automatically generated by gridds
`assembleBam.memory`|Int|32|Memory allocated for this job (GB)
`assembleBam.timeout`|Int|24|Hours before task timeout
`assembleBam.threads`|Int|8|Requested CPU threads
`callSvs.workingDirNorm`|String|"~{basename(normBam)}.gridss.working"|Working dir for normal bam, automatically generated by gridds
`callSvs.workingDirTumr`|String|"~{basename(tumorBam)}.gridss.working"|Working dir for tumor bam, automatically generated by gridds
`callSvs.gridssScript`|String|"$GRIDSS_ROOT/gridss --jar $GRIDSS_ROOT/gridss-2.13.2-gridss-jar-with-dependencies.jar"|Script to run GRIDSS
`callSvs.threads`|Int|8|Requested CPU threads
`callSvs.memory`|Int|50|Memory allocated for this job (GB), split b/w high-demand chunk and everything else
`callSvs.overhead`|Int|8|Memory for other things ran by JVM
`callSvs.timeout`|Int|24|Hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`structuralVcf`|File|Structural Variant .vcf file


## Commands
 This section lists commands run by gridss workflow
 
 * Running GRIDSS
 
 perform joint SV calling on the tumour and normal. The workflow consists of three steps
 described below:
 
 
 ## SV_PREP script runs on inputs before preprocessing and assembly steps
 
 ```
     set -euo pipefail
     mkdir WORKING_DIR
 
     SVPREP_SCRIPT  \
       -sample INPUT_NAME \
       -bam_file INPUT_BAM \
       -output_dir WORKING_DIR \
       -ref_genome REF_FASTA \
       -ref_genome_version FASTA_VERSION \
       -blacklist_bed BLOCKLIST \
       -known_fusion_bed KNOWN_FUSIONS \
       -existing_junction_file JUNCTIONS \ [Optional]
       -partition_size PARTITIONS \
       -apply_downsampling
 
    samtools sort -T WORKING_DIR --reference REF_FASTA WORKING_DIR/INPUT_NAME.sv_prep.bam -o WORKING_DIR/INPUT_NAME.sv_prep.sorted.bam
    rm WORKING_DIR/INPUT_NAME.sv_prep.bam
   
 ```
 ## Preprocessing 
 
 ```
    GRIDSS_SCRIPT
    -b BLACKLIST \
    -r REF_FASTA \
    -s preprocess \
    -t THREADS \
    INPUT_BAM
 ```
 
 ## Assembly
 
 Assembly is performed in parallel on 10Mbase chunks. This step scales to 8 cores with peak memory usage of 31Gb
 and can be distributed across multiple nodes (using scatter)
 
 ```
     set -euo pipefail
     mkdir WORKDIR_NORM WORKDIR_TUMOR
     ln -s PROCESSED_NORMBAM WORKDIR_NORM/basename(NORMBAM).sv.bam
     ln -s PROCESSED_NORMIDX WORKDIR_NORM/basename(NORMBAM).sv.bam.csi
     ln -s PROCESSED_TUMRBAM WORKDIR_TUMOR/basename(TUMRBAM).sv.bam
     ln -s PROCESSED_TUMRIDX WORKDIR_TUMOR/basename(TUMRBAM).sv.bam.csi
 
     GRIDSS_SCRIPT
     -b BLACKLIST \
     -r REF_FASTA \
     -t THREADS \
     -s assemble \
     -a assembly.bam \
     --jobnodes JOB_NODES \
     --jobindex JOB_INDEX \
     NORM_BAM TUMOR_BAM
 ```
 
 ## Call SVs
 
 Using pre-processed data call Structural Variants
 
 ```
     set -euo pipefail
     mkdir WORKDIR_NORM WORKDIR_TUMOR
     ln -s PROCESSED_NORMBAM WORKDIR_NORM/basename(NORMBAM).sv.bam
     ln -s PROCESSED_NORMIDX WORKDIR_NORM/basename(NORMBAM).sv.bam.csi
     ln -s PROCESSED_TUMRBAM WORKDIR_TUMOR/basename(TUMRBAM).sv.bam
     ln -s PROCESSED_TUMRIDX WORKDIR_TUMOR/basename(TUMRBAM).sv.bam.csi
     mkdir assembly.bam.gridss.working
     
     for f in ~{sep=' ' assembleChunks}
     do
       ln -s $f -t assembly.bam.gridss.working
     done
 
     GRIDSS_SCRIPT
     -b BLACKLIST \
     --jvmheap (MEMORY - OVERHEAD)g \
     --otherjvmheap OVERHEADg \
     --reference REF_FASTA \
     -a assembly.bam \
     -s assemble,call \
     -t THREADS \
     -o OUTPUT_PREFIX.allocated.vcf \
     NORM_BAM TUMOR_BAM
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
