# gridss

performs somatic genomic rearrangement detection and classification

## Overview

GRIDSS is a module software suite containing tools useful for the detection of genomic rearrangements. It includes a genome-wide break-end assembler, as well as a structural variation caller for Illumina sequencing data. Variants are called based on alignment-guided positional de Bruijn graph genome-wide break-end assembly, split read, and read pair evidence.
GRIDSS makes extensive use of the standard tags defined by SAM specifications. Due to the modular design, any step (such as split read identification) can be replaced by another implementation that also outputs using the standard tags. It is hoped that this tool can serve as an exemplar modular structural variant pipeline designed for interoperability with other tools.

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
`tumorBam`|File|Input tumor file (bam)
`tumorBai`|File|Input tumor file index (bai)
`normBam`|File|Input normal file (bam)
`normBai`|File|Input normal file index (bai)
`outputFileNamePrefix`|String|Output file prefix


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`assemblyChunks`|Int|4|How many chunks to use for assembly job, may be calculated based on input size
`genomeVersion`|String|"38"|Genome identification, will be used in the output vcf header


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`extractTumorName.memory`|Int|4|Memory allocated for this job (GB)
`extractTumorName.timeout`|Int|4|Hours before task timeout
`extractNormalName.memory`|Int|4|Memory allocated for this job (GB)
`extractNormalName.timeout`|Int|4|Hours before task timeout
`splitFaiToArray.memory`|Int|1|Memory allocated for this job
`splitFaiToArray.timeout`|Int|1|Hours before task timeout
`splitForTumor.memory`|Int|1|Memory allocated for this job
`splitForTumor.timeout`|Int|1|Hours before task timeout
`prepTumor.junctions`|File?|None|Optional exisiting junctions file
`prepTumor.memory`|Int|86|Memory allocated for this job (GB)
`prepTumor.timeout`|Int|30|Hours before task timeout
`prepTumor.threads`|Int|4|Requested CPU threads
`prepTumor.overhead`|Int|6|Overhead for java (GB)
`prepTumor.blocklist`|String|"$HMFTOOLS_DATA_ROOT/sv/gridss_blacklist.38.bed.gz"|bed file with regions to ignore
`prepTumor.refFasta`|String|"$HG38_ROOT/hg38_random.fa"|Reference FASTA file
`prepTumor.refFastaVersion`|String|"38"|version of fasta file, initially it is 38 only
`prepTumor.svprepScript`|String|"java  -Xmx~{round((memory * scaleCoefficient)) - overhead}G -jar $HMFTOOLS_ROOT/svprep.jar"|path to the pre-processing script
`prepTumor.ensembldata`|String|"$HMFTOOLS_DATA_ROOT/ensembl_data"|path to ENSEMBL data (hmftools resource)
`prepTumor.knownfusion`|String|"$HMFTOOLS_DATA_ROOT/sv/known_fusions.38.bedpe"|path to known fusions (hmftools resource)
`prepTumor.workingDir`|String|"~{basename(inputBam)}.gridss.working"|Working directory
`prepTumor.additionalParameters`|String?|None|Any additional parameters to svprep we want to pass
`prepTumor.supportFragCap`|Int?|None|Support frag cap, limit supporting reads per junction
`prepTumor.partition`|Int|10000|Partition size
`mergeTumorData.jobMemory`|Int|1|Memory allocated for this job
`mergeTumorData.timeout`|Int|1|Hours before task timeout
`mergeTumorData.modules`|String|"samtools/1.14"|modules required for the task
`splitForNormal.memory`|Int|1|Memory allocated for this job
`splitForNormal.timeout`|Int|1|Hours before task timeout
`prepNormal.memory`|Int|86|Memory allocated for this job (GB)
`prepNormal.timeout`|Int|30|Hours before task timeout
`prepNormal.threads`|Int|4|Requested CPU threads
`prepNormal.overhead`|Int|6|Overhead for java (GB)
`prepNormal.blocklist`|String|"$HMFTOOLS_DATA_ROOT/sv/gridss_blacklist.38.bed.gz"|bed file with regions to ignore
`prepNormal.refFasta`|String|"$HG38_ROOT/hg38_random.fa"|Reference FASTA file
`prepNormal.refFastaVersion`|String|"38"|version of fasta file, initially it is 38 only
`prepNormal.svprepScript`|String|"java  -Xmx~{round((memory * scaleCoefficient)) - overhead}G -jar $HMFTOOLS_ROOT/svprep.jar"|path to the pre-processing script
`prepNormal.ensembldata`|String|"$HMFTOOLS_DATA_ROOT/ensembl_data"|path to ENSEMBL data (hmftools resource)
`prepNormal.knownfusion`|String|"$HMFTOOLS_DATA_ROOT/sv/known_fusions.38.bedpe"|path to known fusions (hmftools resource)
`prepNormal.workingDir`|String|"~{basename(inputBam)}.gridss.working"|Working directory
`prepNormal.additionalParameters`|String?|None|Any additional parameters to svprep we want to pass
`prepNormal.supportFragCap`|Int?|None|Support frag cap, limit supporting reads per junction
`prepNormal.partition`|Int|10000|Partition size
`mergeNormalData.jobMemory`|Int|1|Memory allocated for this job
`mergeNormalData.timeout`|Int|1|Hours before task timeout
`mergeNormalData.modules`|String|"samtools/1.14"|modules required for the task
`preprocessNormal.gridssScript`|String|"gridss --jar $GRIDSS_ROOT/bin/gridss-2.13.2-gridss-jar-with-dependencies.jar"|Script to run GRIDSS
`preprocessNormal.workingDir`|String|"~{basename(inputBam)}.gridss.working"|Working dir for processing, automatically generated by gridds
`preprocessNormal.memory`|Int|16|Memory allocated for this job (GB)
`preprocessNormal.timeout`|Int|12|Hours before task timeout
`preprocessNormal.threads`|Int|4|Requested CPU threads
`preprocessNormal.additionalParameters`|String?|None|Any additional parameters to svprep we want to pass
`preprocessTumor.gridssScript`|String|"gridss --jar $GRIDSS_ROOT/bin/gridss-2.13.2-gridss-jar-with-dependencies.jar"|Script to run GRIDSS
`preprocessTumor.workingDir`|String|"~{basename(inputBam)}.gridss.working"|Working dir for processing, automatically generated by gridds
`preprocessTumor.memory`|Int|16|Memory allocated for this job (GB)
`preprocessTumor.timeout`|Int|12|Hours before task timeout
`preprocessTumor.threads`|Int|4|Requested CPU threads
`preprocessTumor.additionalParameters`|String?|None|Any additional parameters to svprep we want to pass
`assembleBam.gridssScript`|String|"gridss --jar $GRIDSS_ROOT/bin/gridss-2.13.2-gridss-jar-with-dependencies.jar"|Script to run GRIDSS
`assembleBam.workingDirNorm`|String|"~{basename(normBam)}.gridss.working"|Working dir for normal bam, automatically generated by gridds
`assembleBam.workingDirTumr`|String|"~{basename(tumorBam)}.gridss.working"|Working dir for tumor bam, automatically generated by gridds
`assembleBam.memory`|Int|32|Memory allocated for this job (GB)
`assembleBam.overhead`|Int|8|A number of GB to substract from jobMemory and use as Java Heap (-Xmx)
`assembleBam.timeout`|Int|24|Hours before task timeout
`assembleBam.threads`|Int|8|Requested CPU threads
`assembleBam.additionalParameters`|String?|None|Any additional parameters to svprep we want to pass
`callSvs.workingDirNorm`|String|"~{basename(normBam)}.gridss.working"|Working dir for normal bam, automatically generated by gridds
`callSvs.workingDirTumr`|String|"~{basename(tumorBam)}.gridss.working"|Working dir for tumor bam, automatically generated by gridds
`callSvs.gridssScript`|String|"gridss --jar $GRIDSS_ROOT/bin/gridss-2.13.2-gridss-jar-with-dependencies.jar"|Script to run GRIDSS
`callSvs.threads`|Int|8|Requested CPU threads
`callSvs.memory`|Int|50|Memory allocated for this job (GB), split b/w high-demand chunk and everything else
`callSvs.overhead`|Int|8|Memory for other things ran by JVM
`callSvs.timeout`|Int|24|Hours before task timeout
`callSvs.additionalParameters`|String?|None|Any additional parameters to svprep we want to pass


### Outputs

Output | Type | Description | Labels
---|---|---|---
`structuralVcf`|File|Structural Variant .vcf file|vidarr_label: structuralVcf 


## Commands
This section lists commands run by gridss workflow
 
* Running GRIDSS
 
perform joint SV calling on the tumour and normal. The workflow consists of three steps
described below:
 
### Get all chromosomes
 
Return list of all chromosomes in assembly excluding ALT and RANDOM contigs, as well as mitochondrial chromosome
 
```
     cut -f 1 REFERENCE_FAI | uniq | grep -v _ | grep -v M | grep ^chr
```
 
### Get RAM-scaling coefficient based on chromosome size
 
This task is used to scale RAM allocation for chunked svprep task based on chromosome size
 
```
     grep -w ^CHROM REFERENCE_FAI | cut -f 2 | awk '{print int(($1/LARGEST_CHR_SIZE} + 0.1) * 10)/10}'
```
 
### Extract input name using GATK
 
Analogous to mutect2 workflow, for keeping names used in BAM alignments.
 
```
     set -euo pipefail
 
     if [ -f INPUT_BAM ]; then
       gatk --java-options "-Xmx1g" GetSampleName -R REF_FASTA -I INPUT_BAM -O input_name.txt -encode
     fi
 
     cat input_name.txt
```
 
### SV_PREP script runs on inputs before preprocessing and assembly steps
 
This is a preprocessing task which is split by chromosome. --partition_size and -junction_frags_cap affect the performance
which may be greatly enchanced by setting -junction_frags_cap to low hundreds. -partition_size is 1M by default, lowering it
increases memory requirements.
 
```
     set -euo pipefail
     mkdir WORKING_DIR
 
       SVPREP_SCRIPT  \
       -sample INPUT_NAME \
       -bam_file INPUT_BAM \
       -output_dir WORKING_DIR \
       -ref_genome REF_FASTA \
       -ref_genome_version REF_FASTA_VERSION \
       -blacklist_bed BLOCKLIST \
       -known_fusion_bed KNOWN_FUSIONS \
       -threads THREADS \
       -specific_chr CHROM ~{"-existing_junction_file " + JUNCTIONS} \
       -partition_size PARTITIONS ~{"-junction_frags_cap " + SUPPORT_FRAG_CAP} \
       -apply_downsampling
 
    samtools sort -T WORKING_DIR --reference REF_FASTA WORKING_DIR/INPUT_NAME.sv_prep.bam -o WORKING_DIR/INPUT_NAME.CHROM.sv_prep.sorted.bam
    cp WORKING_DIR/INPUT_NAME.sv_prep.junctions.csv WORKING_DIR/INPUT_NAME.CHROM.sv_prep.junctions.csv
    rm WORKING_DIR/INPUT_NAME.sv_prep.bam
```
 
### Aggregate chunked results from svprep task
 
```
     set -euo pipefail
     samtools merge -o INPUT_NAME.sv_prep.sorted.bam ~{sep=" " prepBamfiles}
     cat ~{sep=" " prepJunctions} | sort -V -u > INPUT_NAME.sv_prep.junctions.csv
```
 
### Preprocessing 
 
```
    GRIDSS_SCRIPT
    -b BLACKLIST \
    -r REF_FASTA \
    -s preprocess \
    -t THREADS \
    INPUT_BAM
```
 
### Assembly
 
Assembly is performed in parallel on 10Mbase chunks. This step scales to 8 cores with peak memory usage of 31Gb
and can be distributed across multiple nodes (using scatter)
 
```
     set -euo pipefail
     mkdir WORKING_DIR_NORM WORKING_DIR_TUMOR
     ln -s PROCESSED_NORM_BAM WORKING_DIR_NORM/BASENAME_NORM_BAM.sv.bam
     ln -s PROCESSED_NORM_CSI WORKING_DIR_NORM/BASENAME_NORM_BAM.sv.bam.csi
     ln -s PROCESSED_TUMR_BAM WORKING_DIR_TUMR/BASENAME_TUMR_BAM.sv.bam
     ln -s PROCESSED_TUMR_CSI WORKING_DIR_TUMR/BASENAME_TUMR_BAM.sv.bam.csi
 
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
 
### Call SVs
 
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
