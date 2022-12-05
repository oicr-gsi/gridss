# gridss

performs somatic genomic rearrangement detection and classification

## Overview

 GRIDSS is a module software suite containing tools useful for the detection of 
 genomic rearrangements. GRIDSS includes a genome-wide break-end assembler, as 
 well as a structural variation caller for Illumina sequencing data. GRIDSS 
 calls variants based on alignment-guided positional de Bruijn graph genome-wide
 break-end assembly, split read, and read pair evidence.
 
 GRIDSS makes extensive use of the standard tags defined by SAM specifications. 
 Due to the modular design, any step (such as split read identification) can be
 replaced by another implementation that also outputs using the standard tags. 
 It is hoped that GRIDSS can serve as an exemplar modular structural variant 
 pipeline designed for interoperability with other tools.

## Dependencies

* [GRIDSS](https://github.com/PapenfussLab/gridss)


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


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|basename("~{tumorBam}",".filter.deduped.realigned.recalibrated.bam")|Output file prefix


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`call_SVs.blacklist`|String?|None|bed file with regions to ignore
`call_SVs.modules`|String|"argparser/2.1.3 stringdist/0.9.8 structuravariantannotation/1.10.1 rtracklayer/1.54.0 gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2/2.0.8 bcftools/1.9 hmftools-data/hg38"|Required environment modules
`call_SVs.refFasta`|String|"$HMFTOOLS_DATA_ROOT/hg38_random.fa"|Reference FASTA file
`call_SVs.gridssScript`|String|"$GRIDSS_ROOT/gridss --jar $GRIDSS_ROOT/gridss-2.13.2-gridss-jar-with-dependencies.jar"|Script to run GRIDSS
`call_SVs.threads`|Int|8|Requested CPU threads
`call_SVs.memory`|Int|50|Memory allocated for this job (GB)
`call_SVs.timeout`|Int|100|Hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`structuralVcf`|File|Structural Variant .vcf file


## Commands

 This section lists command(s) run by GRIDSS workflow
 
 ```
     set -euo pipefail
 
     GRIDSS_SCRIPT
     -b BLACKLIST (optional)
     --reference REF_FASTA
     --jobnodes THREADS
     --output OUTPUT_PREFIX.allocated.vcf \
     NORM_BAM TUMOR_BAM
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
