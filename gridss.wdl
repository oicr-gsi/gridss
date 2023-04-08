version 1.0

struct GenomeResources {
    String svprepModules
    String gridssModules
    String refFasta
    String ensembldata
    String knownfusion
    String blocklist
}

workflow gridss {
  input {
    String tumorName
    String normalName
    File tumorBam
    File tumorBai
    File normBam
    File normBai
    Int assemblyChunks = 4
    String outputFileNamePrefix = basename("~{tumorBam}", ".filter.deduped.realigned.recalibrated.bam")
    String genomeVersion = "38"
  }

  parameter_meta {
    tumorName: "Name of tumor input, will be used in the output vcf header"
    tumorBam: "Input tumor file (bam)"
    tumorBai: "Input tumor file index (bai)"
    normalName: "Name of the normal input, will be used in the output vcf header"
    normBam: "Input normal file (bam)"
    normBai: "Input normal file index (bai)"
    assemblyChunks: "How many chunks to use for assembly job, may be calculated based on input size"
    outputFileNamePrefix: "Output file prefix"
    genomeVersion: "Genome identification, will be used in the output vcf header"
  }

  Map[String,GenomeResources] resources = {
    "38": {
      "svprepModules": "hmftools/1.1 hmftools-data/53138 hg38-gridss-index/1.0",
      "gridssModules": "gridss/2.13.2 hmftools-data/53138 hg38-gridss-index/1.0",
      "refFasta": "$HG38_GRIDSS_INDEX_ROOT/hg38_random.fa",
      "ensembldata": "$HMFTOOLS_DATA_ROOT/ensembl_data",
      "knownfusion": "$HMFTOOLS_DATA_ROOT/sv/known_fusions.38.bedpe",
      "blocklist": "$HMFTOOLS_DATA_ROOT/sv/gridss_blacklist.38.bed.gz"
    }
  }


  call svprep as prepTumor {                                                                 
    input: 
    inputname = tumorName,
    inputBam = tumorBam, 
    inputBai = tumorBai,
    modules = resources [ genomeVersion ].svprepModules
  }                                                                                         


  call svprep as prepNormal {
    input:
    inputname = normalName,
    inputBam = normBam,
    inputBai = normBai,
    junctions = prepTumor.prepd_junctions,
    modules = resources [ genomeVersion ].svprepModules                                    
  }

  call preprocessInputs as preprocessNormal {
    input:
      samplename = normalName,
      inputBam = prepNormal.prepd_bam,
      modules = resources [ genomeVersion ].gridssModules,
      blocklist = resources [ genomeVersion ].blocklist,
      refFasta = resources [ genomeVersion ].refFasta
  }

  call preprocessInputs as preprocessTumor {
    input:
      samplename = tumorName,
      inputBam = prepTumor.prepd_bam,
      modules = resources [ genomeVersion ].gridssModules,
      blocklist = resources [ genomeVersion ].blocklist,
      refFasta = resources [ genomeVersion ].refFasta
  }

  scatter (i in range(assemblyChunks)) {
    call assembleBam {
      input:
        normBam = prepNormal.prepd_bam,
        tumorBam = prepTumor.prepd_bam,
        normalName = normalName,
        tumorName = tumorName,
        processedNormBam = preprocessNormal.preprocessedBam,
        processedNormCsi = preprocessNormal.preprocessedIdx,
        processedTumrBam = preprocessTumor.preprocessedBam,
        processedTumrCsi = preprocessTumor.preprocessedIdx,
        jobNodes = assemblyChunks,
        jobIndex = i,
        modules = resources [ genomeVersion ].gridssModules,
        blocklist = resources [ genomeVersion ].blocklist,
        refFasta = resources [ genomeVersion ].refFasta        
    }
  }

  call callSvs {
    input:
      normalName = normalName,
      tumorName = tumorName,
      tumorBam = tumorBam,
      normBam = normBam,
      processedNormBam = preprocessNormal.preprocessedBam,
      processedNormCsi = preprocessNormal.preprocessedIdx,
      processedTumrBam = preprocessTumor.preprocessedBam,
      processedTumrCsi = preprocessTumor.preprocessedIdx,
      assembleChunks = flatten(assembleBam.chunks),
      outputFileNamePrefix = outputFileNamePrefix,
      modules = resources [ genomeVersion ].gridssModules,
      blocklist = resources [ genomeVersion ].blocklist,
      refFasta = resources [ genomeVersion ].refFasta
  }

  meta {
    author: "Felix Beaudry and Alexander Fortuna"
    email: "fbeaudry@oicr.on.ca"
    description: "performs somatic genomic rearrangement detection and classification"
    dependencies: [
      {
        name: "GRIDSS",
        url: "https://github.com/PapenfussLab/gridss"
      },
      {
        name: "hmftools",
        url: "https://github.com/hartwigmedical/hmftools"
      }
    ]
    output_meta: {
      structuralVcf : "Structural Variant .vcf file"
    }
  }

  output {
    File structuralVcf = callSvs.structuralVcf
  }
}

# =================================
# Job to filter bams before preprocess
# =================================

task svprep {
  input {
    String inputname
    File inputBam
    File inputBai
    File? junctions
    String blocklist = "$HMFTOOLS_DATA_ROOT/sv/gridss_blacklist.38.bed.gz"
    String refFasta = "$HG38_ROOT/hg38_random.fa"
    String refFastaVersion = "38"
    String svprepScript = "java  -Xmx80G -jar $HMFTOOLS_ROOT/svprep.jar"
    String ensembldata = "$HMFTOOLS_DATA_ROOT/ensembl_data"
    String knownfusion = "$HMFTOOLS_DATA_ROOT/sv/known_fusions.38.bedpe"
    String workingDir = "~{basename(inputBam)}.gridss.working"
    String modules
    Int memory = 86
    Int timeout = 30
    Int threads = 4
    Int partition = 10000
  }

  parameter_meta {
    inputname: "input name to be used in the output"
    inputBam: "input .bam file"
    inputBai: "input .bai file"
    junctions: "Optional exisiting junctions file"
    blocklist: "bed file with regions to ignore"
    refFasta: "Reference FASTA file"
    refFastaVersion: "version of fasta file, initially it is 38 only"
    svprepScript: "path to the pre-processing script"
    ensembldata: "path to ENSEMBL data (hmftools resource)"
    knownfusion: "path to known fusions (hmftools resource)"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
    partition: "Partition size"
    workingDir: "Working directory"
  }

  command <<<
    set -euo pipefail
    mkdir ~{workingDir}

    ~{svprepScript}  \
      -sample ~{inputname} \
      -bam_file ~{inputBam} \
      -output_dir ~{workingDir}/ \
      -ref_genome ~{refFasta} \
      -ref_genome_version ~{refFastaVersion} \
      -blacklist_bed ~{blocklist} \
      -known_fusion_bed ~{knownfusion} \
      {"-existing_junction_file " + junctions} \
      -partition_size ~{partition} \
      -apply_downsampling

   samtools sort -T ~{workingDir} --reference ~{refFasta} ~{workingDir}/~{inputname}.sv_prep.bam -o ~{workingDir}/~{inputname}.sv_prep.sorted.bam
   rm ~{workingDir}/~{inputname}.sv_prep.bam
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  meta {
    output_meta: {
      prepd_bam: "processed bam file",
      prepd_junctions: "junctions .csv file"
    }
  }

  output {
    File prepd_bam = "~{workingDir}/~{inputname}.sv_prep.sorted.bam"
    File prepd_junctions  = "~{workingDir}/~{inputname}.sv_prep.junctions.csv"
  }
}

# =================================
# Job to preprocess input bam files
# =================================
task preprocessInputs {
  input {
    File inputBam
    String? blocklist
    String samplename
    String modules
    String refFasta 
    String gridssScript = "$GRIDSS_ROOT/gridss --jar $GRIDSS_ROOT/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    String workingDir = "~{basename(inputBam)}.gridss.working"
    Int memory = 16
    Int timeout = 12
    Int threads = 4
  }

  parameter_meta {
    inputBam: "input .bam file"
    blocklist: "bed file with regions to ignore"
    refFasta: "Reference FASTA file"
    gridssScript: "Script to run GRIDSS"
    workingDir: "Working dir for processing, automatically generated by gridds"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
  }

  command <<<
   ~{gridssScript} ~{"-b" + blocklist} \
   -r ~{refFasta} \
   -s preprocess \
   -t ~{threads} \
   --labels ~{samplename} \
   ~{inputBam}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  meta {
    output_meta: {
      preprocessedBam: "processed bam file",
      preprocessedIdx: "index of the processed bam"
    }
  }

  output {
    File preprocessedBam = "~{workingDir}/~{basename(inputBam)}.sv.bam"
    File preprocessedIdx = "~{workingDir}/~{basename(inputBam)}.sv.bam.csi"
  }
}

# ======================================
# Assembly task, for speed we scatter it
# ======================================

task assembleBam {
  input {
    File tumorBam
    File normBam
    File processedNormBam
    File processedNormCsi
    File processedTumrBam
    File processedTumrCsi
    String? blocklist
    String normalName
    String tumorName
    String modules 
    String refFasta 
    String gridssScript = "$GRIDSS_ROOT/gridss --jar $GRIDSS_ROOT/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    String workingDirNorm = "~{basename(normBam)}.gridss.working"
    String workingDirTumr = "~{basename(tumorBam)}.gridss.working"
    Int jobNodes = 1
    Int jobIndex = 0
    Int memory = 32
    Int timeout = 24
    Int threads = 8
  }
 
  parameter_meta {
    normBam: "normal input .bam file"
    tumorBam: "tumor input .bam file"
    processedNormBam: "processed normal input .bam file"
    processedNormCsi: "processed normal input .csi file"
    processedTumrBam: "processed tumor input .bam file"
    processedTumrCsi: "processed tumor input .bam file"
    blocklist: "bed file with regions to ignore"
    refFasta: "Reference FASTA file"
    gridssScript: "Script to run GRIDSS"
    workingDirNorm: "Working dir for normal bam, automatically generated by gridds"
    workingDirTumr: "Working dir for tumor bam, automatically generated by gridds"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail
    mkdir ~{workingDirNorm} ~{workingDirTumr}
    ln -s ~{processedNormBam} ~{workingDirNorm}/~{basename(normBam)}.sv.bam
    ln -s ~{processedNormCsi} ~{workingDirNorm}/~{basename(normBam)}.sv.bam.csi
    ln -s ~{processedTumrBam} ~{workingDirTumr}/~{basename(tumorBam)}.sv.bam
    ln -s ~{processedTumrCsi} ~{workingDirTumr}/~{basename(tumorBam)}.sv.bam.csi

    ~{gridssScript} ~{"-b" + blocklist} \
    -r ~{refFasta} \
    -t ~{threads} \
    -s assemble \
    -a assembly.bam \
    --jobnodes ~{jobNodes} \
    --jobindex ~{jobIndex} \
    --labels ~{normalName},~{tumorName} \
    ~{normBam} ~{tumorBam} 

  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  meta {
    output_meta: {
      chunks: "Assembly chunks, used in the final step sfter this one"
    }
  }

  output {
    Array[File] chunks = glob("assembly.bam.gridss.working/*.ba*")
  }
}


# ===============================
# Final task for calling variants
# ===============================
task callSvs {
  input {
    File normBam
    File tumorBam
    File processedNormBam
    File processedNormCsi
    File processedTumrBam
    File processedTumrCsi
    String normalName
    String tumorName
    String workingDirNorm = "~{basename(normBam)}.gridss.working"
    String workingDirTumr = "~{basename(tumorBam)}.gridss.working"
    Array[File] assembleChunks
    String? blocklist
    String modules 
    String refFasta 
    String gridssScript = "$GRIDSS_ROOT/gridss --jar $GRIDSS_ROOT/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    String outputFileNamePrefix
    Int threads = 8
    Int memory = 50
    Int overhead = 8
    Int timeout = 24
  }

  parameter_meta {
    normBam: "normal input .bam file"
    tumorBam: "tumor input .bam file"
    processedNormBam: "processed normal input .bam file"
    processedNormCsi: "processed normal input .csi file"
    processedTumrBam: "processed tumor input .bam file"
    processedTumrCsi: "processed tumor input .bam file"
    workingDirNorm: "Working dir for normal bam, automatically generated by gridds"
    workingDirTumr: "Working dir for tumor bam, automatically generated by gridds"
    blocklist: "bed file with regions to ignore"
    refFasta: "Reference FASTA file"
    gridssScript: "Script to run GRIDSS"
    outputFileNamePrefix: "Output prefix for vcf file"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB), split b/w high-demand chunk and everything else"
    overhead: "Memory for other things ran by JVM"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail
    mkdir ~{workingDirNorm} ~{workingDirTumr}
    ln -s ~{processedNormBam} ~{workingDirNorm}/~{basename(normBam)}.sv.bam
    ln -s ~{processedNormCsi} ~{workingDirNorm}/~{basename(normBam)}.sv.bam.csi
    ln -s ~{processedTumrBam} ~{workingDirTumr}/~{basename(tumorBam)}.sv.bam
    ln -s ~{processedTumrCsi} ~{workingDirTumr}/~{basename(tumorBam)}.sv.bam.csi

    mkdir assembly.bam.gridss.working
    for f in ~{sep=' ' assembleChunks}
    do
      ln -s $f -t assembly.bam.gridss.working
    done

    ~{gridssScript} ~{"-b" + blocklist} \
    --jvmheap ~{memory - overhead}g \
    --otherjvmheap ~{overhead}g \
    --reference ~{refFasta} \
    -a assembly.bam \
    -s assemble,call \
    -t ~{threads} \
    -o ~{outputFileNamePrefix}.allocated.vcf \
    --labels ~{normalName},~{tumorName} \
    ~{normBam} ~{tumorBam}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  meta {
    output_meta: {
      structuralVcf: "Final vcf with structural variants"
    }
  }

  output {
    File structuralVcf = "~{outputFileNamePrefix}.allocated.vcf"
  }
}
