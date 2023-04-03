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
    tumorBam: "Input tumor file (bam)"
    tumorBai: "Input tumor file index (bai)"
    normBam: "Input normal file (bam)"
    normBai: "Input normal file index (bai)"
    assemblyChunks: "How many chunks to use for assembly job, may be calculated based on input size"
    outputFileNamePrefix: "Output file prefix"
  }

  Map[String,GenomeResources] resources = {
    "38": {
      "svprepModules": "hmftools/1.1 hmftools-data/53138 hg38-gridss-index/1.0 samtools/1.14",
      "gridssModules": "gridss/2.13.2 hmftools-data/53138 hg38-gridss-index/1.0",
      "refFasta": "$HG38_GRIDSS_INDEX_ROOT/hg38_random.fa",
      "ensembldata": "$HMFTOOLS_DATA_ROOT/ensembl_data",
      "knownfusion": "$HMFTOOLS_DATA_ROOT/sv/known_fusions.38.bedpe",
      "blocklist": "$HMFTOOLS_DATA_ROOT/sv/gridss_blacklist.38.bed.gz"
    }
  }

  call svprep {
    input:
    tumorname = tumorName,
    normalname = normalName,
    tumorBam = tumorBam,
    tumorBai = tumorBai,
    normBam = normBam,
    normBai = normBai,
    refFastaVersion = genomeVersion,
    modules = resources [ genomeVersion ].svprepModules,
    ensembldata = resources [ genomeVersion ].ensembldata,
    knownfusion = resources [ genomeVersion ].knownfusion,
    blocklist = resources [ genomeVersion ].blocklist,
    refFasta = resources [ genomeVersion ].refFasta
  }
  
  call preprocessInputs as preprocessNormal {
    input:
      samplename = normalName,
      inputBam = svprep.prepd_normal,
      modules = resources [ genomeVersion ].gridssModules,
      blocklist = resources [ genomeVersion ].blocklist,
      refFasta = resources [ genomeVersion ].refFasta
  }

  call preprocessInputs as preprocessTumor {
    input:
      samplename = tumorName,
      inputBam = svprep.prepd_tumor,
      modules = resources [ genomeVersion ].gridssModules,
      blocklist = resources [ genomeVersion ].blocklist,
      refFasta = resources [ genomeVersion ].refFasta
  }

  scatter (i in range(assemblyChunks)) {
    call assembleBam {
      input:
        normBam = svprep.prepd_normal,
        tumorBam = svprep.prepd_tumor,
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
      },
      {
        name: "samtools"
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
    String tumorname
    File tumorBam
    File tumorBai
    String normalname
    File normBam
    File normBai
    String refFasta 
    String refFastaVersion 
    String svprepScript = "java  -Xmx80G -jar $HMFTOOLS_ROOT/svprep.jar"
    String ensembldata 
    String knownfusion 
    String blocklist 
    String modules 
    Int memory = 80
    Int timeout = 30
    Int threads = 4
    Int partition = 10000
  }

  parameter_meta {
    tumorBam: "input tumor .bam file"
    normBam: "input normal .bai file"
    blocklist: "bed file with regions to ignore"
    refFasta: "Reference FASTA file"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
  }

  command <<<

    mkdir svprep

    ~{svprepScript}  \
      -sample ~{tumorname} \
      -bam_file ~{tumorBam} \
      -output_dir svprep/ \
      -ref_genome ~{refFasta} \
      -ref_genome_version ~{refFastaVersion} \
      -blacklist_bed ~{blocklist} \
      -known_fusion_bed ~{knownfusion} \
      -partition_size ~{partition} \
      -apply_downsampling 

    ~{svprepScript}  \
      -sample ~{normalname} \
      -bam_file ~{normBam} \
      -output_dir svprep/ \
      -ref_genome ~{refFasta} \
      -ref_genome_version ~{refFastaVersion} \
      -blacklist_bed  ~{blocklist} \
      -known_fusion_bed ~{knownfusion} \
      -existing_junction_file svprep/~{tumorname}.sv_prep.junctions.csv \
      -partition_size ~{partition} \
      -apply_downsampling 

      samtools sort svprep/~{normalname}.sv_prep.bam >svprep/~{normalname}.sv_prep.sort.bam
      samtools sort svprep/~{tumorname}.sv_prep.bam >svprep/~{tumorname}.sv_prep.sort.bam


  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  meta {
    output_meta: {
      prepd_normal: "processed bam file",
      prepd_tumor: "index of the processed bam"
    }
  }

  output {
    File prepd_normal = "svprep/~{normalname}.sv_prep.sort.bam"
    File prepd_tumor = "svprep/~{tumorname}.sv_prep.sort.bam"
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
    File preprocessedBam = "~{workingDir}/~{samplename}.sv_prep.sort.bam.sv.bam"
    File preprocessedIdx = "~{workingDir}/~{samplename}.sv_prep.sort.bam.sv.bam.csi"
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
    String workingDirNorm = "~{basename(processedNormBam)}.gridss.working"
    String workingDirTumr = "~{basename(processedTumrBam)}.gridss.working"
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
    ln -s ~{processedNormBam} -t ~{workingDirNorm}
    ln -s ~{processedNormCsi} -t ~{workingDirNorm}
    ln -s ~{processedTumrBam} -t ~{workingDirTumr}
    ln -s ~{processedTumrCsi} -t ~{workingDirTumr}

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
    ln -s ~{processedNormBam} -t ~{workingDirNorm}
    ln -s ~{processedNormCsi} -t ~{workingDirNorm}
    ln -s ~{processedTumrBam} -t ~{workingDirTumr}
    ln -s ~{processedTumrCsi} -t ~{workingDirTumr}
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
