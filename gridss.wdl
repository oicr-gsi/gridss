version 1.0

struct GenomeResources {
    String svprepModules
    String gridssModules
    String gatkModules
    String refFasta
    String refFai
    String ensembldata
    String knownfusion
    String blocklist
    Int largest
}

workflow gridss {
  input {
    File tumorBam
    File tumorBai
    File normBam
    File normBai
    Int assemblyChunks = 4
    String outputFileNamePrefix 
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
      "gridssModules": "gridss/2.13.2m hmftools-data/53138 hg38-gridss-index/1.0",
      "gatkModules": "hg38-gridss-index/1.0 gatk/4.1.6.0",
      "refFasta": "$HG38_GRIDSS_INDEX_ROOT/hg38_random.fa",
      "refFai": "$HG38_GRIDSS_INDEX_ROOT/hg38_random.fa.fai",
      "ensembldata": "$HMFTOOLS_DATA_ROOT/ensembl_data",
      "knownfusion": "$HMFTOOLS_DATA_ROOT/sv/known_fusions.38.bedpe",
      "blocklist": "$HMFTOOLS_DATA_ROOT/sv/gridss_blacklist.38.bed.gz",
      "largest": "248956422"
    }
  }

  call extractName as extractTumorName {
    input:
    refFasta = resources [ genomeVersion ].refFasta,
    refFai = resources [ genomeVersion ].refFai,
    modules = resources [ genomeVersion ].gatkModules,
    inputBam = tumorBam,
    inputBai = tumorBai
  }

  call extractName as extractNormalName {
    input:
    refFasta = resources [ genomeVersion ].refFasta,
    refFai = resources [ genomeVersion ].refFai,
    modules = resources [ genomeVersion ].gatkModules,
    inputBam = normBam,
    inputBai = normBai
  }

  call splitFaiToArray {
    input:
    modules = resources [ genomeVersion ].gatkModules,
    refFai = resources [ genomeVersion ].refFai
  }

  scatter(c in splitFaiToArray.out) {
    call getChrCoefficient as splitForTumor {
      input: 
      modules = resources [ genomeVersion ].gatkModules,
      refFai = resources [ genomeVersion ].refFai,
      largestChrom = resources [ genomeVersion ].largest,
      chromosome = c
    }
    call svprep as prepTumor {                                                                 
      input:
      chrom = c,
      scaleCoefficient = splitForTumor.coeff,
      inputname = extractTumorName.input_name,
      inputBam = tumorBam, 
      inputBai = tumorBai,
      modules = resources [ genomeVersion ].svprepModules
    }
  }

  # Collector jobs for Normal and tumor
  call aggregateData as mergeTumorData {
    input:
    inputname = extractTumorName.input_name,
    prepBamfiles = prepTumor.prepd_bam,
    prepJunctions = prepTumor.prepd_junctions
  }
 
  scatter(c in splitFaiToArray.out) {
    call getChrCoefficient as splitForNormal {
      input:
      modules = resources [ genomeVersion ].gatkModules,
      refFai = resources [ genomeVersion ].refFai,
      largestChrom = resources [ genomeVersion ].largest,
      chromosome = c
    }
    call svprep as prepNormal {
      input:
      chrom = c,
      scaleCoefficient = splitForNormal.coeff,
      inputname = extractNormalName.input_name,
      inputBam = normBam,
      inputBai = normBai,
      junctions = mergeTumorData.prepd_merged_junctions,
      modules = resources [ genomeVersion ].svprepModules
    }
  }

  call aggregateData as mergeNormalData {
    input: 
    inputname = extractNormalName.input_name,
    prepBamfiles = prepNormal.prepd_bam,
    prepJunctions = prepNormal.prepd_junctions
  }

  call preprocessInputs as preprocessNormal {
    input:
      samplename = extractNormalName.input_name,
      inputBam = mergeNormalData.prepd_merged_bam,
      modules = resources [ genomeVersion ].gridssModules,
      blocklist = resources [ genomeVersion ].blocklist,
      refFasta = resources [ genomeVersion ].refFasta
  }

  call preprocessInputs as preprocessTumor {
    input:
      samplename = extractTumorName.input_name,
      inputBam = mergeTumorData.prepd_merged_bam,
      modules = resources [ genomeVersion ].gridssModules,
      blocklist = resources [ genomeVersion ].blocklist,
      refFasta = resources [ genomeVersion ].refFasta
  }

  scatter (i in range(assemblyChunks)) {
    call assembleBam {
      input:
        normBam = mergeNormalData.prepd_merged_bam,
        tumorBam = mergeTumorData.prepd_merged_bam,
        normalName = extractNormalName.input_name,
        tumorName = extractTumorName.input_name,
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
      normalName = extractNormalName.input_name,
      tumorName = extractTumorName.input_name,
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

# =====================================================================================
#  Borrowed from cfMedipsQc workflow. Build array of chromosome names given a .fai file
# =====================================================================================
task splitFaiToArray {
  input {
    Int memory = 1
    Int timeout = 1
    String modules
    String refFai
  }

  parameter_meta {
    refFai: ".fai file for the reference genome, we use it to extract chromosome ids"
    timeout: "Hours before task timeout"
    memory: "Memory allocated for this job"
    modules: "Names and versions of modules to load"
  }

  command <<<
    cut -f 1 ~{refFai} | uniq | grep -v _ | grep -v M | grep ^chr
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    Array[String] out = read_lines(stdout())
  }

  meta {
    output_meta: {
      out: "Chromosomes to split extractMedipsCounts jobs by, in Array[String] format."
    }
  }
}

# ================================================================
#  Scaling coefficient - use to scale RAM allocation by chromosome
# ================================================================
task getChrCoefficient {
  input {
    Int memory = 1
    Int timeout = 1
    Int largestChrom
    String chromosome
    String modules
    String refFai
  }

  parameter_meta {
    refFai: ".fai file for the reference genome, we use it to extract chromosome ids"
    timeout: "Hours before task timeout"
    chromosome: "Chromosome to check"
    memory: "Memory allocated for this job"
    modules: "Names and versions of modules to load"
    largestChrom: "Length of the largest chromosome in a genome"
  }

  command <<<
    grep -w ^~{chromosome} ~{refFai} | cut -f 2 | awk '{print int(($1/~{largestChrom} + 0.1) * 10)/10}'
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    String coeff = read_string(stdout())
  }

  meta {
    output_meta: {
      coeff: "Length ratio as relative to the largest chromosome."
    }
  }
}

# =========================================
# Job to extract names from input bam files
# =========================================
task extractName {
  input {
    String modules
    String refFasta 
    String refFai 
    File inputBam
    File inputBai
    Int memory = 4
    Int timeout = 4
  }

  parameter_meta {
    inputBam: "input .bam file"
    inputBai: "input .bai file"
    refFasta: "Reference FASTA file"
    refFai: "Reference fai index"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail

    if [ -f "~{inputBam}" ]; then
      gatk --java-options "-Xmx1g" GetSampleName -R ~{refFasta} -I ~{inputBam} -O input_name.txt -encode
    fi

    cat input_name.txt
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  meta {
    output_meta: {
      input_name: "name of the input"
    }
  }

  output {
    String input_name = read_string(stdout()) 
  }
}

# ====================================
# Job to filter bams before preprocess
# ====================================

task svprep {
  input {
    String inputname
    File inputBam
    File inputBai
    File? junctions
    Int memory = 86
    Int timeout = 30
    Int threads = 4
    Int overhead = 6
    String blocklist = "$HMFTOOLS_DATA_ROOT/sv/gridss_blacklist.38.bed.gz"
    String refFasta = "$HG38_ROOT/hg38_random.fa"
    String refFastaVersion = "38"
    String svprepScript = "java  -Xmx~{round(memory * scaleCoefficient) - overhead}G -jar $HMFTOOLS_ROOT/svprep.jar"
    String ensembldata = "$HMFTOOLS_DATA_ROOT/ensembl_data"
    String knownfusion = "$HMFTOOLS_DATA_ROOT/sv/known_fusions.38.bedpe"
    String workingDir = "~{basename(inputBam)}.gridss.working"
    String chrom
    String modules
    Float scaleCoefficient
    String? additionalParameters
    Int? supportFragCap
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
    chrom: "This is a task meant to run in a scattered way, specify chromosome with chrom"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    overhead: "Overhead for java (GB)"
    timeout: "Hours before task timeout"
    scaleCoefficient: "Scaling RAM by the size of chromosome"
    additionalParameters: "Any additional parameters to svprep we want to pass"
    partition: "Partition size"
    supportFragCap: "Support frag cap, limit supporting reads per junction"
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
      -threads ~{threads} ~{additionalParameters}\
      -specific_chr ~{chrom} ~{"-existing_junction_file " + junctions} \
      -partition_size ~{partition} ~{"-junction_frags_cap " + supportFragCap} \
      -apply_downsampling

   samtools sort -T ~{workingDir} --reference ~{refFasta} ~{workingDir}/~{inputname}.sv_prep.bam -o ~{workingDir}/~{inputname}.~{chrom}.sv_prep.sorted.bam
   cp ~{workingDir}/~{inputname}.sv_prep.junctions.csv ~{workingDir}/~{inputname}.~{chrom}.sv_prep.junctions.csv
   rm ~{workingDir}/~{inputname}.sv_prep.bam
  >>>

  runtime {
    memory:  "~{round(memory * scaleCoefficient)} GB"
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
    File prepd_bam = "~{workingDir}/~{inputname}.~{chrom}.sv_prep.sorted.bam"
    File prepd_junctions  = "~{workingDir}/~{inputname}.~{chrom}.sv_prep.junctions.csv"
  }
}

# ======================================================
#  This task is for collecting results from scatter jobs
# ======================================================
task aggregateData {
  input {
    Array[File] prepBamfiles
    Array[File] prepJunctions
    String inputname
    Int jobMemory = 1
    Int timeout = 1
    String modules = "samtools/1.14"
  }

  parameter_meta {
    prepBamfiles: "bam files from svprep scatter job"
    prepJunctions: "junctions files from svprep scatter job"
    inputname: "Prefix to use for file name"
    timeout: "Hours before task timeout"
    jobMemory: "Memory allocated for this job"
    modules: "modules required for the task"
  }

  command <<<
    set -euo pipefail
    samtools merge -o ~{inputname}.sv_prep.sorted.bam ~{sep=" " prepBamfiles}
    cat ~{sep=" " prepJunctions} | sort -V -u > ~{inputname}.sv_prep.junctions.csv
  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  meta {
    output_meta: {
      prepd_merged_bam: "merged processed bam file",
      prepd_merged_junctions: "merged junctions .csv file"
    }
  }

  output {
    File prepd_merged_bam = "~{inputname}.sv_prep.sorted.bam"
    File prepd_merged_junctions  = "~{inputname}.sv_prep.junctions.csv"
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
    String gridssScript = "gridss --jar $GRIDSS_ROOT/bin/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    String workingDir = "~{basename(inputBam)}.gridss.working"
    Int memory = 16
    Int timeout = 12
    Int threads = 4
    String? additionalParameters
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
    additionalParameters: "Any additional parameters to svprep we want to pass"
  }

  command <<<
   ~{gridssScript} ~{"-b" + blocklist} \
   -r ~{refFasta} \
   -s preprocess \
   -t ~{threads} ~{additionalParameters}\
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
    String gridssScript = "gridss --jar $GRIDSS_ROOT/bin/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    String workingDirNorm = "~{basename(normBam)}.gridss.working"
    String workingDirTumr = "~{basename(tumorBam)}.gridss.working"
    Int jobNodes = 1
    Int jobIndex = 0
    Int memory = 32
    Int timeout = 24
    Int threads = 8
    String? additionalParameters
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
    additionalParameters: "Any additional parameters to svprep we want to pass"
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
    -t ~{threads} ~{additionalParameters}\
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
    String gridssScript = "gridss --jar $GRIDSS_ROOT/bin/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    String outputFileNamePrefix
    Int threads = 8
    Int memory = 50
    Int overhead = 8
    Int timeout = 24
    String? additionalParameters
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
    additionalParameters: "Any additional parameters to svprep we want to pass"
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
    -t ~{threads} ~{additionalParameters} \
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
