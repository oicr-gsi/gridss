version 1.0

workflow gridss {
  input {
    File tumorBam
    File tumorBai
    File normBam
    File normBai
    String outputFileNamePrefix = basename("~{tumorBam}", ".filter.deduped.realigned.recalibrated.bam")
  }

  parameter_meta {
    tumorBam: "Input tumor file (bam)"
    tumorBai: "Input tumor file index (bai)"
    normBam: "Input normal file (bam)"
    normBai: "Input normal file index (bai)"
    outputFileNamePrefix: "Output file prefix"
  }

  call call_SVs {
    input:
      tumorBam = tumorBam,
      normBam = normBam,
      tumorBai = tumorBai,
      normBai = normBai,
      outputFileNamePrefix = outputFileNamePrefix
  }

  meta {
    author: "Felix Beaudry and Alexander Fortuna"
    email: "fbeaudry@oicr.on.ca"
    description: "performs somatic genomic rearrangement detection and classification"
    dependencies: [
      {
        name: "GRIDSS",
        url: "https://github.com/PapenfussLab/gridss"
      }
    ]
    output_meta: {
      structuralVcf : "Structural Variant .vcf file"
    }
  }

  output {
    File structuralVcf = call_SVs.structuralVcf
  }
}

task call_SVs {
  input {
    File normBam
    File normBai
    File tumorBam
    File tumorBai
    String? blacklist
    String modules = "argparser/2.1.3 stringdist/0.9.8 structuravariantannotation/1.10.1 rtracklayer/1.54.0 gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2/2.0.8 bcftools/1.9 hmftools-data/hg38"
    String refFasta = "$HMFTOOLS_DATA_ROOT/hg38_random.fa"
    String gridssScript = "$GRIDSS_ROOT/gridss --jar $GRIDSS_ROOT/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    String outputFileNamePrefix
    Int threads = 8
    Int memory = 50
    Int timeout = 100
  }

  parameter_meta {
    normBam: "normal input .bam file"
    tumorBam: "tumor input .bam file"
    normBai: "normal input .bai file"
    tumorBai: "tumor input .bai file"
    blacklist: "bed file with regions to ignore"
    refFasta: "Reference FASTA file"
    gridssScript: "Script to run GRIDSS"
    outputFileNamePrefix: "Output prefix for vcf file"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail

    ~{gridssScript} \
    ~{"-b" + blacklist} \
    --jvmheap ~{memory - 20} \
    --reference ~{refFasta} \
    --jobnodes ~{threads} \
    --output ~{outputFileNamePrefix}.allocated.vcf \
    ~{normBam} ~{tumorBam}
  >>>

  runtime {
    cpu: "~{threads}"
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
