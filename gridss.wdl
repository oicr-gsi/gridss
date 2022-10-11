version 1.0

workflow gridss {
  input {
    File tumorBam
    File tumorBai
    File normBam
    File normBai
    String normName = basename("~{normBam}", ".filter.deduped.realigned.recalibrated.bam")
    String tumorName = basename("~{tumorBam}", ".filter.deduped.realigned.recalibrated.bam")
  }

  parameter_meta {
    tumorBam: "Input tumor file (bam)"
    tumorBai: "Input tumor file index (bai)"
    normBam: "Input normal file (bam)"
    normBai: "Input normal file index (bai)"
  }

  call call_SVs {
    input:
      tumorBam = tumorBam,
      normBam = normBam,
      tumorBai = tumorBai,
      normBai = normBai,
      normName = normName,
      tumorName = tumorName
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
      File structuralVcf = "~{tumorName}.gridss.working/~{tumorName}.allocated.vcf"
  }
}

task call_SVs {
  input {
    File normBam
    File normBai
    File tumorBam
    File tumorBai
    String normName
    String tumorName
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    String refFasta = "$HMFTOOLS_DATA_ROOT/hg38_random.fa"
    String gridssScript = "$GRIDSS_ROOT/gridss --jar $GRIDSS_ROOT/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    Int threads = 8
    Int memory = 50
    Int timeout = 100
  }

  command <<<
    set -euo pipefail

    #mkdir ~{tumorName}

    ~{gridssScript} \
      --reference ~{refFasta} \
      #--output ~{tumorName} \
      --output ~{tumorName}.gridss.working/~{tumorName}.allocated.vcf \
      ~{normBam} ~{tumorBam}

    #mv ~{tumorName}/*.vcf ~{tumorName}.allocated.vcf

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File structuralVcf = "~{tumorName}.gridss.working/~{tumorName}.allocated.vcf"
  }
}