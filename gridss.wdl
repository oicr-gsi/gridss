version 1.0

workflow gridss {
  input {
    File tumorBam
    File tumorBai
    File normalBam
    File normalBai
  }

  parameter_meta {
    tumorBam: "Input tumor file (bam or sam)."
    normalBam: "Input normal file (bam or sam)."
  }

  meta {
    author: "Felix Beaudry and Alexander Fortuna"
    email: "afortuna@oicr.on.ca"
    description: "performs somatic genomic rearrangement detection and classification"
    dependencies: [
    {
      name: "GRIDSS",
      url: "https://github.com/PapenfussLab/gridss"
    }
    ]
  }

  call gridss {
    input:
      tumorBam = tumorBam
      normalBam = normalBam
      tumorBai = tumorBai
      normalBai = normalBai
  }

  output {
  }
}

task gridss {
  input {
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    String refFasta = "${HMFTOOLS_DATA_ROOT}/hg38_random.fa"
    String gridssScript = "${GRIDSS_ROOT}/gridss --jar ${GRIDSS_ROOT}/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    File normBam
    File normBai
    File tumorBam
    File tumorBai
    Int threads = 8
    Int memory = 50
    Int timeout = 100
  }

  command <<<
    set -euo pipefail

    ~{gridssScript} \
    --reference ~{refFasta} \
    --output ./ \
    ~{normBam} ~{tumorBam}

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File structuralVcf = "~{outputVcf}.allocated.vcf"
  }
}
