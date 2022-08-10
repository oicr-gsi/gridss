version 1.0

workflow hmftools {
  input {
    File tumorBam
    File tumorBai
    File? normalBam
    File? normalBai
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
      name: "hmftools",
      url: ""
    },
    {
      name: "GRIDSS",
      url: ""
    }]
  }

  call preprocess as preprocessnormal {
    input:
      sampleBam = normalBam
  }

  call preprocess as preprocesstumor {
    input:
      sampleBam = tumorBam
  }

  call SVcall {
    input:
      tumorBam = preprocesstumor.sampleBam
      normalBam = preprocessnormal.sampleBam
  }

  call cobalt {

  }

  call amber {
  
  }

  call purple {
  
  }

  output {
  }
}

task preprocess {
  input {
    String modules = "gridss/2.13.2 hg38/p12 hmftools/1.0"
    String refFasta = "$HG19_ROOT/hg19_random.fa"
    String refFai = "$HG19_ROOT/hg19_random.fa.fai"
    String refDict = "$HG19_ROOT/hg19_random.dict"
    File sampleBam
    File sampleBai
    Int threads = 4
    Int memory = 32
    Int timeout = 24
  }

  command <<<
    set -euo pipefail

    ~{GRIDSS_ROOT}/gridss --jar ~{GRIDSS_ROOT}/gridss-2.13.2-gridss-jar-with-dependencies.jar \
    --reference ~{GRIDSS_RESOURCE_ROOT}/hg38_random.fa --steps preprocess --output ~{sampleBam}.gridss ~{sampleBam}

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File unfilteredVcf = "~{outputVcf}"
  }
}


task SVcall {
  input {
    String modules = "gatk/4.1.6.0 hg19/p13"
    String refFasta = "$HG19_ROOT/hg19_random.fa"
    String refFai = "$HG19_ROOT/hg19_random.fa.fai"
    String refDict = "$HG19_ROOT/hg19_random.dict"
    File tumorBam
    File tumorBai
    File normalBam
    File normalBai
    Int threads = 4
    Int memory = 32
    Int timeout = 24
  }

  command <<<
    set -euo pipefail

  ~{GRIDSS_ROOT}/gridss --jar ~{GRIDSS_ROOT}/gridss-2.13.2-gridss-jar-with-dependencies.jar -t 1 \
  --steps assemble,call \
  --reference ~{GRIDSS_RESOURCE_ROOT}/hg38_random.fa \
  --output ~{tumorBam}.gridss ~{normalBam} ~{tumorBam} 

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File unfilteredVcf = "~{outputVcf}"
  }
}
