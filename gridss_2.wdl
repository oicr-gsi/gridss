version 1.0

workflow gridss {
  input {
    File tumorBam
    File tumorBai
    File normBam
    File normBai
    String normName = basename("~{normBam}", "_ch21.bam")
    String tumorName = basename("~{tumorBam}", "_ch21.bam")
    Array[file] sampleBAMs = [tumorBam, normBam]
  }

  parameter_meta {
    tumorBam: "Input tumor file (bam)"
    tumorBai: "Input tumor file index (bai)"
    normBam: "Input normal file (bam)"
    normBai: "Input normal file index (bai)"
  }

  #call call_SVs {
    #input:
      #tumorBam = tumorBam,
      #normBam = normBam,
      #tumorBai = tumorBai,
      #normBai = normBai,
      #normName = normName,
      #tumorName = tumorName
  #}



  meta {
    author: "Felix Beaudry and Alexander Fortuna and Wen Tong"
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

  call setupreference {}

  scatter (sample in sampleBAMs) {
    call preprocess { input: bamFile=sample }
  }

  call assemble {
    input:
      tumorBam = tumorBam,
      normBam = normBam
  }

  call call {
    input:
      tumorBam = tumorBam,
      normBam = normBam
  }


  output {
      File structuralVcf = "~{tumorName}.allocated.vcf"
  }
}

task setupreference {
  
  input {
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    String gridssScript = "$GRIDSS_ROOT/gridss --jar $GRIDSS_ROOT/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    Int threads = 8
    Int memory = 50
    Int timeout = 100
  }

  command {
  set -euo pipefail

  ~{gridssScript} \
    -s setupreference 
  }
  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }  

}

task preprocess {
  
  input {
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    String gridssScript = "$GRIDSS_ROOT/gridss --jar $GRIDSS_ROOT/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    Int threads = 8
    Int memory = 50
    Int timeout = 100
    File bamFile
  }

  command {
  set -euo pipefail

  ~{gridssScript} \
    -s preprocess ~{bamFile}
  }
  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }  

  output {

    File bamFile = ${bamFile}
  }

}

task assemble {
  
  input {
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    String gridssScript = "$GRIDSS_ROOT/gridss --jar $GRIDSS_ROOT/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    Int threads = 8
    Int memory = 50
    Int timeout = 100
    File normBam
    File tumorBam
  }

  command {
  set -euo pipefail

  ~{gridssScript} \
    -s assemble \
    -a assemble_tumor_norm.bam \
    ~{normBam} ~{tumorBam}
  }
  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }  

  output {

    File bamFile = assemble_tumor_norm.bam
  }

}

task call {
  
  input {
    String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    String gridssScript = "$GRIDSS_ROOT/gridss --jar $GRIDSS_ROOT/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    String refFasta = "$HMFTOOLS_DATA_ROOT/hg38_random.fa"
    Int threads = 8
    Int memory = 50
    Int timeout = 100
    File normBam
    File tumorBam
    File tumorBam
    File tumorBai
    String normName
    String tumorName
    File assembly = assemble_tumor_norm.bam
  }

  command <<<
  set -euo pipefail

  ~{gridssScript} \
    -s call \
    -a ~{assembly} \
    --reference ~{refFasta} \
    --output ~{tumorName}.allocated.vcf \
    ~{normBam} ~{tumorBam}
  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }  

  output {

    File structuralVcf = "~{tumorName}.allocated.vcf"
  }

}

#task call_SVs {
  #input {
    #File normBam
    #File normBai
    #File tumorBam
    #File tumorBai
    #String normName
    #String tumorName
    #String modules = "argparser stringdist structuravariantannotation rtracklayer gridss/2.13.2 hg38/p12 hmftools/1.0 kraken2 bcftools hmftools-data/hg38"
    #String refFasta = "$HMFTOOLS_DATA_ROOT/hg38_random.fa"
    #String gridssScript = "$GRIDSS_ROOT/gridss --jar $GRIDSS_ROOT/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    #Int threads = 8
    #Int memory = 50
    #Int timeout = 100
  #}

  #command <<<
    #set -euo pipefail

    #~{gridssScript} \
      #--reference ~{refFasta} \
      #--output ~{tumorName}.allocated.vcf \
      #~{normBam} ~{tumorBam}

  #>>>

  #runtime {
    #cpu: "~{threads}"
    #memory:  "~{memory} GB"
    #modules: "~{modules}"
    #timeout: "~{timeout}"
  #}

  #output {
    #File structuralVcf = "~{tumorName}.allocated.vcf"
  #}
#}