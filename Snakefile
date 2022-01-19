configfile: "config.yaml"
SAMPLES = config['SAMPLES']

rule all: 
   input:
       expand("{sample}_Aligned.out.sam", sample = SAMPLES),
       expand("{sample}.RG.sam", sample = SAMPLES),
       expand("{sample}.dedupped.bam", sample = SAMPLES),
       expand("{sample}.split.bam", sample = SAMPLES),
       expand("{sample}.recal_data.table", sample = SAMPLES),
       expand("{sample}.recalibrated.bam", sample = SAMPLES),
       expand("{sample}.g.vcf", sample=SAMPLES), 
       directory(expand("{my_db}", my_db = config['GENOMEDB'])),
       expand("{cohort}.vcf.gz", cohort=config['ALL_VCF'])

if config['PAIRED']:
    rule trim:
       input:
           r1 = "{sample}.r_1.fq.gz",
           r2 = "{sample}.r_2.fq.gz"
       output:
           "galore/{sample}.r_1_val_1.fq.gz",
           "galore/{sample}.r_2_val_2.fq.gz"
       conda: 'env/env-trim.yaml'
       shell:
           """
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
           """
    rule align:
        input:
           "galore/{sample}.r_1_val_1.fq.gz",
           "galore/{sample}.r_2_val_2.fq.gz"
        output:
             "{sample}_Aligned.out.sam"
        params:
             threads = config['THREADS'],
             gtf = config['GTF'],
             prefix = "{sample}_",
             index = config['INDEX']
        shell:
           """
           STAR --genomeDir {params.index} --runThreadN {params.threads} --readFilesIn {input[0]} {input[1]}  --outFileNamePrefix {params.prefix} --sjdbGTFfile {params.gtf}  --twopassMode Basic
           """


else:
     rule trim:
       input:
           "{sample}.fq.gz",

       output:
           "galore/{sample}_trimmed.fq.gz",
       conda: 'env/env-trim.yaml'
       shell:
           """
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore {input}
           """

     rule align:
        input:
              "galore/{sample}_trimmed.fq.gz"
        output:
             "{sample}_Aligned.out.sam"
        params:
             threads = config['THREADS'],
             gtf = config['GTF'],
             prefix = "{sample}_",
             index = config['INDEX']
        shell:
           """
           STAR --genomeDir {params.index} --runThreadN {params.threads} --readFilesIn {input}  --outFileNamePrefix {params.prefix} --sjdbGTFfile {params.gtf}  --twopassMode Basic
           """

rule AddRG:
    input:
       "{sample}_Aligned.out.sam"
    output:
       "{sample}.RG.sam"
    params:
       RG = config['RG']
    conda: "env/env-picard.yaml"
    shell:
        """
        picard AddOrReplaceReadGroups I={input} O={output} SO=coordinate RGID=@{params} RGSM={wildcards.sample} RGPL=Illumina RGLB={wildcards.sample} RGPU={params}_{wildcards.sample} VALIDATION_STRINGENCY=SILENT
        """

rule deduplicate: 
    input: 
        "{sample}.RG.sam" 
    output: 
        "{sample}.dedupped.bam", 
        "{sample}.output.metrics"
    shell: 
        """
        picard MarkDuplicates I={input} O={output[0]}  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={output[1]} 
        """

rule split: 
    input: 
       "{sample}.dedupped.bam" 
    params: 
       genome= config['GENOME']
    output: 
       "{sample}.split.bam" 
    shell: 
       """
       gatk SplitNCigarReads -R {params.genome} -I {input} -O {output} 
       """

rule recalibrate_a: 
    input:
       "{sample}.split.bam"
    params:
        genome= config['GENOME'],
        DBSNP = config['DBSNP'],
        INDELS = config['INDELS'],
        GOLD_STANDARD = config['GOLD_STANDARD'] 
    output:
       "{sample}.recal_data.table"
    shell: 
       """
       gatk BaseRecalibrator -I {input} -R {params.genome}  --known-sites {params.DBSNP} --known-sites {params.INDELS} --known-sites {params.GOLD_STANDARD} -O {output} 
       """

rule recalibrate_b: 
    input: 
       "{sample}.recal_data.table"
    params: 
       genome= config['GENOME'], 
    output: 
       "{sample}.recalibrated.bam"
    shell: 
      """ 
        gatk ApplyBQSR -I {input} -R {params.genome} --bqsr-recal-file {input} --disable-sequence-dictionary-validation -O {output}
      """

rule tovcf:
   input:
      "{sample}.recalibrated.bam",
      expand("{genome}.fasta", genome=config['GENOME'])
   params:
     mem_threads = {"-Xmx100g -XX:ParallelGCThreads=4"}
   log: "logs/{sample}.tovcf.log"
   benchmark: "logs/{sample}.tovcf.benchmark"
   conda: 'env/env-gatk.yaml'
   output:
       "{sample}.g.vcf"
   shell:
       """
       gatk --java-options "{params.mem_threads}" HaplotypeCaller -R {input[1]} -I {input[0]} -ERC GVCF -O {output[0]} -G StandardAnnotation -G AS_StandardAnnotation
       """

rule GenomeDBImport:
     input:
         sample = expand("{sample}.g.vcf", sample = SAMPLES),
         genome = expand("{genome}.fasta", genome=config['GENOME'])
     params:
         INTERVALS = config['INTERVALS_FILE'],
         DB = config['GENOMEDB'],
         mem = {"-Xmx100g"},
         I =  lambda w: "-V " + " -V ".join(expand("{sample}.g.vcf", sample = SAMPLES))
     log: "logs/dbimport.log"
     benchmark: "logs/DBImport.benchmark"
     conda: 'env/env-gatk.yaml'
     output:
        directory(expand("{my_db}", my_db = config['GENOMEDB']))
     shell:
         """
         gatk --java-options {params.mem} GenomicsDBImport {params.I} --genomicsdb-workspace-path {params.DB} -L {params.INTERVALS}
         """

rule jointcall:
    input:
       DB = directory(expand("{my_db}", my_db = config['GENOMEDB'])),
       genome = expand("{genome}.fasta", genome=config['GENOME'])
    params:
       mem = {"-Xmx100g"}
    log: "logs/jointcall.log"
    benchmark: "logs/jointcall.benchmark"
    conda: 'env/env-gatk.yaml'
    output:
       expand("{cohort}.vcf.gz", cohort=config['ALL_VCF'])
    shell:
        """
          gatk --java-options {params.mem} GenotypeGVCFs -R {input.genome} -V gendb://{input.DB} -O {output} -G StandardAnnotation -G AS_StandardAnnotation
        """


