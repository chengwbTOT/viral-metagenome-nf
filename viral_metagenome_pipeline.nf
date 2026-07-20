// Enable DSL2 standard syntax
nextflow.enable.dsl = 2

/* ======================================
   Configurable parameters (match all settings from the original script, adjust as needed)
   ====================================== */
params.threads = 40                     // Number of CPU threads
params.assembly_memory = 500.GB         // Memory allocated for assembly step
params.min_contig_length = 5000         // Minimum contig length threshold for ViWrap
params.pfam_evalue = "1e-5"             // E-value cutoff for Pfam annotation

// Database paths (identical to original script paths, modifiable as needed)
params.viwrap_db = "/opt/user/chengwb/Database/ViWrap_db"
params.pfam_db = "/opt/user/chengwb/Database/Pfam_v35.0/v37.1/Pfam-A.hmm"

// Conda environment paths (correspond to 3 environment switches in the original script)
params.conda_metawrap = "metawrap-env"
params.conda_viwrap = "/opt/user/chengwb/Database/yml_enviroments/ViWrap"
params.conda_hmmer = "/opt/user/chengwb/Database/yml_enviroments/hmmer"

// Root output directory
params.outdir = "results"

/* ======================================
   Process 1: Read quality control (corresponds to metawrap read_qc step in original script)
   ====================================== */
process READ_QC {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/02_READ_QC", mode: 'copy'
    conda params.conda_metawrap  // Auto-activate metawrap environment, replaces source activate

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("final_pure_reads_1.fastq"), path("final_pure_reads_2.fastq")

    script:
    """
    metawrap read_qc -1 ${r1} \
                     -2 ${r2} \
                     -t ${params.threads} \
                     -o ./ \
                     --skip-bmtagger
    """
}

/* ======================================
   Process 2: Metagenomic assembly (corresponds to metaspades step in original script)
   ====================================== */
process METASPADES_ASSEMBLY {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/03_ASSEMBLY", mode: 'copy'
    conda params.conda_metawrap
    memory params.assembly_memory

    input:
    tuple val(sample_id), path(r1_clean), path(r2_clean)

    output:
    tuple val(sample_id), path("contigs.fasta")

    script:
    """
    spades.py --meta -1 ${r1_clean} \
                  -2 ${r2_clean} \
                  -m ${params.assembly_memory.toGiga()} \
                  -t ${params.threads} \
                  -o ./
    """
}

/* ======================================
   Process 3: Viral identification and annotation via ViWrap (corresponds to ViWrap run step in original script)
   ====================================== */
process VIWRAP_ANNOTATION {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/04_ViWrap", mode: 'copy'
    conda params.conda_viwrap

    input:
    tuple val(sample_id), path(contigs), path(r1_clean), path(r2_clean)

    output:
    tuple val(sample_id), path("08_ViWrap_summary_outdir")

    script:
    """
    ViWrap run --input_metagenome ${contigs} \
               --input_reads ${r1_clean},${r2_clean} \
               --out_dir ./ \
               --db_dir ${params.viwrap_db} \
               --identify_method vs \
               --conda_env_dir /opt/user/chengwb/Database/yml_enviroments \
               --threads ${params.threads} \
               --input_length_limit ${params.min_contig_length}
    """
}

/* ======================================
   Process 4: Result curation and merging (corresponds to data processing step in original script)
   ====================================== */
process RESULT_SUMMARY {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/04_ViWrap/08_ViWrap_summary_outdir", mode: 'copy'

    input:
    tuple val(sample_id), path(summary_dir)

    output:
    tuple val(sample_id), path("Virus_summary_simplify.txt"), path("Virus_genomes_files")

    script:
    """
    cd ${summary_dir}

    # Extract high/medium quality viral genomes
    grep -E "Complete|High-quality|Medium-quality" Virus_summary_info.txt > Virus_summary_tmp1.txt

    # Merge normalized abundance information
    gawk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[\$1]=\$5}ARGIND==2{print \$0,a[\$1]}' Virus_normalized_abundance.txt Virus_summary_tmp1.txt > Virus_summary_tmp2.txt

    # Merge taxonomic classification information
    gawk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[\$1]=\$2}ARGIND==2{print \$0,a[\$1]}' Tax_classification_result.txt Virus_summary_tmp2.txt > Virus_summary_tmp3.txt

    # Merge host prediction information
    sed 's/,/\t/g' Host_prediction_to_genus_m90.csv > Host_prediction_to_genus_m90.txt
    gawk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[\$1]=\$3;b[\$1]=\$4}ARGIND==2{print \$0,a[\$1],b[\$1]}' Host_prediction_to_genus_m90.txt Virus_summary_tmp3.txt > Virus_summary_tmp4.txt

    # Sort by abundance
    sort -nr -k 11 -t \$'\t' Virus_summary_tmp4.txt > Virus_summary_simplify.txt

    # Add column headers
    sed -i '1 i Virus\tgenome_size\tscaffold_num\tprotein_count\tAMG_KOs\tlytic_state\tcheckv_quality\tmiuvig_quality\tcompleteness\tcompleteness_method\tMeanCov.Percent\tTax_classification\tHost genus\tHost genus Confidence score' Virus_summary_simplify.txt

    # Rename viral sequences (add sample prefix)
    rename "s/vRhyme/${sample_id}_vRhyme/g" ./Virus_genomes_files/*
    sed -i "s/^>/>${sample_id}/g" ./Virus_genomes_files/*

    # Copy results back to working directory
    cp Virus_summary_simplify.txt ${task.workDir}/
    cp -r Virus_genomes_files ${task.workDir}/
    """
}

/* ======================================
   Process 5: Pfam functional annotation (corresponds to hmmscan step in original script)
   ====================================== */
process PFAM_ANNOTATION {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/04_ViWrap/08_ViWrap_summary_outdir/Annotation/Pfamv37.1", mode: 'copy'
    conda params.conda_hmmer

    input:
    tuple val(sample_id), path(virus_genomes_dir)

    output:
    path("*.hmm.tbl")

    script:
    """
    mkdir -p ./
    for faa_file in ${virus_genomes_dir}/*.faa
    do
        name=\$(basename \$faa_file .faa)
        hmmscan -o ./\${name}.hmm.out \
                --tblout ./\${name}.hmm.tbl \
                --noali -E ${params.pfam_evalue} --cpu ${params.threads} \
                ${params.pfam_db} \
                \$faa_file
    done
    """
}

/* ======================================
   Main workflow: connect all analysis steps
   ====================================== */
workflow {
    // Read all paired-end sequencing files (supports batch samples, auto-matches *_R{1,2}.fastq format)
    read_pairs = Channel.fromFilePairs("raw_data/*_R{1,2}.fastq")

    // Pipeline chaining
    step_qc = READ_QC(read_pairs)
    step_assembly = METASPADES_ASSEMBLY(step_qc)
    step_viwrap = VIWRAP_ANNOTATION(step_assembly.join(step_qc))
    step_summary = RESULT_SUMMARY(step_viwrap)
    step_pfam = PFAM_ANNOTATION(step_summary)

    // Completion notification
    step_pfam.view { "Sample ${it[0]} analysis completed. Summary output: ${params.outdir}/${it[0]}/04_ViWrap/08_ViWrap_summary_outdir/Virus_summary_simplify.txt" }
}
