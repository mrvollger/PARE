"""
Rules for extracting RE sequences from assemblies
"""


rule extract_re_sequences:
    """
    Extract sequences for regulatory elements from their source assembly.
    Creates unslopped bed (exact coordinates) and slopped bed + fasta (with flanking).
    """
    input:
        bed=get_re_bed,
        fasta=get_assembly,
        fai=lambda wc: get_assembly(wc) + ".fai",
    output:
        fasta=temp("temp/sequences/{sample_id}.re_sequences.fa"),
        bed=temp("temp/sequences/{sample_id}.re_sequences.bed"),
        bed_slop=temp("temp/sequences/{sample_id}.re_sequences.slop.bed"),
    log:
        "logs/extract_sequences/{sample_id}.log",
    conda:
        "../envs/bfx.yml"
    threads: 1
    resources:
        mem_mb=get_mem_mb,
        runtime=60,
    params:
        slop=config["filtering"]["slop"],
    shell:
        """
        # Filter bed file by length and create 4-column BED with ID as chrom_start_end
        gunzip -cf {input.bed} | \
        awk -v min={config[filtering][min_re_length]} \
            -v max={config[filtering][max_re_length]} \
            'BEGIN {{OFS="\\t"}} {{
                if($1 !~ /^#/) {{
                    len=$3-$2;
                    if(len>=min && len<=max) {{
                        id=$1"_"$2"_"$3;
                        print $1, $2, $3, id
                    }}
                }}
            }}' \
            > {output.bed} 2> {log}

        # Create slopped version with flanking sequence
        bedtools slop \
            -i {output.bed} \
            -g {input.fai} \
            -b {params.slop} \
            > {output.bed_slop} 2>> {log}

        # Extract slopped sequences for alignment using rb get-fasta (handles bgzip)
        rb get-fasta \
            --fasta {input.fasta} \
            --bed {output.bed_slop} \
            --name \
            > {output.fasta} 2>> {log}
        """


rule merge_assembly_fais:
    """
    Merge all assembly FAI files for chromosome length lookup
    """
    input:
        fais=get_all_assembly_fais(),
    output:
        fai=temp("temp/merged_assemblies.fai"),
    log:
        "logs/merge_fais/all.log",
    threads: 1
    resources:
        mem_mb=1024,
        runtime=5,
    shell:
        """
        cat {input.fais} > {output.fai} 2> {log}
        """


rule merge_unslopped_beds:
    """
    Merge all unslopped (exact) RE BED files for liftover trimming.
    Format for rb liftover --qbed: chrom, start, end, name
    where chrom matches the adjusted PAF query name (chromosome).
    """
    input:
        beds=expand("temp/sequences/{sample_id}.re_sequences.bed", sample_id=get_all_sample_ids()),
    output:
        bed=temp("temp/merged_re_sequences.bed"),
    log:
        "logs/merge_beds/all.log",
    threads: 1
    resources:
        mem_mb=2048,
        runtime=10,
    shell:
        """
        # Just concatenate the unslopped beds - they already have the right format
        # chrom, start, end, name where chrom is the chromosome name
        cat {input.beds} > {output.bed} 2> {log}
        """
