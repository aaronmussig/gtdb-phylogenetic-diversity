rule calculate_pd:
    input:
        gtdb_tree="output/gtdb/{domain}_r220.tree",
        taxonomy="output/gtdb/{domain}_taxonomy_r220.tsv",
        scaled_tree="output/phylorank/{domain}_r220_red_scaled.tree"
    output:
        "output/phylogenetic_diversity_{domain}.tsv"
    conda:
        "../envs/pd.yaml"
    shell:
        "python workflow/scripts/calculate_phylogenetic_diversity.py {input.taxonomy} {input.gtdb_tree} {input.scaled_tree} {output}"
