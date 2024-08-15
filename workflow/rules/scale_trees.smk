rule red_scale_tree:
    input:
        "output/gtdb/{domain}_r220.tree"
    output:
        "output/phylorank/{domain}_r220_red_scaled.tree"
    conda:
        "../envs/phylorank.yaml"
    shell:
        "phylorank scale_tree {input} {output}"
