from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

rule gtdb_tree:
    input:
        HTTP.remote("https://data.gtdb.ecogenomic.org/releases/release220/220.0/{domain}_r220.tree.gz",keep_local=False)
    output:
        "output/gtdb/{domain}_r220.tree"
    shell:
        "gunzip -c {input} > {output}"


rule gtdb_taxonomy:
    input:
        HTTP.remote("https://data.gtdb.ecogenomic.org/releases/release220/220.0/{domain}_taxonomy_r220.tsv.gz",keep_local=False)
    output:
        "output/gtdb/{domain}_taxonomy_r220.tsv"
    shell:
        "gunzip -c {input} > {output}"
