configfile: "data/configs/pdb_config.yml"

outdir=config["outdir"]
ids=config["ids"]

# SAMPLES WITH PHYID
# read phyids and get 000-code
with open(ids) as input_ids:
    phyids = ["{:04d}".format(int(code)) for code in input_ids]


rule all:
    input:
        expand(outdir+"data/info/info_{code}.tsv", code=phyids),
        expand(outdir+"data/aln/aln_{code}.tar.gz", code=phyids),
        expand(outdir+"data/txts/best_trees_{code}.txt", code=phyids)

 

# maybe user can just input phyid and then you download from ftp
pdb_ftp="ftp://phylomedb.org/phylomedb/phylomes/"

# alignments??

rule Get_phylome_data:
    params:
        trees=pdb_ftp+"phylome_{code}/best_trees.txt.gz",
        aln=pdb_ftp+"phylome_{code}/all_algs.tar.gz",
        info=pdb_ftp+"phylome_{code}/phylome_info.txt.gz"
    output:
        info=outdir+"data/info/info_{code}.tsv",
        aln=outdir+"data/aln/aln_{code}.tar.gz",
        txt=outdir+"data/txts/best_trees_{code}.txt"
    shell:
        """
wget -O - {params.info} | zcat | \
awk '/TaxaID/,0' | awk 'NR>2' | sed 's/ \+\t/\t/g' |  \
taxonkit reformat -I 1 -r "Unassigned"  -f "{{k}}\t{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}" > {output.info}

wget {params.aln} -O {output.aln}

wget {params.trees} -O {output.txt}.gz
gzip -d {output.txt}
"""
