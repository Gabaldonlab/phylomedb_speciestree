configfile: "data/configs/pdb_config.yml"

outdir=config["outdir"]
to_download=config["download_ids"]
ids=config["ids"]
failed=config["failed_ids"]

# SAMPLES WITH PHYID
# read phyids and get 000-code
with open(to_download) as input_ids:
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
        sptree_ids=ids,
        failed=failed,
        trees=pdb_ftp+"phylome_{code}/best_trees.txt.gz",
        aln=pdb_ftp+"phylome_{code}/all_algs.tar.gz",
        info=pdb_ftp+"phylome_{code}/phylome_info.txt.gz",
        phy_id=lambda wcs: wcs.code
    output:
        info=outdir+"data/info/info_{code}.tsv",
        aln=outdir+"data/aln/aln_{code}.tar.gz",
        txt=outdir+"data/txts/best_trees_{code}.txt"
    shell:
        """
wget -q -O - {params.info} | zcat | awk '/TaxaID/,0' | awk 'NR>2' | sed 's/ \+\t/\t/g' |  \
taxonkit reformat -I 1 -r "Unassigned"  -f "{{k}}\t{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}" > {output.info}

wget {params.aln} -O {output.aln}

wget -q {params.trees} -O {output.txt}.gz
gzip -d {output.txt}

n_alns=$(tar tvf {output.aln} | grep -c clean)
n_trees=$(wc -l < {output.txt})

id={params.phy_id}
phyid=$(expr "$id" + 0)

if (( $n_alns < $n_trees )); then

    echo "less alignment ($n_alns) than trees ($n_trees)"
    echo "check phylome data $phyid manually!"
    echo $phyid >> {params.failed}
    sort -u {params.failed} -o {params.failed}
    echo "added $phyid in failed ids"

    rm {output.aln} {output.txt} {output.info}

    exit 1

else

    echo $phyid >> {params.sptree_ids}
    sort -u {params.sptree_ids} -o {params.sptree_ids}
    echo "added $phyid in sptree reconstuction id"

fi
"""
