configfile: "data/configs/pdb_config.yml"

outdir=config["outdir"]
to_download=config["download_ids"]
ids=config["ids"]
failed=config["failed_ids"]
minimum_occupancy=config['minimum_occupancy']

pdb_ftp="ftp://phylomedb.org/phylomedb/phylomes/"

# SAMPLES WITH PHYID
# read phyids and get 000-code
with open(to_download) as input_ids:
    phyids = ["{:04d}".format(int(code)) for code in input_ids]


rule all:
    input:
        expand(outdir+"data/info/info_{code}.tsv", code=phyids),
        # expand(outdir+"data/aln/aln_{code}.tar.gz", code=phyids),
        # expand(outdir+"data/txts/best_trees_{code}.txt", code=phyids),
        expand(outdir+"data/txts/treestats_{code}.tsv", code=phyids),
        expand(outdir+"data/all_trees/best_trees_{code}_labs.nwk", code=phyids),
        expand(outdir+"data/all_trees/best_trees_{code}.nwk", code=phyids),
        expand(outdir+"data/all_trees/best_trees_{code}_labs.nwk", code=phyids),
        expand(outdir+"data/aln/aln_{code}_sc", code=phyids),
        expand(outdir+"data/single_copy/ids_{code}.txt", code=phyids),
        expand(outdir+"data/single_copy/best_trees_{code}_sc.nwk", code=phyids),
        expand(outdir+"data/single_copy/best_trees_{code}_sc_labs.nwk", code=phyids)

 
rule Download_phylome_data:
    params:
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

    echo "FTP files looks alright, proceeding"

fi
"""

rule Get_tree_stats:
    input:
        outdir+"data/txts/best_trees_{code}.txt"
    output:
        outdir+"data/txts/treestats_{code}.tsv"
    shell:
        """
nsp=$(cut -f4 {input} | nw_labels -I - | cut -f2 -d'_' | sort -u | wc -l)
python scripts/get_stats.py -i {input} -o {output} -n $nsp
"""


rule Get_all_trees:
    input:
        outdir+"data/txts/best_trees_{code}.txt"
    output:
        nwk=outdir+"data/all_trees/best_trees_{code}.nwk",
        nwk_l=outdir+"data/all_trees/best_trees_{code}_labs.nwk",
    shell:
        """
cut -f4 {input} > {output.nwk}
sed -E "s/Phy[A-Z0-9]{{7}}_//g" {output.nwk} > {output.nwk_l}
"""


rule Get_single_copy:
    input:
        trees=outdir+"data/txts/best_trees_{code}.txt",
        stats=outdir+"data/txts/treestats_{code}.tsv",
        alns=outdir+"data/aln/aln_{code}.tar.gz"
    output:
        alns=directory(outdir+"data/aln/aln_{code}_sc"),
        ids=outdir+"data/single_copy/ids_{code}.txt",
        alns_paths=temp(outdir+"data/single_copy/ids_toextract_{code}.txt"),
        sc=outdir+"data/single_copy/best_trees_{code}_sc.nwk",
        sc_l=outdir+"data/single_copy/best_trees_{code}_sc_labs.nwk"
    log:
        outdir+"log/phy_{code}_sc.log",
    params:
        sptree_ids=ids,
        failed=failed,
        phy_id=lambda wcs: wcs.code,
        min_oc=minimum_occupancy
    shell:
        """
awk -v min={params.min_oc} '$12>=min && $13=="True"' {input.stats} | cut -f1 > {output.ids}
awk '{{print "all_algs/"$1".clean.fasta"}}' {output.ids} > {output.alns_paths}


join -1 1 -2 1 <(sort {input.trees}) <(sort {output.ids}) | cut -f4 -d' ' > {output.sc}
sed -E 's/Phy[A-Z0-9]{{7}}_//g' {output.sc} > {output.sc_l}

n_allsp=$(cut -f4 {input.trees} | nw_labels -I - | cut -f2 -d'_' | sort -u | wc -l)
n_scsp=$(cat {output.sc_l} | nw_labels -I - | cut -f2 -d'_' | sort -u | wc -l)

id={params.phy_id}
phyid=$(expr "$id" + 0)

if (( $n_scsp < $n_allsp )); then

    echo "less species in single copy trees ($n_scsp) than all species ($n_allsp)"
    echo $phyid >> {params.failed}
    sort -u {params.failed} -o {params.failed}
    echo "added $phyid in failed ids"

    exit 1

else

    echo $phyid >> {params.sptree_ids}
    sort -u {params.sptree_ids} -o {params.sptree_ids}
    echo "added $phyid in sptree reconstuction id"

fi

echo "$(wc -l < {output.sc}) out of $(wc -l < {input.trees}) are single copy" > {log}

echo "Extracting alignments for single copy with occupancy greater than {params.min_oc}." >> {log}
echo "Gene names will be substituted with species names." >> {log}

mkdir -p {output.alns}

tar -C {output.alns} --strip-components 1 -xzf {input.alns} --files-from {output.alns_paths}

python scripts/check_alns.py -i {output.alns}

sed -i -E 's/^>Phy[A-Z0-9]{{7}}_/>/g' {output.alns}/*fasta

echo "Done!" >> {log}

if [ -f {output.alns}/.snakemake_timestamp ]; then
    rm {output.alns}/.snakemake_timestamp
fi
"""
# ./scripts/get_121.sh {input.trees} {output.sc}
# num_sp=$(nw_labels {output.sc_l} -I | sort | uniq | wc -l)
# min_sp=$(echo "{params.min_oc}*$num_sp/1" | bc)
#  nw_topology -I output/data/single_copy/best_trees_0003_sc_labs.nwk | nw_order -c a - | sort | uniq -c | sort -n -k1 | awk '{$1=$1};1' | sed 's/ /\t/' | awk '$1>1'
