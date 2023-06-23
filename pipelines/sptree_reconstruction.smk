configfile: "data/configs/pdb_config.yml"

# DEPENDENCIES
# newick_utils
# ete3 
# treeswift
# astral
# duptree
# astral pro


# to run activate snakemake 

# if error cran repo
# options(repos = c(CRAN = "http://cran.rstudio.com"))

# if weird error pthreads 
# export OMP_NUM_THREADS=1
# export USE_SIMPLE_THREADED_LEVEL3=1

# for MN4 do 
# conda activate snakemake

outdir=config["outdir"]
disco=config["disco"]
minimum_occupancy=config['minimum_occupancy']

# Species 2 age dictionary
with open(config["s2a"]) as s2a_file:
    s2a = eval(s2a_file.read())

# SAMPLES WITH PHYID
# read phyids and get 000-code
with open(config['ids']) as input_ids:
    phyids = ["{:04d}".format(int(code)) for code in input_ids]

# species tree reconstruction method
methods = ["ad", "apro", "asteroid", "astral", "dt", "sprax", "wastral", "fmulrfs"] # "cat"

rule results:
    input:
        expand(outdir+"plots/phy_{code}_all_presence.png", code=phyids),
        # THIS MAY BE USELESS ALTHOUGH VERY COOL! some times out of mem error
        # expand(outdir+"plots/phy_{code}_redundancy.pdf", code=phyids),
        # expand(outdir+"data/all_trees/rooted_mv_best_trees_{code}.nwk", code=phyids),
        # expand(outdir+"sptrees/phylome_{code}_sptree_dt.nwk", code=phyids),
        # expand(outdir+"sptrees/phylome_{code}_sptree_apro.nwk", code=phyids),
        # expand(outdir+"sptrees/phylome_{code}_sptree_ad.nwk", code=phyids),
        # expand(outdir+"sptrees/phylome_{code}_sptree_astral.nwk", code=phyids),
        # expand(outdir+"sptrees/phylome_{code}_sptree_wastral.nwk", code=phyids),
        # expand(outdir+"sptrees/phylome_{code}_sptree_asteroid.nwk", code=phyids),
        # expand(outdir+"sptrees/phylome_{code}_sptree_cat.nwk", code=phyids),
        # expand(outdir+"sptrees/phylome_{code}_sptree_sprax.nwk", code=phyids),
        # expand(outdir+"sptrees/phylome_{code}_sptree_fmulrfs.nwk", code=phyids),
        expand(outdir+"data/txts/treestats_{code}.tsv", code=phyids),
        expand(outdir+"sptrees/phylome_{code}_sptree.nwk", code=phyids)


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
        stats=outdir+"data/txts/treestats_{code}.tsv",
        min_oc=minimum_occupancy
    shell:
        """
awk -v min={params.min_oc} '$12>=min && $13=="True"' {params.stats} | cut -f1 > {output.ids}
awk '{{print "all_algs/"$1".clean.fasta"}}' {output.ids} > {output.alns_paths}

join -1 1 -2 1 <(sort {input.trees}) <(sort {output.ids}) | cut -f4 -d' ' > {output.sc}
sed -E 's/Phy[A-Z0-9]{{7}}_//g' {output.sc} > {output.sc_l}

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

# SPECIES TREE

## SINGLE COPY

### CONCAT

rule Concat:
    input:
        outdir+"data/aln/aln_{code}_sc"
    output:
        st=outdir+"sptrees/phylome_{code}_sptree_cat.nwk"
    log:
        outdir+"log/phy_{code}_cat.log"
    params:
        prefix=outdir+"input/cat/phylome_{code}_cat"
    threads:
        8
    shell:
        """
mkdir -p $(dirname {params.prefix})

if [ -f {input}/.snakemake_timestamp ]; then
    rm {input}/.snakemake_timestamp
fi

iqtree2 -T {threads} -s {input} --prefix {params.prefix} -bb 1000 -m LG+G4+I --undo > {log}

cp {params.prefix}.contree {output.st}
"""
# -mset DCmut,JTTDCMut,LG,WAG,VT
# with -s it assumes a single model, with -p it will infer a concatenated partitioned alignment, this can get very slow with
# > 100 genes so I'll use -s for now


### ASTRAL

rule Astral:
    input:
        outdir+"data/single_copy/best_trees_{code}_sc_labs.nwk"
    output:
        gt=outdir+"input/astral/phylome_{code}_input_astral.nwk",
        st=outdir+"sptrees/phylome_{code}_sptree_astral.nwk"
    log:
        outdir+"log/phy_{code}_astral.log"
    threads:
        8
    params:
        collapse=config['collapse']
    shell:
        """
mkdir -p $(dirname {log})

nw_ed {input} 'i & b<={params.collapse}' o > {output.gt}
astral -t {threads} -i {output.gt} -o {output.st}  2>&1 | tee {log}
"""
# echo "collapsing nodes with support < than {params.collapse}"
# this is wASTRAL unweighted from ASTER suite 

# rule astral:
#     input:
#         outdir+"data/single_copy/best_trees_{code}_sc_labs.nwk"
#     output:
#         gt=outdir+"input/astral/phylome_{code}_input_astral.nwk",
#         st=outdir+"sptrees/phylome_{code}_sptree_astral.nwk"
#     log:
#         outdir+"log/phy_{code}_astral.log"
#     params:
#         astral=astral_java,
#         collapse=config['collapse']
#     shell:
#         """
# mkdir -p $(dirname {log})

# echo "collapsing nodes with support < than {params.collapse}"
# nw_ed {input} 'i & b<={params.collapse}' o > {output.gt}
# java -jar {params.astral} -i {output.gt} -o {output.st}  2>&1 | tee {log}
# """
# this is astralIII 5.7.8
# -t 32 for discovista stuff

rule wASTRAL:
    input:
        outdir+"data/single_copy/best_trees_{code}_sc_labs.nwk"
    output:
        # gt=outdir+"input/astral/phylome_{code}_input_wastral.nwk",
        outdir+"sptrees/phylome_{code}_sptree_wastral.nwk"
    log:
        outdir+"log/phy_{code}_wastral.log"
    # params:
    #     collapse=config['collapse']
    threads:
        8
    shell:
        """
mkdir -p $(dirname {log})
astral-hybrid -x 1 -t {threads} -L -o {output} {input} 2> {log}
"""


## MULTI COPY STUFF


### ASTEROID

rule Asteroid:
    input:
        mc=outdir+"data/all_trees/best_trees_{code}.nwk",
        sc=outdir+"data/single_copy/best_trees_{code}_sc_labs.nwk"
    output:
        mapping=outdir+"input/ad/phylome_{code}_mapping_ad.txt",
        st=outdir+"sptrees/phylome_{code}_sptree_asteroid.nwk"
    params:
        prefix=outdir+"sptrees/phylome_{code}_asteroid_tmp"
    threads:
        8
    log:
        outdir+"log/phy_{code}_asteroid.log"
    shell:
        """
set +e
python scripts/prepare_trees.py -i {input.mc} -o {output.mapping} -m asteroid
mpiexec --oversubscribe -np {threads} asteroid -i {input.mc} -b 100 -p {params.prefix} -m {output.mapping} 2>&1 | tee {log}
exitcode=$?

if [ $exitcode -eq 134 ]; then
    echo "Asteroid failed with multi copy gene trees! The resulting species tree will be based on single copy gene trees" >> {log}
    mpiexec --oversubscribe -np {threads} asteroid -i {input.sc} -b 100 -p {params.prefix} 2>&1 | tee -a {log}
fi

mv {params.prefix}.bestTree.newick {output.st}
rm {params.prefix}.allTrees.newick
rm {params.prefix}.scores.txt
echo "bs trees" >> {log}
cat {params.prefix}.bsTrees.newick >> {log}
rm {params.prefix}.bsTrees.newick
"""

# bool AsteroidOptimizer::computeAndApplyBestSPR(): Assertion `newScore > _lastScore' failed.
# /bin/bash: line 2: 279212 Aborted
# mpiexec -np {threads} this does segfault even with one thread
# --use-gene-bl ??????? I think it's smarter comparde to topological distance 
# "In my experience, using gene branch lengths results (on average) in less accurate trees, so I don't necessarily recommend using this option."
# you can parallelize with mpiexec
# ASTRID does not make sense if using asteroid as its worse in all cases


### ASTRAL-PRO

rule Astral_PRO:
    input:
        outdir+"data/all_trees/best_trees_{code}_labs.nwk"
    output:
        st=outdir+"sptrees/phylome_{code}_sptree_apro.nwk"
    log:
        log=outdir+"log/phy_{code}_apro.log"
    threads:
        8
    shell:
        """
mkdir -p $(dirname {log})
astral-pro -u 1 -t {threads} -o {output.st} {input} 2>&1 | tee {log}
"""

### DISCO

rule DISCO:
    input:
        outdir+"data/all_trees/best_trees_{code}.nwk"
    output:
        outdir+"data/all_trees/best_trees_{code}_disco.nwk",
    shell:
        "python scripts/disco.py -i {input} -o {output} -d \"_\" -n 1"


### ASTRAL-DISCO

rule Astral_DISCO:
    input:
        outdir+"data/all_trees/best_trees_{code}_disco.nwk",
    output:
        outdir+"sptrees/phylome_{code}_sptree_ad.nwk"
    log:
        outdir+"log/phy_{code}_ad.log"
    params:
        # astral=astral_java
        collapse=config['collapse']
    threads:
        8
    shell:
        """
mkdir -p $(dirname {log})
nw_ed {input} 'i & b<={params.collapse}' o | \
astral -u 1 -t {threads} -i /dev/stdin -o {output}  2>&1 | tee {log}
"""


### SpeciesRax

rule SpeciesRax:
    input:
        outdir+"data/all_trees/best_trees_{code}.nwk"
    output:
        gt=directory(outdir+"input/sprax/phylome_{code}"),
        gt2=directory(outdir+"input/sprax/phylome_{code}/results"),
        st=outdir+"sptrees/phylome_{code}_sptree_sprax.nwk"
    params:
        outdir+"data/txts/best_trees_{code}.txt"
    log:
        outdir+"log/phy_{code}_sprax.log"
    threads:
        8
    shell:
        """
python scripts/prepare_trees.py -i {params} -o {output.gt} -m speciesrax

mpiexec --oversubscribe -np {threads} generax --families {output.gt}/family.txt \
--strategy SKIP --si-strategy HYBRID --species-tree MiniNJ \
--rec-model UndatedDTL --per-family-rates --prune-species-tree \
--si-estimate-bl --si-quartet-support  --prefix {output.gt2} > {log}

cp {output.gt2}/species_trees/inferred_species_tree.newick {output.st}
"""

# an optional starting species tree
# a set of unrooted gene trees OR a set of multiple sequence alignments.
# the mapping between gene taxa and species taxa.

### duptree
rule Duptree:
    input:
        outdir+"data/all_trees/best_trees_{code}_labs.nwk"
    output:
        gt=outdir+"input/dt/phylome_{code}_input_dt.nwk",
        st=outdir+"sptrees/phylome_{code}_sptree_dt.nwk"
    log:
        outdir+"log/phy_{code}_dt.log"
    shell:
        """
mkdir -p $(dirname {log})

python scripts/prepare_trees.py -i {input} -o {output.gt} -m duptree
duptree -i {output.gt} -o {output.st} --nogenetree > {log}
"""
# echo "this is a test, only using 100 trees, remember to remove this"

### fastmulrfs

rule FastMulRFS:
    input:
        outdir+"data/all_trees/best_trees_{code}_labs.nwk"
    output:
        gt=outdir+"input/fmulrfs/phylome_{code}_input_fmulrfs.nwk",
        st=outdir+"sptrees/phylome_{code}_sptree_fmulrfs.nwk"
    log:
        outdir+"log/phy_{code}_fmulrfs.log"
    shell:
        """
prefix=$(echo {output.gt} | sed 's/_input_fmulrfs.nwk/_results_fmulrfs/')
./scripts/run_fastmulrfs.sh {input} {output.gt} $prefix &> {log}
cp ${{prefix}}.single {output.st}
"""

# QC


### Final sptrees

rule Get_consensus:
    input:
        sptrees=expand(outdir+"sptrees/phylome_{{code}}_sptree_{method}.nwk", method=methods),
        alns=outdir+"data/aln/aln_{code}_sc",
        sc_l=outdir+"data/single_copy/best_trees_{code}_sc_labs.nwk"
    output:
        contree=temp(outdir+"sptrees/phylome_{code}_contree.nwk"),
        st=outdir+"sptrees/phylome_{code}_sptree.nwk"
    params:
        con_prefix=outdir+"sptrees/phylome_{code}_consensus",
    log:
        outdir+"log/phy_{code}_consensus.log"
    threads:
        8
    shell:
        """
cat {input.sptrees} | nw_reroot -d - | nw_topology -I - > {output.contree}

iqtree2 -con -t {output.contree}
rm {output.contree}.log
mv {output.contree}.contree {output.contree}

if [ -f {input.alns}/.snakemake_timestamp ]; then
    rm {input.alns}/.snakemake_timestamp
fi

iqtree2 -te {output.contree} --gcf {input.sc_l} -p {input.alns} --scf 100 --prefix {params.con_prefix} -T {threads}

mv {params.con_prefix}.cf.tree {output.st}
mv {params.con_prefix}.log {log}
rm {params.con_prefix}.cf.tree.nex {params.con_prefix}.cf.branch {params.con_prefix}.cf.stat
"""
# iqtree2 -p output/data/aln/aln_0096_sc/ -t output/results/phylome_0096_sptree.nwk
# this to infer branch lengths

rule Plot_sptrees:
    input:
        sts=expand(outdir+"sptrees/phylome_{{code}}_sptree_{method}.nwk", method=methods),
        consensus=outdir+"sptrees/phylome_{code}_sptree.nwk"
    output:
        # directory(outdir+"results/phylome_{code}")
        outdir+"results/phylome_{code}/phylome_{code}_rooted_sptree.nwk"
    shell:
        """
Rscript scripts/compare_sptrees.R -i $(echo {input.sts} | sed -e 's/\s\+/,/g') -c {wildcards.code} -o {output} -s {input.consensus}
"""

rule Plot_occupancy:
    input:
        trees=outdir+"data/all_trees/best_trees_{code}.nwk",
        sc=outdir+"data/single_copy/best_trees_{code}_sc.nwk",
        info=outdir+"data/info/info_{code}.tsv",
        consensus=outdir+"results/phylome_{code}/phylome_{code}_rooted_sptree.nwk"
    output:
        all_t=outdir+"plots/phy_{code}_all_presence.png",
        sc=outdir+"plots/phy_{code}_sc_presence.png",
        occ=outdir+"plots/phy_{code}_all_oc.png"
    params:
        mo=minimum_occupancy
    shell:
        """
mkdir -p $(dirname {output.all_t})

Rscript scripts/plot_hm.R -g {input.trees} -o {output.all_t} -i {input.info} --occupancy {output.occ} -s {input.consensus}
Rscript scripts/plot_hm.R -a -g {input.sc} -o {output.sc} -i {input.info} -s {input.consensus} -m {params.mo}
"""



# Get support and branch lengths for consensus species tree

# rule Root_trees:
#     input:
#         outdir+"data/all_trees/best_trees_{code}.nwk"
#     output:
#         mv=outdir+"data/all_trees/rooted_mv_best_trees_{code}.nwk",
#         mid=outdir+"data/all_trees/rooted_mid_best_trees_{code}.nwk",
#         mad=outdir+"data/all_trees/rooted_mad_best_trees_{code}.nwk"
#     log:
#         mv=outdir+"log/phy_{code}_mv.log",
#         mid=outdir+"log/phy_{code}_mid.log"
#         mad=outdir+"log/phy_{code}_mad.log"
#     shell:
#         """
# FastRoot.py -i {input} -o {output.mv} -f {log.mv}
# FastRoot.py -i {input} -o {output.mid} -f {log.mid} -m MP
# ./scripts/mad {input} -n > {log.mad}
# mv {input}.rooted {output.mad}
# """

# https://github.com/larabreithaupt/seastaR covariance matrix between species

# SPTREE
# CA_disco may also be cool although super comp demanding

# GT ROOTING
# STRIDE does not work properly and would need significant recoding, not very different from duptree rooting also
# QR rooting, only single copy???
# which tree to root???

