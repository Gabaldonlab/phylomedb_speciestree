configfile: "data/configs/pdb_config.yml"


# DEPENDENCIES
# newick_utils
# ete3 
# treeswift
# astral
# duptree
# astral pro

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
        expand(outdir+"results/phylome_{code}/phylome_{code}_all_presence.pdf", code=phyids),
        # THIS MAY BE USELESS ALTHOUGH VERY COOL! some times out of mem error
        # expand(outdir+"plots/phy_{code}_redundancy.pdf", code=phyids),
        expand(outdir+"results/phylome_{code}/phylome_{code}_rooted_sptree.nwk", code=phyids)

# SPECIES TREE

## SINGLE COPY

### CONCAT

# rule Concat:
#     input:
#         outdir+"data/aln/aln_{code}_sc"
#     output:
#         st=outdir+"sptrees/phylome_{code}_sptree_cat.nwk"
#     log:
#         outdir+"log/phy_{code}_cat.log"
#     params:
#         prefix=outdir+"input/cat/phylome_{code}_cat"
#     threads:
#         8
#     shell:
#         """
# mkdir -p $(dirname {params.prefix})

# if [ -f {input}/.snakemake_timestamp ]; then
#     rm {input}/.snakemake_timestamp
# fi

# iqtree2 -T {threads} -s {input} --prefix {params.prefix} -bb 1000 -m LG+G4+I --undo > {log}

# cp {params.prefix}.contree {output.st}
# """
# -mset DCmut,JTTDCMut,LG,WAG,VT
# with -s it assumes a single model, with -p it will infer a concatenated partitioned alignment, this can get very slow with
# > 100 genes so I'll use -s for now


### ASTRAL

rule prepare_Astral:
    input:
        sc_labs=outdir+"data/single_copy/best_trees_{code}_sc_labs.nwk"
    output:
        astral_gt=outdir+"input/astral/phylome_{code}_input_astral.nwk"
    benchmark:
        outdir+"benchmarks/{code}_astral_prep.txt"
    params:
        collapse=config['collapse']
    shell:
        "nw_ed {input.sc_labs} 'i & b<={params.collapse}' o > {output.astral_gt}"


rule Astral:
    input:
        gt=outdir+"input/astral/phylome_{code}_input_astral.nwk"
    output:
        st=outdir+"sptrees/phylome_{code}_sptree_astral.nwk"
    log:
        outdir+"log/phy_{code}_astral.log"
    benchmark:
        outdir+"benchmarks/{code}_astral.txt"
    threads:
        8
    shell:
        """
astral -t {threads} -i {input.gt} -o {output.st} 2> {log}
"""
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

# echo "collapsing nodes with support < than {params.collapse}"
# nw_ed {input} 'i & b<={params.collapse}' o > {output.gt}
# java -jar {params.astral} -i {output.gt} -o {output.st}  2>&1 | tee {log}
# """
# this is astralIII 5.7.8
# -t 32 for discovista stuff

rule wASTRAL:
    input:
        gt=outdir+"data/single_copy/best_trees_{code}_sc_labs.nwk"
    output:
        # gt=outdir+"input/astral/phylome_{code}_input_wastral.nwk",
        st=outdir+"sptrees/phylome_{code}_sptree_wastral.nwk"
    log:
        outdir+"log/phy_{code}_wastral.log"
    benchmark:
        outdir+"benchmarks/{code}_wastral.txt"
    threads:
        8
    shell:
        """
astral-hybrid -x 1 -t {threads} -L -o {output.st} {input.gt} 2> {log}
"""


## MULTI COPY

### ASTEROID

rule prepare_Asteroid:
    input:
        mc=outdir+"data/all_trees/best_trees_{code}.nwk",
    output:
        asteroid_map=outdir+"input/asteroid/phylome_{code}_mapping_ad.txt"
    benchmark:
        outdir+"benchmarks/{code}_ad_prep.txt"
    shell:
        "python scripts/prepare_trees.py -i {input.mc} -o {output.asteroid_map} -m asteroid"

rule Asteroid:
    input:
        mc=outdir+"data/all_trees/best_trees_{code}.nwk",
        asteroid_map=outdir+"input/asteroid/phylome_{code}_mapping_ad.txt",
        sc=outdir+"data/single_copy/best_trees_{code}_sc_labs.nwk"
    output:
        st=outdir+"sptrees/phylome_{code}_sptree_asteroid.nwk"
    params:
        prefix=outdir+"sptrees/phylome_{code}_asteroid_tmp"
    benchmark:
        outdir+"benchmarks/{code}_asteroid.txt"
    threads:
        8
    log:
        outdir+"log/phy_{code}_asteroid.log"
    shell:
        """
set +e
mpiexec --oversubscribe -np {threads} asteroid -i {input.mc} -b 100 -p {params.prefix} -m {input.asteroid_map} &> {log}
exitcode=$?

if [ $exitcode -eq 134 ]; then
    echo "Asteroid failed with multi copy gene trees! The resulting species tree will be based on single copy gene trees" >> {log}
    mpiexec --oversubscribe -np {threads} asteroid -i {input.sc} -b 100 -p {params.prefix} &>> {log}
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
    benchmark:
        outdir+"benchmarks/{code}_apro.txt"
    threads:
        8
    shell:
        """
astral-pro -u 1 -t {threads} -o {output.st} {input} 2> {log}
"""


### ASTRAL-DISCO
rule prepare_Astral_DISCO:
    input:
        mc=outdir+"data/all_trees/best_trees_{code}.nwk"
    output:
        sc_disco=outdir+"data/all_trees/best_trees_{code}_disco.nwk",
        ad_gt=outdir+"input/ad/phylome_{code}_input_ad.nwk"
    benchmark:
        outdir+"benchmarks/{code}_ad_prep.txt"
    params:
        collapse=config['collapse']
    shell:
        """
python scripts/disco.py -i {input.mc} -o {output.sc_disco} -d \"_\" -n 1
nw_ed {output.sc_disco} 'i & b<={params.collapse}' o > {output.ad_gt}
"""

rule Astral_DISCO:
    input:
        gt=outdir+"input/ad/phylome_{code}_input_ad.nwk"
    output:
        st=outdir+"sptrees/phylome_{code}_sptree_ad.nwk"
    log:
        outdir+"log/phy_{code}_ad.log"
    benchmark:
        outdir+"benchmarks/{code}_ad.txt"
    threads:
        8
    shell:
        """
astral -u 1 -t {threads} -i {input.gt} -o {output.st} 2> {log}
"""

### SpeciesRax

rule prepare_SpeciesRax:
    input:
        mc_raw=outdir+"data/txts/best_trees_{code}.txt"
    output:
        sprax_fams=temp(directory(outdir+"input/sprax/phylome_{code}")),
    benchmark:
        outdir+"benchmarks/{code}_sprax_prep.txt"
    shell:
        "python scripts/prepare_trees.py -i {input.mc_raw} -o {output.sprax_fams} -m speciesrax"


rule SpeciesRax:
    input:
        sprax_fams=outdir+"input/sprax/phylome_{code}"
    output:
        sprax_res=temp(directory(outdir+"input/sprax/phylome_{code}_results")),
        st=outdir+"sptrees/phylome_{code}_sptree_sprax.nwk"
    log:
        outdir+"log/phy_{code}_sprax.log"
    benchmark:
        outdir+"benchmarks/{code}_sprax.txt"
    threads:
        8
    shell:
        """
mpiexec --oversubscribe -np {threads} generax --families {input.sprax_fams}/family.txt \
--strategy SKIP --si-strategy HYBRID --species-tree MiniNJ \
--rec-model UndatedDTL --per-family-rates --prune-species-tree \
--si-estimate-bl --si-quartet-support --prefix {output.sprax_res} &> {log}

cp {output.sprax_res}/species_trees/inferred_species_tree.newick {output.st}
"""

# an optional starting species tree
# a set of unrooted gene trees OR a set of multiple sequence alignments.
# the mapping between gene taxa and species taxa.

### duptree

rule prepare_Duptree:
    input:
        mc_labs=outdir+"data/all_trees/best_trees_{code}_labs.nwk"
    output:
        dt_gt=outdir+"input/dt/phylome_{code}_input_dt.nwk",
    benchmark:
        outdir+"benchmarks/{code}_dt_prep.txt"
    shell:
        "python scripts/prepare_trees.py -i {input.mc_labs} -o {output.dt_gt} -m duptree"


rule Duptree:
    input:
        gt=outdir+"input/dt/phylome_{code}_input_dt.nwk",
    output:
        st=outdir+"sptrees/phylome_{code}_sptree_dt.nwk"
    log:
        outdir+"log/phy_{code}_dt.log"
    benchmark:
        outdir+"benchmarks/{code}_dt.txt"
    shell:
        """
duptree -i {input.gt} -o {output.st} --nogenetree &> {log}
"""

### fastmulrfs

rule prepare_FastMulRFS:
    input:
        mc_labs=outdir+"data/all_trees/best_trees_{code}_labs.nwk"
    output:
        fmulrfs_gt=outdir+"input/fmulrfs/phylome_{code}_input_fmulrfs.nwk"
    benchmark:
        outdir+"benchmarks/{code}_fmulrfs_prep.txt"
    shell:
        "python scripts/fastmulrfs/python-tools/preprocess_multrees_v3.py -i {input.mc_labs} -o {output.fmulrfs_gt}"


rule FastMulRFS:
    input:
        gt=outdir+"input/fmulrfs/phylome_{code}_input_fmulrfs.nwk",
    output:
        st=outdir+"sptrees/phylome_{code}_sptree_fmulrfs.nwk"
    log:
        outdir+"log/phy_{code}_fmulrfs.log"
    benchmark:
        outdir+"benchmarks/{code}_fmulrfs.txt"
    shell:
        """
sp_prefix=$(echo {input.gt} | sed 's/_input_fmulrfs.nwk/_results_fmulrfs/')
./scripts/fastmulrfs/external/FastRFS/build/FastRFS -i {input.gt} -o $sp_prefix &> {log}
cp ${{sp_prefix}}.single {output.st}
"""

# QC

### Final sptrees

rule Get_consensus:
    input:
        sptrees=expand(outdir+"sptrees/phylome_{{code}}_sptree_{method}.nwk", method=methods),
        alns=outdir+"data/aln/aln_{code}_sc",
        sc_l=outdir+"data/single_copy/best_trees_{code}_sc_labs.nwk",
        mc_raw=outdir+"data/txts/best_trees_{code}.txt"
    output:
        contree=outdir+"contrees/phylome_{code}_contree.nwk",
        partition_file=outdir+"contrees/phylome_{code}_part.txt",
        st=outdir+"contrees/phylome_{code}_sptree.nwk"
    params:
        con_prefix=outdir+"contrees/phylome_{code}_consensus",
        bl_prefix=outdir+"contrees/phylome_{code}_bl",
    benchmark:
        outdir+"benchmarks/{code}_consensus.txt"
    log:
        outdir+"log/phy_{code}_consensus.log"
    threads:
        8
    shell:
        """
mkdir -p $(dirname {output.contree})
cat {input.sptrees} | nw_reroot -d - | nw_topology -I - > {output.contree}

iqtree2 -con -t {output.contree} --quiet
mv {output.contree}.contree {output.contree}
cat {output.contree}.log > {log}
rm {output.contree}.log

if [ -f {input.alns}/.snakemake_timestamp ]; then
    rm {input.alns}/.snakemake_timestamp
fi

iqtree2 -te {output.contree} --gcf {input.sc_l} -p {input.alns} --scf 100 --prefix {params.con_prefix} -T {threads} --quiet
cat {params.con_prefix}.log >> {log}
rm {params.con_prefix}.log {params.con_prefix}.cf.tree.nex {params.con_prefix}.cf.branch {params.con_prefix}.cf.stat

./scripts/create_part.sh {input.alns} {input.mc_raw} > {output.partition_file}
n_tosub=$(grep -v "end;" {output.partition_file} | wc -l)
sed -i "${{n_tosub}}s/,/;/" {output.partition_file}

iqtree2 -te {params.con_prefix}.cf.tree -s {input.alns} -p {output.partition_file} --prefix {params.bl_prefix} -T {threads} --safe --quiet
mv {params.bl_prefix}.treefile {output.st}
cat {params.bl_prefix}.log >> {log}
rm {params.con_prefix}.cf.tree {params.bl_prefix}.log {params.bl_prefix}.ckp.gz {params.bl_prefix}.iqtree
"""
# iqtree2 -p output/data/aln/aln_0096_sc/ -t output/results/phylome_0096_sptree.nwk
# this to infer branch lengths


rule Plot_results:
    input:
        sts=expand(outdir+"sptrees/phylome_{{code}}_sptree_{method}.nwk", method=methods),
        consensus=outdir+"contrees/phylome_{code}_sptree.nwk",
        info=outdir+"data/info/info_{code}.tsv"
    output:
        # directory(outdir+"results/phylome_{code}")
        outdir+"results/phylome_{code}/phylome_{code}_rooted_sptree.nwk"
    shell:
        """
Rscript scripts/compare_sptrees.R -i $(echo {input.sts} | sed -e 's/\s\+/,/g') -c {wildcards.code} -t {input.info} -o {output} -s {input.consensus}
"""

rule Plot_occupancy:
    input:
        trees=outdir+"data/all_trees/best_trees_{code}.nwk",
        sc=outdir+"data/single_copy/best_trees_{code}_sc.nwk",
        info=outdir+"data/info/info_{code}.tsv",
        consensus=outdir+"results/phylome_{code}/phylome_{code}_rooted_sptree.nwk"
    output:
        all_t=outdir+"results/phylome_{code}/phylome_{code}_all_presence.pdf",
        sc=outdir+"results/phylome_{code}/phylome_{code}_sc_presence.pdf",
        occ=outdir+"results/phylome_{code}/phylome_{code}_all_oc.pdf"
    params:
        mo=minimum_occupancy
    shell:
        """
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