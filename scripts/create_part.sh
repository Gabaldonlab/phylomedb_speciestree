aln_dir=$1
bt_file=$2

echo "#nexus"
echo "begin sets;"

for aln in $aln_dir/*fasta; do
    # if echo $aln | grep -q 584; then
    #     type="CODON"
    # else
    type="AA"
    # fi
    len=$(seqkit fx2tab --length --name --header-line $aln | awk 'NR>1' | cut  -f2 | sort -u)
    echo "  charset $(basename $aln) = $aln: $type, 1-$len;"
done

echo "  charpartition mymodels ="

for aln in $aln_dir/*fasta; do
    gn=$(basename $aln '.clean.fasta')
    model=$(awk -v gn=$gn '$1==gn' $bt_file | cut -f2)
    if [ "$model" == "GY_F3X4" ]; then
        model=GY+F3X4
    fi
    echo "    $model: $(basename $aln),"
done

echo "end;"

# #nexus
# begin sets;
#   charset Phy004A5S0_224308.clean.fasta = alns/Phy004A5S0_224308.clean.fasta: , ;
#   charset Phy004A5S1_224308.clean.fasta = alns/Phy004A5S1_224308.clean.fasta: , ;
#   charset Phy004A5S3_224308.clean.fasta = alns/Phy004A5S3_224308.clean.fasta: , ;
#   charset Phy004A5S5_224308.clean.fasta = alns/Phy004A5S5_224308.clean.fasta: , ;
#   charset Phy004A5S6_224308.clean.fasta = alns/Phy004A5S6_224308.clean.fasta: , ;
#   charset Phy004A5SF_224308.clean.fasta = alns/Phy004A5SF_224308.clean.fasta: , ;
#   charset Phy004A5SL_224308.clean.fasta = alns/Phy004A5SL_224308.clean.fasta: , ;
#   charset Phy004A5SR_224308.clean.fasta = alns/Phy004A5SR_224308.clean.fasta: , ;
#   charset Phy004A5ST_224308.clean.fasta = alns/Phy004A5ST_224308.clean.fasta: , ;
#   charset Phy004A5SW_224308.clean.fasta = alns/Phy004A5SW_224308.clean.fasta: , ;
#   charset Phy004A5SX_224308.clean.fasta = alns/Phy004A5SX_224308.clean.fasta: , ;
#   charpartition mymodels =
#     LG: Phy004A5S0_224308.clean.fasta{8.41182},
#     LG: Phy004A5S1_224308.clean.fasta{9.67173},
#     LG: Phy004A5S3_224308.clean.fasta{2.73076},
#     LG: Phy004A5S5_224308.clean.fasta{2.85386},
#     LG: Phy004A5S6_224308.clean.fasta{5.99817},
#     LG: Phy004A5SF_224308.clean.fasta{3.2653},
#     LG: Phy004A5SL_224308.clean.fasta{3.34878},
#     LG: Phy004A5SR_224308.clean.fasta{8.12475},
#     LG: Phy004A5ST_224308.clean.fasta{3.30847},
#     LG: Phy004A5SW_224308.clean.fasta{3.58615},
#     LG: Phy004A5SX_224308.clean.fasta{6.00453};
# end;
