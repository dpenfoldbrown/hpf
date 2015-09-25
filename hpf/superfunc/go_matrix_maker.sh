#!/bin/bash

echo "Creating lots of matrices: GO parent and MGI->GO annotation"
echo ""

DBS="mygolite_012006
mygolite_012007
mygolite_012008
mygolite_012009
mygolite_012010
mygolite_012011
mygolite_012012"

PARENT_PY="/Users/dpb/bonneau-dev/hpf/trunk/src/hpf/superfunc/go_parent_matrix.py"
MGI_PY="/Users/dpb/bonneau-dev/hpf/trunk/src/hpf/superfunc/go_mgi_matrix.py"

MGI_FILE="/Users/dpb/bonneau-dev/hpf/trunk/src/projects/superfunc/mgiIds.txt"

BP_FILE="/Users/dpb/bonneau-dev/hpf/trunk/src/projects/superfunc/GO/goBP.txt"
CC_FILE="/Users/dpb/bonneau-dev/hpf/trunk/src/projects/superfunc/GO/goCC.txt"
MF_FILE="/Users/dpb/bonneau-dev/hpf/trunk/src/projects/superfunc/GO/goMF.txt"

OUT_LOC="/Users/dpb/Documents/superfunc/go/sparse_matrices"

# For each DB, create parent and mgi matrices for each GO cat. (BP, CC, and MF)
for db in $DBS
do
    echo "-<> Creating GO parent matrices for $db"
    python $PARENT_PY -g $BP_FILE -d $db --sparse -o $OUT_LOC/$db"_goBP_parent_sparse.txt"
    python $PARENT_PY -g $CC_FILE -d $db --sparse -o $OUT_LOC/$db"_goCC_parent_sparse.txt"
    python $PARENT_PY -g $MF_FILE -d $db --sparse -o $OUT_LOC/$db"_goMF_parent_sparse.txt"

    echo "-<> Creating MGI-GO annotation matrices for $db"
    python $MGI_PY -m $MGI_FILE -g $BP_FILE -d $db --sparse -o $OUT_LOC/$db"_goBP_mgi_sparse.txt"
    python $MGI_PY -m $MGI_FILE -g $CC_FILE -d $db --sparse -o $OUT_LOC/$db"_goCC_mgi_sparse.txt"
    python $MGI_PY -m $MGI_FILE -g $MF_FILE -d $db --sparse -o $OUT_LOC/$db"_goMF_mgi_sparse.txt"
done

echo ""
echo "Creating all matrices complete"
exit 0
