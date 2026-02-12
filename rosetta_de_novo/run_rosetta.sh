#!/usr/bin/env bash
set -euo pipefail

ROSETTA_MAIN="/home/vedansh/Downloads/rosetta_binary_ubuntu_3.15_bundle/rosetta.binary.ubuntu.release-408/main"
ROSETTA_DB="/home/vedansh/Downloads/rosetta_binary_ubuntu_3.15_bundle/rosetta.binary.ubuntu.release-408/main/database"
FLEXPEP="/home/vedansh/Downloads/rosetta_binary_ubuntu_3.15_bundle/rosetta.binary.ubuntu.release-408/main/source/bin/FlexPepDocking.static.linuxgccrelease"
ROSETTA_SCRIPTS="/home/vedansh/Downloads/rosetta_binary_ubuntu_3.15_bundle/rosetta.binary.ubuntu.release-408/main/source/bin/rosetta_scripts.static.linuxgccrelease"

RELAX_XML="/home/vedansh/projects/random/git-install/lib/relax_script.xml"
MMP8_END=161
PEP_START=162
PEP_END=179

mkdir -p mutated prepack refine relax

# Apply D-residue mutations to each start pose
# for f in starts/start_*.pdb; do
#   name=$(basename "$f" .pdb)
#   $ROSETTA_SCRIPTS -database "$ROSETTA_DB" -parser:protocol mutate_peptide.xml -s "$f" -out:path:all mutated -out:suffix _mut
# done

# # Prepack and refine (adjust -nstruct as needed)
# for f in starts/start_*.pdb; do
#   name=$(basename "$f" .pdb)
#   cst=constraints/constraints_${name#start_}.cst
#   input=$f
#   if [ -f mutated/${name}_mut_0001.pdb ]; then input=mutated/${name}_mut_0001.pdb; fi
#   input_base=$(basename "$input" .pdb)
#   $FLEXPEP -database "$ROSETTA_DB" -s "$input" -flexPepDocking:flexpep_prepack -out:path:all prepack -out:suffix _prepack
#   prepack=prepack/${input_base}_prepack_0001.pdb
#   $FLEXPEP -database "$ROSETTA_DB" -s "$prepack" -flexPepDocking:pep_refine -nstruct 200 -constraints:cst_file "$cst" -constraints:cst_fa_weight 1.0 -out:path:all refine -out:suffix _refine
# done

# Relax best models (edit input glob as needed)
for f in refine/*_refine_*.pdb; do
  $ROSETTA_SCRIPTS -database "$ROSETTA_DB" -s "/home/vedansh/projects/random/git-install/lib/rosetta_de_novo/refine/start_0010_mut_0001_prepack_0001_refine_0186.pdb" -parser:protocol "$RELAX_XML" -overwrite -parser:script_vars MMP8_END=$MMP8_END PEP_START=$PEP_START PEP_END=$PEP_END -out:path:all relax -out:suffix _relax
done
