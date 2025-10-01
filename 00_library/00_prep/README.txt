Preparation of H-2Db structure for RIFdock

1. Download 1S7U crystal structure from RCSB
2. Remove crystallographic waters
3. Run ./bcov_fastrelax.py on the resulting structure (to add hydrogens, fix non-ideal bond lengths and angles, and potentially adjust sidechains to more favorable rotamers)
4. Delete B2m and a3 domains
5. renumber residues to start from 1, rename chain to chain A
6. Result is ./1S7U_relaxed_head.pdb