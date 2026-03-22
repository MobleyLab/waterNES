Individual folders for each of the systems and the each folder is named as "system_ligand1_ligand2"

The systems are:
1. Scytalone Dehydratase (ligands C3D and C5D) Ben Shalom et al. JCTC 2020, 16, 12, 7883-7894 (https://pubs.acs.org/doi/full/10.1021/acs.jctc.0c00785)

2. Thrombin (Ligands B5 and B1a) Ben Shalom et al. JCTC 2020, 16, 12, 7883-7894 (https://pubs.acs.org/doi/full/10.1021/acs.jctc.0c00785)

3. Factor Xa (Ligands IID and IIE) Abel et al. JACS 2008, 130, 9, 2817-2831 (https://pubs.acs.org/doi/full/10.1021/ja0771033?casa_token=l9GPHbyOwykAAAAA%3A1lLA1RXM6JyFC1g5ulxeol-cin6VcLe99Z60thXad_xiKDPEq_TOHF9nf-VXUL937-On4X37JwerAQ)

4. BACE1 (Ligands C4j and C4d) Wahl et al. JCIM 2019, 59, 2, 754-765 (https://pubs.acs.org/doi/full/10.1021/acs.jcim.8b00826)

5. EGFR Kinase (Ligands 7 and 8) Michel et al. JACS 2009, 131, 42, 15403-15411 (https://pubs.acs.org/doi/10.1021/ja906058w)


In each of the folder, the files are:

1. system.pdb : This file contains the structure of the complex (protein and the hybrid ligand) as the starting point for my simulations. The complex structure is solvated in TIP3P water and NaCl ions are added for a concentration of 0.2 molar. The buried water is also present in the structure.

2. minimized.gro : The .gro file of an energy-minimized structure of the complex (protein and the hybrid ligand). It contains the attachment site for the buried water (residue ATT) as well as the buried water (residue MOL).

3. system.top : This file contains the GROMACS topology for the structure "minimized.gro"

4. hybrid_ligand.pdb : The structure of  the hybrid ligand.

5. hybrid_ligand.top : The topology of the hybrid ligand.

6. ligand1.pdb : The PDB structure for the ligand1 (that binds with the buried water).

7. ligand1.sdf: The SDF file of the structure for the ligand1 (that binds with the buried water).

8. ligand1.top : The GROMACS topology file for ligand1 (that binds with the buried water).

9. ligand2.pdb : The PDB structure for the ligand2 (that dispaces the buried water).

10. ligand2.sdf : The SDF file of the structure for the ligand2 (that dispaces the buried water).

11. ligand2.top : The GROMACS topology file for ligand2 (that dispaces the buried water).

12. posre_protein.itp : Position restraints file for the protein.

13. posre_ligand.itp : osition restraints file for the hybrid ligand.


