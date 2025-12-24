var toolsData = [
  {
    "IDX": 1,
    "Tool Name": "query_uniprot",
    "Description": "Query UniProt for protein information.",
    "category": "Databases",
    "Server Name": "VenusFactory"
  },
  {
    "IDX": 2,
    "Tool Name": "query_interpro",
    "Description": "Query InterPro for protein domain information.",
    "category": "Databases",
    "Server Name": "VenusFactory"
  },
  {
    "IDX": 3,
    "Tool Name": "download_pdb_structure",
    "Description": "Download PDB structure file.",
    "category": "Databases",
    "Server Name": "VenusFactory"
  },
  {
    "IDX": 4,
    "Tool Name": "download_ncbi_sequence",
    "Description": "Download NCBI sequence file.",
    "category": "Databases",
    "Server Name": "VenusFactory"
  },
  {
    "IDX": 5,
    "Tool Name": "download_alphafold_structure",
    "Description": "Download AlphaFold structure file.",
    "category": "Databases",
    "Server Name": "VenusFactory"
  },
  {
    "IDX": 6,
    "Tool Name": "extract_pdb_sequence",
    "Description": "Extract sequence from PDB file.",
    "category": "Computational Tools",
    "Server Name": "VenusFactory"
  },
  {
    "IDX": 7,
    "Tool Name": "predict_zero_shot_sequence",
    "Description": "Predict zero-shot sequence.",
    "category": "Model Services",
    "Server Name": "VenusFactory"
  },
  {
    "IDX": 8,
    "Tool Name": "predict_zero_shot_structure",
    "Description": "Predict zero-shot structure.",
    "category": "Model Services",
    "Server Name": "VenusFactory"
  },
  {
    "IDX": 9,
    "Tool Name": "predict_protein_function",
    "Description": "Predict protein function.",
    "category": "Computational Tools",
    "Server Name": "VenusFactory"
  },
  {
    "IDX": 10,
    "Tool Name": "predict_functional_residue",
    "Description": "Predict functional residue.",
    "category": "Computational Tools",
    "Server Name": "VenusFactory"
  },
  {
    "IDX": 11,
    "Tool Name": "predict_protein_properties",
    "Description": "Predict protein properties.",
    "category": "Model Services",
    "Server Name": "VenusFactory"
  },
  {
    "IDX": 12,
    "Tool Name": "search_literature",
    "Description": "Search literature.",
    "category": "Literature Search",
    "Server Name": "VenusFactory"
  },
  {
    "IDX": 13,
    "Tool Name": "is_valid_protein_sequence",
    "Description": "Check if the input protein sequence string is valid.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 14,
    "Tool Name": "is_valid_smiles",
    "Description": "Check if the input SMILES string is valid",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 15,
    "Tool Name": "convert_smiles_to_other_format",
    "Description": "Convert a list of SMILES strings or a list of SMI file paths into other molecular formats.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 16,
    "Tool Name": "convert_pdb_to_pdbqt_dock",
    "Description": "Convert a protein file from PDB format to PDBQT format for docking preparation.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 17,
    "Tool Name": "convert_complex_cif_to_pdb",
    "Description": "Convert a protein-ligand complex file from CIF format to PDB format.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 18,
    "Tool Name": "visualize_protein",
    "Description": "Visualize the protein structure and save as an image.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 19,
    "Tool Name": "visualize_molecule",
    "Description": "Visualize the molecular structure and save as an image.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 20,
    "Tool Name": "visualize_complex",
    "Description": "Visualize the protein-ligand complex structure and save as an image.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 21,
    "Tool Name": "retrieve_protein_data_by_pdbcode",
    "Description": "Retrieve and download the protein sequence (.fasta) and structure (.pdb) files from RCSB using pdb code.",
    "category": "Databases",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 22,
    "Tool Name": "retrieve_smiles_from_name",
    "Description": "Retrieve SMILES strings from PubChem using compound names.",
    "category": "Databases",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 23,
    "Tool Name": "fix_pdb_dock",
    "Description": "Use PDBFixer to repair the protein structure PDB file in preparation for molecular docking. Note that use this tool to process the protein PDB file before performing molecular docking.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 24,
    "Tool Name": "read_smi_file",
    "Description": "Read the input smi file and extract the SMILES strings along with their corresponding compound names.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 25,
    "Tool Name": "read_fasta_file",
    "Description": "Parse protein sequences from the input fasta file.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 26,
    "Tool Name": "calculate_mol_basic_info",
    "Description": "Compute a set of basic molecular properties for each SMILES.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 27,
    "Tool Name": "calculate_mol_hydrophobicity",
    "Description": "Compute hydrophobicity-related molecular descriptors for each SMILES.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 28,
    "Tool Name": "calculate_mol_hbond",
    "Description": "Compute hydrogen bonding-related properties for each SMILES.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 29,
    "Tool Name": "calculate_mol_structure_complexity",
    "Description": "Compute a set of molecular complexity descriptors for each SMILES.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 30,
    "Tool Name": "calculate_mol_topology",
    "Description": "Compute a set of topological descriptors for each SMILES.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 31,
    "Tool Name": "calculate_mol_drug_chemistry",
    "Description": "Compute key drug-likeness metrics for each SMILES.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 32,
    "Tool Name": "calculate_mol_charge",
    "Description": "Compute Gasteiger partial charges and formal charge for each SMILES.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 33,
    "Tool Name": "calculate_mol_complexity",
    "Description": "Compute custom molecular complexity-related descriptors for each SMILES.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 34,
    "Tool Name": "calculate_protein_sequence_properties",
    "Description": "Compute a set of physicochemical properties for the input protein sequence.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 35,
    "Tool Name": "calculate_pdb_basic_info",
    "Description": "Read a protein pdb file and compute basic structural statistics.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 36,
    "Tool Name": "calculate_pdb_structural_geometry",
    "Description": "Read a protein pdb file and compute key geometric properties based on C¦Á atom coordinates.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 37,
    "Tool Name": "calculate_pdb_quality_metrics",
    "Description": "Read a protein pdb file and compute three key quality indicators.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 38,
    "Tool Name": "calculate_pdb_composition_info",
    "Description": "Read a protein pdb file and analyze compositional details by counting occurrences of each atom name, residue name, and the number of atoms per chain.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 39,
    "Tool Name": "calculate_smiles_similarity",
    "Description": "Compute the Tanimoto similarities between a target molecule and a list of candidate molecules using Morgan fingerprints.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 40,
    "Tool Name": "save_main_chain_pdb",
    "Description": "Extract the specified chain or the longest amino-acid chain from the input protein structure file and save as a new PDB file.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 41,
    "Tool Name": "extract_pdb_chains",
    "Description": "Extract the amino acid sequence of each chain from the PDB file.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 42,
    "Tool Name": "search_uniprot_id",
    "Description": "Search UniProt ID by gene name.",
    "category": "Databases",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 43,
    "Tool Name": "download_alphafold_structure",
    "Description": "Download predicted protein structures from the AlphaFold database.",
    "category": "Databases",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 44,
    "Tool Name": "small_file_to_base64",
    "Description": "Convert files smaller than 10MB to Base64 encoding.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 45,
    "Tool Name": "base64_to_file",
    "Description": "Convert Base64 encoding back to a file.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Tool"
  },
  {
    "IDX": 46,
    "Tool Name": "pred_protein_structure_esmfold",
    "Description": "Use the ESMFold model for protein 3D structure prediction.",
    "category": "Model Services",
    "Server Name": "DrugSDA-Model"
  },
  {
    "IDX": 47,
    "Tool Name": "pred_molecule_admet",
    "Description": "Predict the ADMET properties of a molecule.",
    "category": "Model Services",
    "Server Name": "DrugSDA-Model"
  },
  {
    "IDX": 48,
    "Tool Name": "pred_pocket_prank",
    "Description": "Use P2Rank to predict ligand binding pockets in the input protein.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Model"
  },
  {
    "IDX": 49,
    "Tool Name": "quick_molecule_docking",
    "Description": "Perform molecular docking using QuickVina2-GPU.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Model"
  },
  {
    "IDX": 50,
    "Tool Name": "calculate_dleps_score",
    "Description": "Enter a list of candidate small molecules. Based on the input disease name, identify upregulated and downregulated genes associated with the disease state, and predict a reversal score for each small molecule. Generally, a score above 0.2 indicates effectiveness, with higher scores being better.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Model"
  },
  {
    "IDX": 51,
    "Tool Name": "pred_mutant_sequence",
    "Description": "Given a protein sequence and its structure, employ the ProSST model to predict mutation effects and obtain the top-k mutated sequences based on their scores.",
    "category": "Model Services",
    "Server Name": "DrugSDA-Model"
  },
  {
    "IDX": 52,
    "Tool Name": "boltz_binding_affinity",
    "Description": "Use Boltz to predict binding affinity between protein (receptor) and small molecules (ligands).",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Model"
  },
  {
    "IDX": 53,
    "Tool Name": "run_fpocket",
    "Description": "Use fpocket to predict protein pockets and set it as the default tool for pocket prediction.",
    "category": "Computational Tools",
    "Server Name": "DrugSDA-Model"
  },
  {
    "IDX": 54,
    "Tool Name": "get_activity_by_id",
    "Description": "Retrieve the details of a single bioactivity entry from the ChEMBL database by its unique activity ID.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 55,
    "Tool Name": "get_activity_by_ids",
    "Description": "Retrieves a list of bioactivity entries from the ChEMBL database. Note: This tool currently does not function as expected. It does not correctly filter by the provided list of IDs and instead returns a default list of activities from the database.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 56,
    "Tool Name": "search_activity",
    "Description": "Performs a full-text search for bioactivity data in the ChEMBL database using a query string. This tool can search across various fields in the activity records, such as assay descriptions, target names, or molecule information.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 57,
    "Tool Name": "get_activity_supplementary_data_by_activity",
    "Description": "Retrieves a default list of supplementary bioactivity data from the ChEMBL database. Note: This tool is flawed and does not accept any filtering parameters, ignoring any provided input. To retrieve data for a specific activity, use the 'get_activity_supplementary_data_by_activity_by_id' tool instead.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 58,
    "Tool Name": "get_activity_supplementary_data_by_activity_by_id",
    "Description": "Retrieve single activitysupplementarydatabyactivity object details by ID.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 59,
    "Tool Name": "get_activity_supplementary_data_by_activity_by_ids",
    "Description": "Retrieve multiple activitysupplementarydatabyactivity objects by IDs.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 60,
    "Tool Name": "get_assay_by_name",
    "Description": "Retrieves the details for a single assay (experimental procedure) from the ChEMBL database using its name.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 61,
    "Tool Name": "get_assay_by_names",
    "Description": "Retrieves detailed information for multiple assays from the ChEMBL database using a list of their names.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 62,
    "Tool Name": "search_assay",
    "Description": "Performs a full-text search for assays (experimental procedures) in the ChEMBL database using a query string. This can search across various fields like the assay description.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 63,
    "Tool Name": "get_assay_class",
    "Description": "Retrieves a default list of assay classifications from the ChEMBL database. Note: This tool is flawed as it does not accept any filtering parameters and ignores any provided input.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 64,
    "Tool Name": "get_atc_class_by_level5",
    "Description": "Retrieves the details for a single ATC (Anatomical Therapeutic Chemical) classification from the ChEMBL database using the level 5 ATC code.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 65,
    "Tool Name": "get_atc_class_by_level5s",
    "Description": "Retrieves the details for multiple ATC (Anatomical Therapeutic Chemical) classifications from the ChEMBL database using a list of their level 5 ATC codes.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 66,
    "Tool Name": "get_binding_site_by_id",
    "Description": "Retrieves the details for a single binding site from the ChEMBL database using its unique integer ID.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 67,
    "Tool Name": "get_binding_site_by_ids",
    "Description": "Retrieves detailed information for multiple binding sites from the ChEMBL database using a list of their unique integer IDs.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 68,
    "Tool Name": "get_biotherapeutic_by_name",
    "Description": "Retrieves the details for a single biotherapeutic from the ChEMBL database using its name.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 69,
    "Tool Name": "get_biotherapeutic_by_names",
    "Description": "Retrieves detailed information for multiple biotherapeutics from the ChEMBL database using a list of their names.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 70,
    "Tool Name": "get_cell_line_by_id",
    "Description": "Retrieves the details for a single cell line from the ChEMBL database using its unique integer ID.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 71,
    "Tool Name": "get_cell_line_by_ids",
    "Description": "Retrieves detailed information for multiple cell lines from the ChEMBL database using a list of their unique integer IDs.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 72,
    "Tool Name": "get_compound_structural_alert",
    "Description": "Retrieve compound structural alert object list.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 73,
    "Tool Name": "get_compound_structural_alert_by_id",
    "Description": "Retrieve compound structural alert object details by ID.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 74,
    "Tool Name": "get_compound_structural_alert_by_ids",
    "Description": "Retrieve multiple compound structural alert objects by IDs.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 75,
    "Tool Name": "get_drug_by_name",
    "Description": "Retrieve single drug object details by name.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 76,
    "Tool Name": "get_drug_by_names",
    "Description": "Retrieve multiple drug objects by names.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 77,
    "Tool Name": "get_drug_indication_by_id",
    "Description": "Retrieve drug indication object details by ID",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 78,
    "Tool Name": "get_drug_indication_by_ids",
    "Description": "Retrieve multiple drug indication objects by IDs.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 79,
    "Tool Name": "get_drug_warning_by_id",
    "Description": "Retrieve single drug_warning object details by ID",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 80,
    "Tool Name": "get_drug_warning_by_ids",
    "Description": "Retrieve multiple drug_warning objects by IDs.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 81,
    "Tool Name": "get_go_slim_by_id",
    "Description": "Retrieves the details for a single GO (Gene Ontology) slim classification from the ChEMBL database using its unique GO ID.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 82,
    "Tool Name": "get_go_slim_by_ids",
    "Description": "Retrieves the details for multiple GO (Gene Ontology) slim classifications from the ChEMBL database using a list of their unique GO IDs.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 83,
    "Tool Name": "get_mechanism_by_id",
    "Description": "Retrieves the details for a single drug mechanism of action from the ChEMBL database using its unique integer ID.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 84,
    "Tool Name": "get_mechanism_by_ids",
    "Description": "Retrieves detailed information for multiple drug mechanisms of action from the ChEMBL database using a list of their unique integer IDs.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 85,
    "Tool Name": "get_metabolism_by_id",
    "Description": "Retrieve single metabolism object details by ID.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 86,
    "Tool Name": "get_metabolism_by_ids",
    "Description": "Retrieves detailed information for multiple drug metabolism records from the ChEMBL database using a list of their unique integer IDs.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 87,
    "Tool Name": "get_molecule_by_name",
    "Description": "Retrieves the details for a single molecule from the ChEMBL database using its name.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 88,
    "Tool Name": "get_molecule_by_names",
    "Description": "Retrieves detailed information for multiple molecules from the ChEMBL database using a list of their names.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 89,
    "Tool Name": "get_molecule_form",
    "Description": "Retrieves a default list of molecule forms from the ChEMBL database. Note: This tool is flawed as it does not accept any filtering parameters and ignores any provided input.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 90,
    "Tool Name": "get_molecule_form_by_name",
    "Description": "Retrieves molecule form information for a given molecule from the ChEMBL database using its name.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 91,
    "Tool Name": "get_molecule_form_by_names",
    "Description": "Retrieves molecule form information for multiple molecules from the ChEMBL database using a list of their names.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 92,
    "Tool Name": "get_organism_by_id",
    "Description": "Retrieves the details for a single organism from the ChEMBL database using its unique integer ID.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 93,
    "Tool Name": "get_organism_by_ids",
    "Description": "Retrieves detailed information for multiple organisms from the ChEMBL database using a list of their unique integer IDs.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 94,
    "Tool Name": "get_protein_classification_by_id",
    "Description": "Retrieves the details for a single protein classification from the ChEMBL database using its unique integer ID.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 95,
    "Tool Name": "get_protein_classification_by_ids",
    "Description": "Retrieves detailed information for multiple protein classifications from the ChEMBL database using a list of their unique integer IDs.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 96,
    "Tool Name": "search_protein_classification",
    "Description": "Search protein_classification object by query string.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 97,
    "Tool Name": "get_similarity_by_smiles",
    "Description": "Retrieve similarity details for compounds based on SMILES.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 98,
    "Tool Name": "get_source_by_id",
    "Description": "Retrieve single source object details by ID.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 99,
    "Tool Name": "get_source_by_ids",
    "Description": "Retrieve multiple source object details by IDs.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 100,
    "Tool Name": "get_substructure_by_smiles",
    "Description": "Retrieve substructure matches using SMILES.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 101,
    "Tool Name": "get_target_by_name",
    "Description": "Retrieve single target object details by name.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 102,
    "Tool Name": "get_target_by_names",
    "Description": "Retrieve multiple target objects by names",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 103,
    "Tool Name": "search_target",
    "Description": "Search target using query string.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 104,
    "Tool Name": "get_target_component",
    "Description": "Retrieve target_component object list",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 105,
    "Tool Name": "get_target_relation_by_related_name",
    "Description": "Retrieve single target_relation object details by related target name",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 106,
    "Tool Name": "get_target_relation_by_related_names",
    "Description": "Retrieve multiple target_relation objects by related target names",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 107,
    "Tool Name": "get_tissue_by_id",
    "Description": "Retrieve single tissue object details by ID.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 108,
    "Tool Name": "get_tissue_by_ids",
    "Description": "Retrieve multiple tissue object details by IDs.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 109,
    "Tool Name": "get_compound_chembl_id_by_name",
    "Description": "Get compound chembl id by name.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 110,
    "Tool Name": "get_target_chembl_id_by_name",
    "Description": "Get target chembl id by name.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 111,
    "Tool Name": "get_assay_chembl_id_by_name",
    "Description": "Get assay chembl id by name.",
    "category": "Databases",
    "Server Name": "Origene-ChEMBL"
  },
  {
    "IDX": 112,
    "Tool Name": "kegg_info",
    "Description": "This operation displays the database release information with statistics for the databases.\nExcept for kegg, genes and ligand, this operation also displays the list of linked databases that can be used in the link operation.\nDisplays the current statistics of a given database.",
    "category": "Databases",
    "Server Name": "Origene-KEGG"
  },
  {
    "IDX": 113,
    "Tool Name": "kegg_find",
    "Description": "KEGG find - Data search. Finds entries with matching query keywords or other query data in a given database.",
    "category": "Databases",
    "Server Name": "Origene-KEGG"
  },
  {
    "IDX": 114,
    "Tool Name": "kegg_list",
    "Description": "This operation can be used to obtain a list of all entries in each database.When the organism code is known, the second form can be used to obtain a list of organism-specific pathways.The third form is a similar option for brite hierarchies.The fourth form may be used to obtain a list of definitions for a given set of database entry identifiers. The maximum number of identifiers that can be given is 10.",
    "category": "Databases",
    "Server Name": "Origene-KEGG"
  },
  {
    "IDX": 115,
    "Tool Name": "kegg_get",
    "Description": "This operation retrieves given database entries in a flat file format or in other formats with options. Flat file formats are available for all KEGG databases except brite. The input is limited up to 10 entries.",
    "category": "Databases",
    "Server Name": "Origene-KEGG"
  },
  {
    "IDX": 116,
    "Tool Name": "kegg_conv",
    "Description": "This operation can be used to convert entry identifiers (accession numbers) of outside databases to KEGG identifiers, and vice versa. The first form allows database to database mapping, while the second form allows conversion of a selected number of entries. The database name \"genes\" may be used only in the second form.",
    "category": "Databases",
    "Server Name": "Origene-KEGG"
  },
  {
    "IDX": 117,
    "Tool Name": "kegg_link",
    "Description": "KEGG link - find related entries by using database cross-references.",
    "category": "Databases",
    "Server Name": "Origene-KEGG"
  },
  {
    "IDX": 118,
    "Tool Name": "mapping_identifiers",
    "Description": "Maps common protein names, synonyms and UniProt identifiers into STRING identifiers",
    "category": "Databases",
    "Server Name": "Origene-STRING"
  },
  {
    "IDX": 119,
    "Tool Name": "get_string_network_interaction",
    "Description": "Retrieve STRING interaction network for one or multiple proteins in various text formats.\nIt will tell you the combined score and all the channel specific scores for the set of proteins.\nYou can also extend the network neighborhood by setting \"add_nodes\", which will add, to your network, new interaction partners in order of their confidence.",
    "category": "Databases",
    "Server Name": "Origene-STRING"
  },
  {
    "IDX": 120,
    "Tool Name": "get_all_interaction_partners_of_the_protein_set",
    "Description": "This method provides the interactions between your provided set of proteins and all the other STRING proteins.\nAs STRING network usually has a lot of low scoring interactions, you may want to limit the number of retrieved interaction per protein using \"limit\" parameter.",
    "category": "Databases",
    "Server Name": "Origene-STRING"
  },
  {
    "IDX": 121,
    "Tool Name": "get_similarity_scores_of_the_protein_set",
    "Description": "STRING internally uses the Smith-Waterman bit scores as a proxy for protein homology.\nUsing this API you can retrieve these scores between the proteins in a selected species.\nThey are symmetric,meaning A->B is equal to B->A.\nThe bit score cut-off below which we do not store or report homology is 50.",
    "category": "Databases",
    "Server Name": "Origene-STRING"
  },
  {
    "IDX": 122,
    "Tool Name": "get_best_similarity_hits_between_species",
    "Description": "Retrieve the similarity from your input protein(s) to the best (most) similar protein from each STRING species.",
    "category": "Databases",
    "Server Name": "Origene-STRING"
  },
  {
    "IDX": 123,
    "Tool Name": "get_functional_enrichment",
    "Description": "STRING maps several databases onto its proteins, this includes: Gene Ontology, KEGG pathways, UniProt Keywords, PubMed publications, Pfam domains, InterPro domains, and SMART domains.\nThe STRING enrichment API method allows you to retrieve functional enrichment for any set of input proteins.\nIt will tell you which of your input proteins have an enriched term and the term's description.\nThe API provides the raw p-values, as well as, False Discovery Rate (B-H corrected p-values).",
    "category": "Databases",
    "Server Name": "Origene-STRING"
  },
  {
    "IDX": 124,
    "Tool Name": "get_functional_annotation",
    "Description": "STRING maps several databases onto its proteins, this includes: Gene Ontology, KEGG pathways, UniProt Keywords, PubMed publications, Pfam domains, InterPro domains, and SMART domains.",
    "category": "Databases",
    "Server Name": "Origene-STRING"
  },
  {
    "IDX": 125,
    "Tool Name": "get_ppi_enrichment",
    "Description": "Get protein-protein interaction enrichment for list of genes denoted by their STRING identifiers",
    "category": "Databases",
    "Server Name": "Origene-STRING"
  },
  {
    "IDX": 126,
    "Tool Name": "tavily_search",
    "Description": "Run the search engine with a given query, retrieving and filtering results.",
    "category": "Literature Search",
    "Server Name": "Origene-Search"
  },
  {
    "IDX": 127,
    "Tool Name": "jina_search",
    "Description": "Run the search engine with a given query, retrieving and filtering results.\nThis implements a two-phase retrieval approach:\n1. Get preview information for many results\n2. Filter the previews for relevance\n3. Get full content for only the relevant results",
    "category": "Literature Search",
    "Server Name": "Origene-Search"
  },
  {
    "IDX": 128,
    "Tool Name": "clinvar_search",
    "Description": "Providing Information of variants of human proteins, and whether these variants are benign, mutual or harmful.\nQuery Format: Accept semantic query. When querying Clinvar, the question must contain a complete protein sequence. Then all the variants related to it will be retrieved.\nFor example, you can query 'According to ClinVar, which variants to the following sequence (bracketed by xml tags) are likely to be benign? <sequence>MAAPILKDVVAYVEVWSSNGTENYS</sequence>'",
    "category": "Databases",
    "Server Name": "Origene-Search"
  },
  {
    "IDX": 129,
    "Tool Name": "gsea_search",
    "Description": "Database containing gene, gene set and collection of gene sets information.\nQuery Format: Accept semantic query. When querying GSEA, the question must contain a gene set name and a gene sets collection name. Then all the gene information in that set will be retrieved.\nFor example, you can query 'Which of the following genes is most likely contained in the gene set GCNP_SHH_UP_LATE.V1_UP,\n    which contains genes up-regulated in granule cell neuron precursors (GCNPs) after stimulation with Shh for 24h.\n    This gene set is a part of the C6 collection: oncogenic signature gene sets.'",
    "category": "Databases",
    "Server Name": "Origene-Search"
  },
  {
    "IDX": 130,
    "Tool Name": "mousemine_search",
    "Description": "A database containing the relationship between mouse gene sets and diseases.\nQuery Format: Accept semantic query. When querying MouseMine, the question must contain a disease code staring with 'MP:'. Then all the gene information related to it will be retrieved.\nFor example, you can query 'Which genes are most likely contained in the gene set MP_ABNORMAL_TUMOR_VASCULARIZATION, which contains mouse genes annotated to abnormal tumor vascularization (MP:0010144) retrieved from the Mouse Genome Informatics database via MouseMine'",
    "category": "Databases",
    "Server Name": "Origene-Search"
  },
  {
    "IDX": 131,
    "Tool Name": "ensemble_search",
    "Description": "Ensemble is a database containing the relationship between human genes and chr locations.\nQuery Format: Accept semantic query. When querying Ensemble, the question must contain a chr position name. Then all the gene information related to it will be retrieved.\nFor example, you can query 'Which human genes are located at chr10q22 according to Ensembl Release 110?'",
    "category": "Databases",
    "Server Name": "Origene-Search"
  },
  {
    "IDX": 132,
    "Tool Name": "pubmed_search",
    "Description": "Search PubMed for academic articles and retrieve abstracts.\nQuery Format: Accept semantic query. For example, you can query 'FGF23 excess cardiovascular complications'",
    "category": "Literature Search",
    "Server Name": "Origene-Search"
  },
  {
    "IDX": 133,
    "Tool Name": "search_pubchem_by_name",
    "Description": "Search PubChem for compounds matching a chemical name.Be aware though that matching chemical names to structure is an inexact science at best, and a name may often refer to more than one record.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 134,
    "Tool Name": "search_pubchem_by_smiles",
    "Description": "Search PubChem for compounds matching a SMILES string.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 135,
    "Tool Name": "get_pubchem_compound_by_cid",
    "Description": "Get detailed compound information by PubChem CID.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 136,
    "Tool Name": "search_pubchem_advanced",
    "Description": "Perform an advanced search on PubChem using a complex query.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 137,
    "Tool Name": "get_substance_by_sid",
    "Description": "Get substance information by PubChem SID.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 138,
    "Tool Name": "get_compound_by_cid",
    "Description": "Get compound information by PubChem CID.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 139,
    "Tool Name": "get_compound_by_name",
    "Description": "Get compound information by chemical name. Be aware though that matching chemical names to structure is an inexact science at best, and a name may often refer to more than one record.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 140,
    "Tool Name": "get_substance_by_name",
    "Description": "Get substance information by name.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 141,
    "Tool Name": "get_compound_property_by_name",
    "Description": "Retrieve a specific chemical property for a compound by its name from PubChem.\nThe response returns a table containing the requested property for the matching compound.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 142,
    "Tool Name": "get_compound_synonyms_by_name",
    "Description": "Retrieve all known synonyms for a given compound by its chemical name from PubChem.\nThe response returns a list of synonyms, including registry numbers, alternate names,\ntrade names, database identifiers, and systematic names.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 143,
    "Tool Name": "get_description_by_sid",
    "Description": "Get detailed description information for a PubChem substance given its SID.\nThe response includes the full Record structure, containing sections such as\n2D Structure, Identity, Source, External ID, Synonyms, Deposit/Modify Dates,\nStatus, and related compounds standardized from this substance.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 144,
    "Tool Name": "get_description_by_cid",
    "Description": "Retrieve detailed description information for a PubChem compound given its CID.\nThe response returns the full Record structure, including sections such as:\n? Structures (2D/3D depictions)\n? Identity (sources, external IDs, synonyms, versioning)\n? Deposit/Modify/Available dates\n? Status of the record\n? Related Records (related compounds, crystal data, articles, etc.)",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 145,
    "Tool Name": "get_general_info_by_compound_name",
    "Description": "Get detailed description of a compound by name, including overall information, drug and medication information, pharmacology and biochemistry information.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 146,
    "Tool Name": "get_description_by_aid",
    "Description": "Retrieve detailed description information for a PubChem bioassay given its AID.\nThe response returns the full Record structure, including sections such as:\n? Record Description (summary provided by depositor)\n? Description (bioassay overview and assay context)\n? Protocol (experimental protocol details)\n? Comment (depositor comments and references)\n? Result Definitions (data table column definitions)\n? Data Table (bioactivity results and flags)\n? Target (protein/gene targets)\n? Related Targets, Entrez Crosslinks\n? Identity (assay metadata: name, source, type, dates)\n? BioAssay Annotations (format, detection method, etc.)",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 147,
    "Tool Name": "get_assay_summary_by_cid",
    "Description": "Retrieve a summary of bioassay activities for a given PubChem compound CID.\nThe response includes a table of assays where the compound has been tested,\nwith details on activity outcome, target information, assay metadata, and references.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 148,
    "Tool Name": "get_assay_summary_by_sid",
    "Description": "Retrieve a summary of bioassay activities for a given PubChem substance SID.\nThe response includes a table listing assays in which the substance was tested,\nwith details on outcomes, target information, assay metadata, and references.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 149,
    "Tool Name": "get_gene_summary_by_geneid",
    "Description": "Get summary information for a gene by Gene ID.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 150,
    "Tool Name": "get_protein_summary_by_accession",
    "Description": "Retrieve summary information for a protein from PubChem by its UniProt or other accession number.\nThe response includes basic annotations such as the protein¡¯s name, taxonomy, and known synonyms.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 151,
    "Tool Name": "get_taxonomy_summary_by_taxonomyid",
    "Description": "Retrieve summary information for a biological taxonomy entry given its NCBI Taxonomy ID.\nThe response includes scientific and common names, rank, lineage, and synonyms.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 152,
    "Tool Name": "get_conformers_by_cid",
    "Description": "Retrieve available conformer identifiers for a given PubChem compound CID.\nConformers represent different 3D geometries computed or provided for the compound.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 153,
    "Tool Name": "get_compounds_by_smiles",
    "Description": "Retrieve compound objects from PubChem based on a given SMILES string.\nEach returned object represents a matching compound entry.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 154,
    "Tool Name": "get_compounds_by_formula",
    "Description": "Retrieve compound objects from PubChem based on a molecular formula.\nThe response returns a list of matching compound entries.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 155,
    "Tool Name": "get_molecular_formula",
    "Description": "Get the molecular formula of a compound.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 156,
    "Tool Name": "get_molecular_weight",
    "Description": "Get the molecular weight of a compound.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 157,
    "Tool Name": "get_isomeric_smiles",
    "Description": "Get the isomeric SMILES of a compound.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 158,
    "Tool Name": "get_xlogp",
    "Description": "Get the XLogP value of a compound.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 159,
    "Tool Name": "get_iupac_name",
    "Description": "Get the IUPAC name of a compound.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 160,
    "Tool Name": "get_synonyms",
    "Description": "Get the synonyms of a compound.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 161,
    "Tool Name": "get_cids_by_smiles",
    "Description": "Obtain the CID corresponding to the drug smiles",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 162,
    "Tool Name": "get_cids_by_formula",
    "Description": "Get a list of CIDs by molecular formula.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 163,
    "Tool Name": "get_sids_by_name",
    "Description": "Get a list of SIDs by name.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 164,
    "Tool Name": "get_substance_by_sid_pcp",
    "Description": "Get a Substance object by SID using PubChemPy.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 165,
    "Tool Name": "get_substances_by_name_pcp",
    "Description": "Get a list of Substance objects by name using PubChemPy.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 166,
    "Tool Name": "get_substances_source_id",
    "Description": "Get the source ID (Unique identifier assigned to the compound or substance by the original database (e.g. DrugBank, ChEMBL, etc.)) of a substance by SID.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 167,
    "Tool Name": "get_substances_synonyms",
    "Description": "Get the synonyms (Different names or identifiers for the same chemical substance) of a substance by SID.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 168,
    "Tool Name": "get_compound_dict",
    "Description": "Get a dictionary of a compound's properties.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 169,
    "Tool Name": "get_compounds_3d",
    "Description": "Get a list of compound objects with 3D structures.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 170,
    "Tool Name": "get_compounds_dict",
    "Description": "Get a dictionary representation of a compound by CID.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 171,
    "Tool Name": "get_substructure_cas",
    "Description": "Get CAS Registry Numbers for compounds containing a specified substructure.",
    "category": "Databases",
    "Server Name": "Origene-PubChem"
  },
  {
    "IDX": 172,
    "Tool Name": "get_gene_metadata_by_gene_name",
    "Description": "Get a gene summary by gene symbol. By default, in paged JSON format.\n\nThe returned JSON contains detailed gene information including:\n- Basic gene info: gene ID, symbol, description, tax ID, species name\n- Gene type and orientation\n- Reference standards and genomic locations\n- Chromosome location\n- External database IDs (HGNC, Swiss-Prot, Ensembl, OMIM)\n- Gene synonyms\n- Transcript and protein counts\n- Gene summary/description\n- Gene Ontology annotations:\n    - Molecular functions (e.g. DNA binding, transcription regulation)\n    - Biological processes (e.g. apoptosis, cell cycle regulation)\n    - Cellular components (e.g. nucleus, cytoplasm)",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 173,
    "Tool Name": "get_gene_by_ids",
    "Description": "Get gene information by gene IDs.",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 174,
    "Tool Name": "get_gene_by_accession",
    "Description": "Get gene information by accession.",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 175,
    "Tool Name": "get_gene_by_accession_dataset_report",
    "Description": "Get dataset reports by accession IDs",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 176,
    "Tool Name": "get_gene_by_accession_product_report",
    "Description": "Get gene product reports by accession IDs",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 177,
    "Tool Name": "get_gene_download_summary_by_id",
    "Description": "Get gene download summary by GeneID",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 178,
    "Tool Name": "get_gene_links_by_id",
    "Description": "Get gene links by gene ID",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 179,
    "Tool Name": "get_gene_dataset_report_by_locus_tag",
    "Description": "Get gene dataset reports by locus tag",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 180,
    "Tool Name": "get_gene_product_report_by_locus_tag",
    "Description": "Get gene product reports by locus tags",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 181,
    "Tool Name": "get_gene_by_symbol_dataset_report",
    "Description": "Get dataset reports by taxons",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 182,
    "Tool Name": "get_gene_by_symbol_product_report",
    "Description": "Get product reports by taxon",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 183,
    "Tool Name": "get_gene_by_taxon_dataset_report",
    "Description": "Get gene dataset reports by taxonomic identifier",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 184,
    "Tool Name": "get_gene_by_taxon_product_report",
    "Description": "Get gene product reports by taxonomic identifier",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 185,
    "Tool Name": "get_gene_dataset_report_by_id",
    "Description": "Get gene information by dataset report",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 186,
    "Tool Name": "get_genome_annotation_report",
    "Description": "Get genome annotation reports by genome accession.",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 187,
    "Tool Name": "get_genome_annotation_summary",
    "Description": "Get genome annotation report summary information.",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 188,
    "Tool Name": "get_genome_revision_history",
    "Description": "Get revision history for assembly by accession",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 189,
    "Tool Name": "get_genome_sequence_reports",
    "Description": "Get sequence reports by accessions",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 190,
    "Tool Name": "check_genome_accessions",
    "Description": "Check the validity of genome accessions",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 191,
    "Tool Name": "get_genome_dataset_report_by_accession",
    "Description": "Get dataset reports by accessions",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 192,
    "Tool Name": "get_genome_download",
    "Description": "Get a genome dataset by accession",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 193,
    "Tool Name": "get_genome_download_summary",
    "Description": "Preview genome dataset download",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 194,
    "Tool Name": "get_genome_links",
    "Description": "Get assembly links by accessions",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 195,
    "Tool Name": "get_genome_dataset_report_by_assembly_name",
    "Description": "Get dataset reports by assembly name",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 196,
    "Tool Name": "get_genome_dataset_report_by_bioproject",
    "Description": "Get dataset reports by bioproject",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 197,
    "Tool Name": "get_genome_dataset_report_by_biosample",
    "Description": "Get dataset reports by biosample id",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 198,
    "Tool Name": "get_sequence_assemblies",
    "Description": "Get assembly accessions for a sequence accession",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 199,
    "Tool Name": "get_genome_dataset_report_by_taxon",
    "Description": "Get dataset reports by taxons",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 200,
    "Tool Name": "get_genome_dataset_report_by_wgs",
    "Description": "Get dataset reports by wgs accession",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 201,
    "Tool Name": "get_virus_annotation_report",
    "Description": "Get virus annotation report by accessions.",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 202,
    "Tool Name": "check_virus_accessions",
    "Description": "Check virus accessions validity",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 203,
    "Tool Name": "get_virus_dataset_report",
    "Description": "Get virus dataset report by accessions",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 204,
    "Tool Name": "get_virus_by_taxon_annotation_report",
    "Description": "Get virus annotation report by taxon",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 205,
    "Tool Name": "get_virus_by_taxon_genome",
    "Description": "Get virus genome by taxon",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 206,
    "Tool Name": "get_virus_by_taxon_genome_table",
    "Description": "Get virus genome table by taxon",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 207,
    "Tool Name": "get_taxonomy_related_ids",
    "Description": "Get related taxonomy IDs",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 208,
    "Tool Name": "get_taxonomy_links",
    "Description": "Get taxonomy links",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 209,
    "Tool Name": "get_taxonomy",
    "Description": "Get taxonomy information.",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 210,
    "Tool Name": "get_taxonomy_dataset_report",
    "Description": "Get taxonomy dataset report.",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 211,
    "Tool Name": "get_taxonomy_filtered_subtree",
    "Description": "Get filtered taxonomy subtree",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 212,
    "Tool Name": "get_taxonomy_name_report",
    "Description": "Get taxonomy name report",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 213,
    "Tool Name": "get_taxonomy_taxon_suggest",
    "Description": "Get taxonomy suggestions",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 214,
    "Tool Name": "get_biosample_report",
    "Description": "Get biosample report.",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 215,
    "Tool Name": "get_organelle_dataset_report",
    "Description": "Get organelle dataset report.",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 216,
    "Tool Name": "get_organelle_by_taxon_dataset_report",
    "Description": "Get organelle dataset report by taxon.",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 217,
    "Tool Name": "get_gene_product_report_by_id",
    "Description": "Get gene product report by gene ID",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 218,
    "Tool Name": "get_gene_orthologs",
    "Description": "Get gene orthologs by gene ID",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 219,
    "Tool Name": "get_gene_by_taxon",
    "Description": "Get gene information by taxon",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 220,
    "Tool Name": "get_gene_counts_by_taxon",
    "Description": "Get gene counts by taxon",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 221,
    "Tool Name": "get_chromosome_summary",
    "Description": "Get chromosome summary by taxon and annotation name",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 222,
    "Tool Name": "get_genome_by_accession",
    "Description": "Get genome information by accession",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 223,
    "Tool Name": "get_prokaryote_gene_dataset_by_refseq_protein_accession",
    "Description": "Get a prokaryote gene dataset by RefSeq protein accession",
    "category": "Databases",
    "Server Name": "Origene-NCBI"
  },
  {
    "IDX": 224,
    "Tool Name": "get_general_info_by_protein_or_gene_name",
    "Description": "Get general information of a protein or gene by name from UniProt database.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 225,
    "Tool Name": "get_protein_sequence_by_name",
    "Description": "Get the human protein sequence by protein name.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 226,
    "Tool Name": "get_uniprotkb_entry_by_accession",
    "Description": "Search UniProtKB by protein entry accession to return all data associated with that entry.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 227,
    "Tool Name": "stream_uniprotkb_entries",
    "Description": "Stream all UniProtKB entries associated with the search term in a single download.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 228,
    "Tool Name": "search_uniprotkb_entries",
    "Description": "Search UniProtKB entries using a query, returns paginated list.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 229,
    "Tool Name": "get_uniref_cluster_by_id",
    "Description": "Search UniRef entry by id to return all data associated with that entry.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 230,
    "Tool Name": "get_uniref_cluster_members_by_id",
    "Description": "Search UniRef entry by member id to return all data associated with that entry.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 231,
    "Tool Name": "get_uniref_light_cluster_by_id",
    "Description": "Search light UniRef entry by id to return all data associated with that entry.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 232,
    "Tool Name": "stream_uniref_clusters",
    "Description": "Stream all UniRef clusters associated with the search term in a single download.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 233,
    "Tool Name": "search_uniref_clusters",
    "Description": "Search UniRef clusters using a query, returns paginated list.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 234,
    "Tool Name": "get_uniparc_entry_by_upi",
    "Description": "Search UniParc entry by id (UPI) to return all data associated with that entry.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 235,
    "Tool Name": "get_uniparc_light_entry_by_upi",
    "Description": "Search UniParc entry by id (UPI) to return all data associated with that entry (light version).",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 236,
    "Tool Name": "get_uniparc_cross_references_by_upi",
    "Description": "Get a page of database cross-reference entries by a UPI.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 237,
    "Tool Name": "stream_uniparc_cross_references_by_upi",
    "Description": "Stream database cross-reference entries for a specified UniParc UPI.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 238,
    "Tool Name": "stream_uniparc_entries",
    "Description": "Stream all UniParc entries associated with the specified search term in a single download.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 239,
    "Tool Name": "search_uniparc_entries",
    "Description": "Search UniParc entries using a query, returns paginated list.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 240,
    "Tool Name": "get_gene_centric_by_accession",
    "Description": "Retrieve a GeneCentric entry by UniProtKB accession.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 241,
    "Tool Name": "get_gene_centric_by_proteome",
    "Description": "Search GeneCentric entry by Proteome ID to return all data associated with that entry.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 242,
    "Tool Name": "stream_gene_centric",
    "Description": "Stream GeneCentric entries matching a query (max 10M entries).",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 243,
    "Tool Name": "search_gene_centric",
    "Description": "Search GeneCentric entries with pagination.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 244,
    "Tool Name": "get_proteome_by_id",
    "Description": "Retrieve a proteome by UniProt Proteome ID.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 245,
    "Tool Name": "stream_proteomes",
    "Description": "Stream Proteome entries matching a query (max 10M entries).",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 246,
    "Tool Name": "search_proteomes",
    "Description": "Search Proteome entries with pagination.",
    "category": "Databases",
    "Server Name": "Origene-UniProt"
  },
  {
    "IDX": 247,
    "Tool Name": "get_gene_specific_expression_in_cancer_type",
    "Description": "Analyze the tissue-specific expression pattern of a given gene across different cancer types\nusing the Firebrowse API (TCGA mRNASeq). It computes mean expression per cancer (cohort),\ncalculates z-score of mean expression, and returns cancer types where the gene is highly\nor lowly expressed.",
    "category": "Databases",
    "Server Name": "Origene-TCGA"
  },
  {
    "IDX": 248,
    "Tool Name": "get_lookup_symbol",
    "Description": "Look up Ensembl gene information by external gene symbol.\n\nThis function allows you to find Ensembl gene records using standard gene symbols\n(e.g., HGNC symbols for human genes).",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 249,
    "Tool Name": "get_homology_symbol",
    "Description": "Find evolutionary homologs (orthologs and paralogs) for a gene identified by symbol.\n\nRetrieve homologous genes across different species, with alignment statistics and\ntaxonomic information. Essential for comparative genomics and evolutionary studies.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 250,
    "Tool Name": "get_sequence_region",
    "Description": "Retrieve genomic DNA sequence for a specific chromosomal region.\n\nExtract the raw nucleotide sequence from a specific genomic location, which can be\nused for primer design, variant analysis, or sequence feature identification.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 251,
    "Tool Name": "get_vep_hgvs",
    "Description": "Predict the functional effects of variants using Variant Effect Predictor (VEP) with HGVS notation.\n\nAnalyzes the molecular consequences of genetic variants on genes, transcripts, and protein sequences.\nVEP provides comprehensive annotation including protein changes, regulatory effects, conservation scores,\nand pathogenicity predictions.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 252,
    "Tool Name": "get_genetree_id",
    "Description": "Retrieve a phylogenetic gene tree by its Ensembl stable identifier.\n\nGene trees represent the evolutionary history of genes across species, showing\northologous and paralogous relationships. These trees are constructed using\nprotein sequence alignments and phylogenetic algorithms.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 253,
    "Tool Name": "get_info_assembly",
    "Description": "Retrieve genome assembly information for a species.\n\nProvides details about the reference genome assembly used in Ensembl, including\nassembly version, accession numbers, and overall structure. Essential for\nunderstanding coordinate systems and genome organization.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 254,
    "Tool Name": "get_xrefs_symbol",
    "Description": "Get cross references for a gene symbol.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 255,
    "Tool Name": "get_archive_id",
    "Description": "Get the latest version of an Ensembl stable identifier.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 256,
    "Tool Name": "post_archive_id",
    "Description": "Get the latest version for multiple Ensembl stable identifiers.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 257,
    "Tool Name": "get_cafe_genetree_id",
    "Description": "Retrieve a CAFE (Computational Analysis of gene Family Evolution) gene tree by ID.\n\nCAFE analyzes the evolution of gene family size across a phylogenetic tree,\nidentifying expansions and contractions of gene families throughout evolution.\nThis helps understand adaptation, functional diversification, and species-specific traits.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 258,
    "Tool Name": "get_cafe_genetree_member_symbol",
    "Description": "Retrieve a CAFE gene tree for a gene identified by symbol.\n\nGet gene family evolution analysis (expansions/contractions) for a gene family\nthat contains the specified gene. Identifies evolutionary patterns without\nneeding to know the specific gene tree ID.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 259,
    "Tool Name": "get_cafe_genetree_member_id",
    "Description": "Retrieve the gene tree containing a gene identified by its Ensembl ID.\n\nFind the phylogenetic tree showing evolutionary relationships for a gene of interest,\nidentified using its Ensembl stable identifier. Useful for understanding gene evolution\nwhen you have the specific gene ID.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 260,
    "Tool Name": "get_genetree_member_symbol",
    "Description": "Get gene tree by symbol.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 261,
    "Tool Name": "get_alignment_region",
    "Description": "Retrieve genomic alignments between species for a specific region.\n\nGet multiple sequence alignments of genomic regions across species, showing\nevolutionary conservation and divergence. Crucial for identifying conserved\nfunctional elements like enhancers or detecting selection pressure.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 262,
    "Tool Name": "get_homology_id",
    "Description": "Find evolutionary homologs (orthologs and paralogs) for a gene identified by Ensembl ID.\n\nRetrieve homologous genes across different species, with alignment statistics and\ntaxonomic information. Essential for comparative genomics and evolutionary studies.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 263,
    "Tool Name": "get_xrefs_id",
    "Description": "Get cross references by ID.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 264,
    "Tool Name": "get_xrefs_name",
    "Description": "Get cross references by name.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 265,
    "Tool Name": "get_info_analysis",
    "Description": "List the analyses and data processing pipelines used for a species genome.\n\nProvides information about the computational methods and analyses used to generate\nEnsembl data for a species, including gene annotation methods, comparative genomics\nanalyses, and variation data processing.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 266,
    "Tool Name": "get_assembly_region_info",
    "Description": "Retrieve detailed information about a specific genomic region or chromosome.\n\nGet assembly metadata for a particular sequence within a genome assembly,\nsuch as chromosome length, scaffold composition, or contig information.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 267,
    "Tool Name": "get_info_biotypes",
    "Description": "Retrieve the catalog of gene and transcript biotypes for a species.\n\nBiotypes classify genes and transcripts according to their biological nature,\nsuch as protein-coding, pseudogene, or various non-coding RNA categories.\nThis information is crucial for filtering and interpreting genomic data.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 268,
    "Tool Name": "get_info_external_dbs",
    "Description": "Get external databases for a species.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 269,
    "Tool Name": "get_map_cdna",
    "Description": "Map cDNA coordinates to genomic coordinates.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 270,
    "Tool Name": "get_map_cds",
    "Description": "Map CDS coordinates to genomic coordinates.",
    "category": "Computational Tools",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 271,
    "Tool Name": "get_map_translation",
    "Description": "Map protein coordinates to genomic coordinates.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 272,
    "Tool Name": "get_ontology_ancestors",
    "Description": "Get ontology ancestors. Note: This tool is sensitive to the format of the input ID and may return a 400 Bad Request error for some valid-looking IDs. It is recommended to use IDs obtained directly from other Ensembl tools.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 273,
    "Tool Name": "get_ontology_descendants",
    "Description": "Get ontology descendants. Note: This tool is sensitive to the format of the input ID and may return a 400 Bad Request error for some valid-looking IDs. It is recommended to use IDs obtained directly from other Ensembl tools.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 274,
    "Tool Name": "get_ontology_id",
    "Description": "Get ontology by ID. Note: This tool is sensitive to the format of the input ID and may return a 400 Bad Request error for some valid-looking IDs. It is recommended to use IDs obtained directly from other Ensembl tools.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 275,
    "Tool Name": "get_ontology_name",
    "Description": "Get ontology by name.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 276,
    "Tool Name": "get_overlap_id",
    "Description": "Get features overlapping a region defined by an identifier. Note: This tool is currently non-functional and returns a 400 Bad Request error for valid Ensembl IDs.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 277,
    "Tool Name": "get_overlap_region",
    "Description": "Get features overlapping a genomic region. Note: This tool may fail with a 400 Bad Request error for valid queries.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 278,
    "Tool Name": "get_overlap_translation",
    "Description": "Get features overlapping a translation.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 279,
    "Tool Name": "get_phenotype_region",
    "Description": "Retrieve phenotype associations for variants in a genomic region.\n\nFind diseases, traits, and phenotypes associated with genetic variants\nlocated within a specific genomic region. Useful for exploring disease\nassociations in GWAS loci or candidate regions.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 280,
    "Tool Name": "get_phenotype_gene",
    "Description": "Retrieve phenotype associations for a specific gene.\n\nFind diseases, traits, and phenotypes associated with a gene of interest.\nThese associations come from various sources including literature curation,\nGWAS studies, and clinical databases.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 281,
    "Tool Name": "get_phenotype_accession",
    "Description": "Retrieve genomic features associated with a specific phenotype ontology term.\n\nFind genes and variants linked to a specific disease or trait identified by\nan ontology accession (e.g., from the Human Phenotype Ontology or Experimental\nFactor Ontology).",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 282,
    "Tool Name": "get_sequence_id",
    "Description": "Retrieve sequence associated with an Ensembl identifier.\n\nGet the nucleotide sequence for a gene or transcript, or the amino acid sequence for a protein.\nUseful for analyzing gene structure, transcript variants, or protein domains.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 283,
    "Tool Name": "get_vep_id",
    "Description": "Predict the functional effects of variants using Variant Effect Predictor (VEP) with variant identifier.\n\nRetrieves comprehensive variant annotation using known variant IDs (e.g., dbSNP rs identifiers).\nProvides molecular consequences, population frequencies, and pathogenicity predictions.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 284,
    "Tool Name": "get_vep_region",
    "Description": "Predict the functional effects of variants using Variant Effect Predictor (VEP) with genomic coordinates.\n\nAnalyzes variants specified by chromosome location and alternate allele.\nParticularly useful for novel variants or those without established identifiers.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 285,
    "Tool Name": "get_variation",
    "Description": "Retrieve detailed information about a genetic variant by its identifier.\n\nProvides comprehensive data about a known genetic variant, including its\ngenomic location, alleles, frequency in populations, phenotype associations,\nand links to external databases.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 286,
    "Tool Name": "get_variant_recoder",
    "Description": "Translate between different variant nomenclature systems and representations.\n\nConvert variant identifiers between different formats (e.g., rsID, HGVS notation,\ngenomic coordinates). Useful for integrating variant data from different sources\nor analysis tools.",
    "category": "Computational Tools",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 287,
    "Tool Name": "get_info_genomes",
    "Description": "Find information about a given genome.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 288,
    "Tool Name": "get_info_genomes_accession",
    "Description": "Find information about genomes containing a specified INSDC accession. Note: The underlying data is sparse, and many valid accessions may result in a null return.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 289,
    "Tool Name": "get_info_genomes_assembly",
    "Description": "Find information about a genome with a specified assembly. Note: This tool may fail with a 400 Bad Request error for valid assembly IDs.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 290,
    "Tool Name": "get_info_genomes_division",
    "Description": "Find information about all genomes in a given division.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 291,
    "Tool Name": "get_info_genomes_taxonomy",
    "Description": "Find information about all genomes beneath a given node of the taxonomy.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 292,
    "Tool Name": "get_info_variation",
    "Description": "List all variation data sources used for a species in Ensembl.\n\nProvides information about the databases, studies, and projects that contributed\nvariation data (SNPs, indels, structural variants) to Ensembl for a species.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 293,
    "Tool Name": "get_info_variation_consequence_types",
    "Description": "Lists all variant consequence types used by Ensembl.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 294,
    "Tool Name": "get_info_variation_populations",
    "Description": "List all variation populations for a species, or list all individuals in a specific population.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 295,
    "Tool Name": "get_ld",
    "Description": "Computes and returns LD values between the given variant and all other variants in a window.",
    "category": "Computational Tools",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 296,
    "Tool Name": "get_ld_pairwise",
    "Description": "Computes and returns LD values between the given variants.",
    "category": "Computational Tools",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 297,
    "Tool Name": "get_ld_region",
    "Description": "Computes and returns LD values between all pairs of variants in the defined region.",
    "category": "Computational Tools",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 298,
    "Tool Name": "get_lookup_id",
    "Description": "Look up details for any Ensembl stable identifier.\n\nRetrieve comprehensive information about any Ensembl entity (gene, transcript, protein, etc.)\nusing its stable identifier.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 299,
    "Tool Name": "post_lookup_id",
    "Description": "Look up details for multiple Ensembl stable identifiers in a single request.\n\nBatch retrieval of information for multiple Ensembl entities (genes, transcripts, proteins, etc.).",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 300,
    "Tool Name": "post_lookup_symbol",
    "Description": "Look up multiple gene symbols in a single request.\n\nBatch retrieval of Ensembl gene information for multiple external gene symbols.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 301,
    "Tool Name": "get_map",
    "Description": "Map coordinates between assemblies.",
    "category": "Computational Tools",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 302,
    "Tool Name": "get_ontology_ancestors_chart",
    "Description": "Reconstruct the entire ancestry of a term from is_a and part_of relationships.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 303,
    "Tool Name": "get_taxonomy_classification",
    "Description": "Return the taxonomic classification of a taxon node.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 304,
    "Tool Name": "get_taxonomy_id",
    "Description": "Search for a taxonomic term by its identifier or name",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 305,
    "Tool Name": "get_taxonomy_name",
    "Description": "Search for a taxonomic id by a non-scientific name.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 306,
    "Tool Name": "get_species_binding_matrix",
    "Description": "Return the specified binding matrix",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 307,
    "Tool Name": "post_sequence_id",
    "Description": "Request multiple types of sequence by a stable identifier list.\n\nEfficiently fetch sequences for multiple genes, transcripts, or proteins in a single request.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 308,
    "Tool Name": "post_sequence_region",
    "Description": "Get sequences by multiple regions.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 309,
    "Tool Name": "get_transcript_haplotypes",
    "Description": "Computes observed transcript haplotype sequences based on phased genotype data.",
    "category": "Computational Tools",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 310,
    "Tool Name": "post_vep_hgvs",
    "Description": "Batch predict the functional effects of multiple variants using VEP with HGVS notation.\n\nEfficiently analyze multiple variants in a single request using the Variant Effect Predictor.\nIdeal for analyzing sets of variants from sequencing data or genetic studies.",
    "category": "Computational Tools",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 311,
    "Tool Name": "post_vep_id",
    "Description": "Batch predict the functional effects of multiple variants using VEP with variant identifiers.\n\nEfficiently analyze multiple known variants in a single request using the Variant Effect Predictor.\nIdeal for analyzing sets of common variants or SNP panel data.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 312,
    "Tool Name": "post_vep_region",
    "Description": "Get variant effect predictions by multiple regions.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 313,
    "Tool Name": "post_variant_recoder",
    "Description": "Translate a list of variant identifiers, HGVS notations or genomic SPDI notations to all possible variant IDs, HGVS and genomic SPDI",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 314,
    "Tool Name": "get_variation_pmcid",
    "Description": "Fetch variants by publication using PubMed Central reference number (PMCID)",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 315,
    "Tool Name": "get_variation_pmid",
    "Description": "Fetch variants by publication using PubMed reference number (PMID)",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 316,
    "Tool Name": "post_variation",
    "Description": "Uses a list of variant identifiers (e.g. rsID) to return the variation features including optional genotype, phenotype and population data",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 317,
    "Tool Name": "get_ga4gh_beacon",
    "Description": "Get Beacon information.\n\nReturns:\n    Dictionary containing Beacon information.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 318,
    "Tool Name": "get_ga4gh_beacon_query",
    "Description": "Query Beacon.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 319,
    "Tool Name": "post_ga4gh_beacon_query",
    "Description": "Query Beacon with POST.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 320,
    "Tool Name": "get_ga4gh_features",
    "Description": "Get GA4GH features by ID.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 321,
    "Tool Name": "post_ga4gh_features_search",
    "Description": "Get a list of sequence annotation features in GA4GH format",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 322,
    "Tool Name": "post_ga4gh_callsets_search",
    "Description": "Search GA4GH callsets.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 323,
    "Tool Name": "get_ga4gh_callsets",
    "Description": "Get the GA4GH record for a specific CallSet given its identifier",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 324,
    "Tool Name": "post_ga4gh_datasets_search",
    "Description": "Get a list of datasets in GA4GH format",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 325,
    "Tool Name": "get_ga4gh_datasets",
    "Description": "Get the GA4GH record for a specific dataset given its identifier",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 326,
    "Tool Name": "post_ga4gh_featuresets_search",
    "Description": "Search GA4GH feature sets.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 327,
    "Tool Name": "get_ga4gh_featuresets",
    "Description": "Return the GA4GH record for a specific featureSet given its identifier",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 328,
    "Tool Name": "get_ga4gh_variants",
    "Description": "Get GA4GH variant by ID.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 329,
    "Tool Name": "post_ga4gh_variantannotations_search",
    "Description": "Return variant annotation information in GA4GH format for a region on a reference sequence",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 330,
    "Tool Name": "post_ga4gh_variants_search",
    "Description": "Return variant call information in GA4GH format for a region on a reference sequence",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 331,
    "Tool Name": "post_ga4gh_variantsets_search",
    "Description": "Search GA4GH variant sets.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 332,
    "Tool Name": "get_ga4gh_variantsets",
    "Description": "Return the GA4GH record for a specific VariantSet given its identifier",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 333,
    "Tool Name": "post_ga4gh_references_search",
    "Description": "Return a list of reference sequences in GA4GH format",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 334,
    "Tool Name": "get_ga4gh_references",
    "Description": "Return data for a specific reference in GA4GH format by id",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 335,
    "Tool Name": "post_ga4gh_referencesets_search",
    "Description": "Search GA4GH reference sets.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 336,
    "Tool Name": "get_ga4gh_referencesets",
    "Description": "Search data for a specific reference set in GA4GH format by ID",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 337,
    "Tool Name": "post_ga4gh_variantannotationsets_search",
    "Description": "Return a list of annotation sets in GA4GH format",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 338,
    "Tool Name": "get_ga4gh_variantannotationsets",
    "Description": "Return meta data for a specific annotation set in GA4GH format by ID",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 339,
    "Tool Name": "get_genetree_member_id",
    "Description": "Retrieve the gene tree containing a gene identified by its Ensembl ID.Find the phylogenetic tree showing evolutionary relationships for a gene of interest,\nidentified using its Ensembl stable identifier. Useful for understanding gene evolution\nwhen you have the specific gene ID.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 340,
    "Tool Name": "get_info_biotypes_groups",
    "Description": "With :group argument provided, list the properties of biotypes within that group. Object type (gene or transcript) can be provided for filtering.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 341,
    "Tool Name": "get_info_biotypes_name",
    "Description": "List the properties of biotypes with a given name. Object type (gene or transcript) can be provided for filtering.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 342,
    "Tool Name": "get_info_compara_species_sets",
    "Description": "List all collections of species analysed with the specified compara method.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 343,
    "Tool Name": "get_info_comparas",
    "Description": "Get all available comparative genomics databases and their data release.",
    "category": "Databases",
    "Server Name": "Origene-Ensembl"
  },
  {
    "IDX": 344,
    "Tool Name": "list_genomes",
    "Description": "Get all supported genome assemblies from UCSC Genome Browser.",
    "category": "Databases",
    "Server Name": "Origene-UCSC"
  },
  {
    "IDX": 345,
    "Tool Name": "list_tracks",
    "Description": "List all tracks for a specific genome assembly.",
    "category": "Databases",
    "Server Name": "Origene-UCSC"
  },
  {
    "IDX": 346,
    "Tool Name": "list_hub_tracks",
    "Description": "List all tracks in a specific track hub for a genome.",
    "category": "Databases",
    "Server Name": "Origene-UCSC"
  },
  {
    "IDX": 347,
    "Tool Name": "list_chromosomes",
    "Description": "List all chromosomes for a genome assembly.",
    "category": "Databases",
    "Server Name": "Origene-UCSC"
  },
  {
    "IDX": 348,
    "Tool Name": "list_public_hubs",
    "Description": "Get list of all public UCSC track hubs.",
    "category": "Databases",
    "Server Name": "Origene-UCSC"
  },
  {
    "IDX": 349,
    "Tool Name": "get_chromosome_sequence",
    "Description": "Get sequence for an entire chromosome.",
    "category": "Databases",
    "Server Name": "Origene-UCSC"
  },
  {
    "IDX": 350,
    "Tool Name": "get_sequence",
    "Description": "Get DNA sequence for a genomic region.",
    "category": "Databases",
    "Server Name": "Origene-UCSC"
  },
  {
    "IDX": 351,
    "Tool Name": "get_track_data",
    "Description": "Get data from a specific track for a genomic region.",
    "category": "Computational Tools",
    "Server Name": "Origene-UCSC"
  },
  {
    "IDX": 352,
    "Tool Name": "get_cytoband",
    "Description": "Get cytoband (chromosome banding) information for a specified genome and chromosome.",
    "category": "Databases",
    "Server Name": "Origene-UCSC"
  },
  {
    "IDX": 353,
    "Tool Name": "get_active_ingredient_info_by_drug_name",
    "Description": "Fetch a list of active ingredients in a specific drug product.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 354,
    "Tool Name": "get_dosage_and_storage_information_by_drug_name",
    "Description": "Retrieve dosage and storage information for a specific drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 355,
    "Tool Name": "get_drug_names_by_abuse_info",
    "Description": "Retrieve drug names based on information about types of abuse and adverse reactions pertinent to those types of abuse.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 356,
    "Tool Name": "get_abuse_info_by_drug_name",
    "Description": "Retrieve information about types of abuse based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 357,
    "Tool Name": "get_drug_names_by_accessories",
    "Description": "Retrieve drug names based on the accessories field information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 358,
    "Tool Name": "get_accessories_info_by_drug_name",
    "Description": "Retrieve information about accessories based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 359,
    "Tool Name": "get_drug_names_by_active_ingredient",
    "Description": "Retrieve drug names based on the active ingredient information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 360,
    "Tool Name": "get_active_ingredient_application_number_manufacturer_name_NDC_number_administration_route_by_drug_name",
    "Description": "Retrieve detailed information about a drug's active ingredient, FDA application number, manufacturer name, National Drug Code (NDC) number, and route of administration; all based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 361,
    "Tool Name": "get_drug_names_by_application_number_manufacturer_name_NDC_number",
    "Description": "Retrieve drug names based on the specified FDA application number, manufacturer name, or National Drug Code (NDC) number.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 362,
    "Tool Name": "get_drug_name_by_adverse_reaction",
    "Description": "Retrieve the drug name based on specific adverse reactions reported.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 363,
    "Tool Name": "get_adverse_reactions_by_drug_name",
    "Description": "Retrieve adverse reactions information based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 364,
    "Tool Name": "get_drug_names_by_alarm",
    "Description": "Retrieve drug names based on the presence of specific alarms, which are related to adverse reaction events.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 365,
    "Tool Name": "get_alarms_by_drug_name",
    "Description": "Retrieve alarms based on the specified drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 366,
    "Tool Name": "get_drug_names_by_animal_pharmacology_info",
    "Description": "Retrieve drug names based on animal pharmacology and toxicology information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 367,
    "Tool Name": "get_animal_pharmacology_info_by_drug_name",
    "Description": "Retrieve animal pharmacology and toxicology information based on drug names.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 368,
    "Tool Name": "get_drug_name_by_info_on_conditions_for_doctor_consultation",
    "Description": "Retrieve the drug names that require asking a doctor before use due to a patient's specific conditions and symptoms.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 369,
    "Tool Name": "get_info_on_conditions_for_doctor_consultation_by_drug_name",
    "Description": "Get information about when a doctor should be consulted before using a specific drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 370,
    "Tool Name": "get_drug_names_by_info_on_consulting_doctor_pharmacist_for_drug_interactions",
    "Description": "Retrieve drug names based on information about when a doctor or pharmacist should be consulted regarding drug interactions.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 371,
    "Tool Name": "get_info_on_consulting_doctor_pharmacist_for_drug_interactions_by_drug_name",
    "Description": "Get information about when a doctor or pharmacist should be consulted regarding drug interactions for a specific drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 372,
    "Tool Name": "get_drug_names_by_assembly_installation_info",
    "Description": "Retrieve drug names based on assembly or installation instructions.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 373,
    "Tool Name": "get_assembly_installation_info_by_drug_name",
    "Description": "Retrieve assembly or installation instructions based on drug names.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 374,
    "Tool Name": "get_drug_names_by_boxed_warning",
    "Description": "Retrieve drug names that have specific boxed warnings (The most serious risk alerts in drug labeling) and adverse effects.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 375,
    "Tool Name": "get_boxed_warning_info_by_drug_name",
    "Description": "Retrieve boxed warning (The most serious risk alerts in drug labeling) and adverse effects information for a specific drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 376,
    "Tool Name": "get_drug_name_by_calibration_instructions",
    "Description": "Retrieve the drug name based on the calibration instructions provided.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 377,
    "Tool Name": "get_calibration_instructions_by_drug_name",
    "Description": "Retrieve calibration instructions based on the specified drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 378,
    "Tool Name": "get_drug_names_by_carcinogenic_mutagenic_fertility_impairment_info",
    "Description": "Retrieve drug names based on the presence of carcinogenic, mutagenic, or fertility impairment information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 379,
    "Tool Name": "get_carcinogenic_mutagenic_fertility_impairment_info_by_drug_name",
    "Description": "Retrieve carcinogenic, mutagenic, or fertility impairment information based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 380,
    "Tool Name": "get_drug_name_by_application_number_NUI_identifier_SPL_document_ID_SPL_set_ID",
    "Description": "Retrieves the drug name based on various identifiers such as the FDA application number, NUI, SPL document ID, or SPL set ID. The tool takes a single string containing the identifier.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 381,
    "Tool Name": "get_drug_names_by_clinical_pharmacology",
    "Description": "Retrieves drug names based on a search query within their clinical pharmacology information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 382,
    "Tool Name": "get_clinical_pharmacology_by_drug_name",
    "Description": "Retrieves clinical pharmacology information for a specific drug from the FDA database by its name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 383,
    "Tool Name": "get_drug_names_by_clinical_studies",
    "Description": "Retrieves drug names based on a search query within their clinical studies information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 384,
    "Tool Name": "get_clinical_studies_info_by_drug_name",
    "Description": "Retrieves clinical studies information for a specific drug from the FDA database by its name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 385,
    "Tool Name": "get_drug_names_by_contraindications",
    "Description": "Retrieves drug names based on a search query within their contraindications information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 386,
    "Tool Name": "get_contraindications_by_drug_name",
    "Description": "Retrieve contraindications information based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 387,
    "Tool Name": "get_drug_names_by_controlled_substance_DEA_schedule",
    "Description": "Retrieves drug names based on a search query within their DEA controlled substance schedule information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 388,
    "Tool Name": "get_controlled_substance_DEA_schedule_info_by_drug_name",
    "Description": "Retrieves information about the controlled substance DEA (Drug Enforcement Administration) schedule for a specific drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 389,
    "Tool Name": "get_drug_name_by_dependence_info",
    "Description": "Retrieve the drug name based on information about dependence characteristics.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 390,
    "Tool Name": "get_dependence_info_by_drug_name",
    "Description": "Retrieves information about dependence characteristics for a specific drug from the FDA database by its name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 391,
    "Tool Name": "get_drug_names_by_disposal_info",
    "Description": "Retrieves drug names based on a search query within their disposal and waste handling information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 392,
    "Tool Name": "get_disposal_info_by_drug_name",
    "Description": "Retrieves disposal and waste handling information for a specific drug from the FDA database by its name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 393,
    "Tool Name": "get_drug_name_by_dosage_info",
    "Description": "Retrieve the drug name based on dosage and administration information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 394,
    "Tool Name": "get_drug_names_by_dosage_forms_and_strengths_info",
    "Description": "Retrieves drug names based on a search query within their dosage forms and strengths information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 395,
    "Tool Name": "get_dosage_forms_and_strengths_by_drug_name",
    "Description": "Retrieves dosage forms and strengths information for a specific drug from the FDA database by its name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 396,
    "Tool Name": "get_drug_name_by_abuse_types_and_related_adverse_reactions_and_controlled_substance_status",
    "Description": "Retrieves the drug name based on information about drug abuse and dependence, including whether the drug is a controlled substance, the types of possible abuse, and adverse reactions relevant to those abuse types.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 397,
    "Tool Name": "get_abuse_types_and_related_adverse_reactions_and_controlled_substance_status_by_drug_name",
    "Description": "Retrieves information about drug abuse and dependence for a specific drug from the FDA database by its name, including controlled substance status, abuse types, and related adverse reactions.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 398,
    "Tool Name": "get_drug_names_by_lab_test_interference",
    "Description": "Retrieves drug names based on a search query within their laboratory test interference information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 399,
    "Tool Name": "get_lab_test_interference_info_by_drug_name",
    "Description": "Retrieves information about laboratory test interferences for a specific drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 400,
    "Tool Name": "get_drug_names_by_drug_interactions",
    "Description": "Retrieves a list of drug names that have known interactions with a specified term.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 401,
    "Tool Name": "get_drug_interactions_by_drug_name",
    "Description": "Retrieve drug interactions based on the specified drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 402,
    "Tool Name": "get_drug_names_by_effective_time",
    "Description": "Retrieve drug names based on the effective time of the labeling document.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 403,
    "Tool Name": "get_effective_time_by_drug_name",
    "Description": "Retrieve effective time of the labeling document based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 404,
    "Tool Name": "get_drug_name_by_environmental_warning",
    "Description": "Retrieve the drug name based on the specified environmental warnings.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 405,
    "Tool Name": "get_environmental_warning_by_drug_name",
    "Description": "Fetch environmental warnings for a specific drug based on its name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 406,
    "Tool Name": "get_drug_names_by_food_safety_warnings",
    "Description": "Retrieve drug names based on specific food safety warnings.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 407,
    "Tool Name": "get_drug_names_by_general_precautions",
    "Description": "Retrieve drug names based on specific general precautions information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 408,
    "Tool Name": "get_general_precautions_by_drug_name",
    "Description": "Retrieve general precautions information based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 409,
    "Tool Name": "get_drug_names_by_geriatric_use",
    "Description": "Retrieve drug names that have specific information about geriatric use.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 410,
    "Tool Name": "get_geriatric_use_info_by_drug_name",
    "Description": "Retrieve information about geriatric use based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 411,
    "Tool Name": "get_dear_health_care_provider_letter_info_by_drug_name",
    "Description": "Fetch information about dear health care provider letters for a specific drug. The letters are sent by drug manufacturers to provide new or updated information about the drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 412,
    "Tool Name": "get_drug_names_by_dear_health_care_provider_letter_info",
    "Description": "Fetch drug names based on information about dear health care provider letters. The letters are sent by drug manufacturers to provide new or updated information about the drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 413,
    "Tool Name": "get_drug_names_by_health_claim",
    "Description": "Retrieve drug names based on specific health claims.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 414,
    "Tool Name": "get_health_claims_by_drug_name",
    "Description": "Retrieve health claims associated with a specific drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 415,
    "Tool Name": "get_drug_name_by_document_id",
    "Description": "Retrieve the drug name based on the document ID.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 416,
    "Tool Name": "get_document_id_by_drug_name",
    "Description": "Retrieve the document ID based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 417,
    "Tool Name": "get_drug_name_by_inactive_ingredient",
    "Description": "Retrieve the drug name based on the inactive ingredient information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 418,
    "Tool Name": "get_inactive_ingredient_info_by_drug_name",
    "Description": "Fetch a list of inactive ingredients in a specific drug product based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 419,
    "Tool Name": "get_drug_names_by_indication",
    "Description": "Retrieve a list of drug names based on a specific indication or usage.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 420,
    "Tool Name": "get_indications_by_drug_name",
    "Description": "Retrieve indications and usage information based on a specific drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 421,
    "Tool Name": "get_drug_names_by_information_for_owners_or_caregivers",
    "Description": "Retrieve drug names based on information for owners or caregivers.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 422,
    "Tool Name": "get_information_for_owners_or_caregivers_by_drug_name",
    "Description": "Retrieves information for owners or caregivers for a specific drug. Note: The underlying data is sparse, and many drugs may not have this information, resulting in a null return.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 423,
    "Tool Name": "get_info_for_patients_by_drug_name",
    "Description": "Retrieves patient information for a specific drug from the FDA database by its name. Note: The underlying data is sparse, and many drugs may not have this information, resulting in a null return.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 424,
    "Tool Name": "get_drug_names_by_instructions_for_use",
    "Description": "Retrieves drug names based on a search query within their instructions for use.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 425,
    "Tool Name": "get_instructions_for_use_by_drug_name",
    "Description": "Retrieves instructions for use information for a specific drug from the FDA database by its name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 426,
    "Tool Name": "retrieve_drug_name_by_device_use",
    "Description": "Retrieves drug names based on a search query within the intended use of their associated medical device.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 427,
    "Tool Name": "retrieve_device_use_by_drug_name",
    "Description": "Retrieves the intended use of the device associated with a specific drug. Note: The underlying data is sparse, and many drugs may not have this information, resulting in a null return.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 428,
    "Tool Name": "get_drug_names_by_child_safety_info",
    "Description": "Retrieves drug names based on a search query within their child safety information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 429,
    "Tool Name": "get_child_safety_info_by_drug_name",
    "Description": "Retrieves child safety information for a specific drug based on its name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 430,
    "Tool Name": "get_drug_name_by_labor_and_delivery_info",
    "Description": "Retrieve the drug name based on information about the drug¡¯s use during labor or delivery.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 431,
    "Tool Name": "get_labor_and_delivery_info_by_drug_name",
    "Description": "Retrieves information about the drug¡¯s use during labor or delivery based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 432,
    "Tool Name": "get_drug_names_by_lab_tests",
    "Description": "Retrieves drug names based on a search query within their laboratory tests information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 433,
    "Tool Name": "get_lab_tests_by_drug_name",
    "Description": "Retrieves laboratory tests information for a specific drug from the FDA database by its name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 434,
    "Tool Name": "get_mechanism_of_action_by_drug_name",
    "Description": "Retrieves the mechanism of action information for a specific drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 435,
    "Tool Name": "get_drug_names_by_mechanism_of_action",
    "Description": "Retrieves drug names based on a search query within their mechanism of action information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 436,
    "Tool Name": "get_drug_name_by_microbiology",
    "Description": "Retrieve the drug name based on microbiology field information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 437,
    "Tool Name": "get_microbiology_info_by_drug_name",
    "Description": "Retrieves microbiology information for a specific drug from the FDA database by its name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 438,
    "Tool Name": "get_drug_names_by_nonclinical_toxicology_info",
    "Description": "Retrieves drug names based on a search query within their nonclinical toxicology information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 439,
    "Tool Name": "get_nonclinical_toxicology_info_by_drug_name",
    "Description": "Retrieves nonclinical toxicology information for a specific drug from the FDA database by its name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 440,
    "Tool Name": "get_drug_names_by_nonteratogenic_effects",
    "Description": "Retrieve drug names based on the presence of nonteratogenic effects information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 441,
    "Tool Name": "get_nonteratogenic_effects_by_drug_name",
    "Description": "Retrieves information about nonteratogenic effects for a specific drug from the FDA database by its name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 442,
    "Tool Name": "get_drug_names_by_info_for_nursing_mothers",
    "Description": "Retrieves drug names based on a search query within their information for nursing mothers.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 443,
    "Tool Name": "get_info_for_nursing_mothers_by_drug_name",
    "Description": "Retrieves information about nursing mothers for a specific drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 444,
    "Tool Name": "get_drug_name_by_other_safety_info",
    "Description": "Retrieves drug names based on a search query within their 'other safety information' section.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 445,
    "Tool Name": "get_other_safety_info_by_drug_name",
    "Description": "Retrieves other safety information for a specific drug from the FDA database by its name. This information may not be specified in other fields.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 446,
    "Tool Name": "get_drug_names_by_overdosage_info",
    "Description": "Retrieves drug names based on a search query within their overdosage information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 447,
    "Tool Name": "get_overdosage_info_by_drug_name",
    "Description": "Retrieves information about signs, symptoms, and laboratory findings of acute overdosage based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 448,
    "Tool Name": "get_drug_name_by_principal_display_panel",
    "Description": "Retrieve the drug name based on the content of the principal display panel of the product package.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 449,
    "Tool Name": "get_principal_display_panel_by_drug_name",
    "Description": "Retrieve the content of the principal display panel of the product package based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 450,
    "Tool Name": "retrieve_drug_names_by_patient_medication_info",
    "Description": "Retrieve drug names based on patient medication information, which is about safe use of the drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 451,
    "Tool Name": "retrieve_patient_medication_info_by_drug_name",
    "Description": "Retrieve patient medication information (which is about safe use of the drug) based on drug names.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 452,
    "Tool Name": "get_drug_names_by_pediatric_use",
    "Description": "Retrieve drug names based on pediatric use information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 453,
    "Tool Name": "get_pediatric_use_info_by_drug_name",
    "Description": "Retrieve pediatric use information based on drug names.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 454,
    "Tool Name": "get_drug_name_by_pharmacodynamics",
    "Description": "Retrieve the drug name based on pharmacodynamics information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 455,
    "Tool Name": "get_pharmacodynamics_by_drug_name",
    "Description": "Retrieve pharmacodynamics information based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 456,
    "Tool Name": "get_pharmacogenomics_info_by_drug_name",
    "Description": "Retrieve pharmacogenomics information based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 457,
    "Tool Name": "get_drug_names_by_pharmacokinetics",
    "Description": "Retrieve drug names based on specific pharmacokinetics information, such as absorption, distribution, elimination, metabolism, drug interactions, and specific patient populations.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 458,
    "Tool Name": "get_pharmacokinetics_by_drug_name",
    "Description": "Retrieve pharmacokinetics information (e.g. absorption, distribution, elimination, metabolism, drug interactions, and specific patient populations) for a specific drug based on its name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 459,
    "Tool Name": "get_drug_name_by_precautions",
    "Description": "Retrieve the drug name based on the precautions field information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 460,
    "Tool Name": "get_precautions_by_drug_name",
    "Description": "Retrieve precautions information based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 461,
    "Tool Name": "get_drug_names_by_pregnancy_effects_info",
    "Description": "Retrieves drug names based on a search query within their pregnancy effects information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 462,
    "Tool Name": "get_pregnancy_effects_info_by_drug_name",
    "Description": "Retrieves information about the effects on pregnancy for a specific drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 463,
    "Tool Name": "get_drug_name_by_pregnancy_or_breastfeeding_info",
    "Description": "Retrieve the drug names based on pregnancy or breastfeeding information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 464,
    "Tool Name": "get_pregnancy_or_breastfeeding_info_by_drug_name",
    "Description": "Retrieves pregnancy or breastfeeding information for a specific drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 465,
    "Tool Name": "get_contact_for_questions_info_by_drug_name",
    "Description": "Retrieves information on who to contact with questions about the drug based on the provided drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 466,
    "Tool Name": "get_recent_changes_by_drug_name",
    "Description": "Retrieves recent major changes in labeling for a specific drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 467,
    "Tool Name": "get_drug_name_by_reference",
    "Description": "Retrieves the drug name based on the reference information provided in the drug labeling. Note: The underlying data is sparse, and many search terms may result in a null return.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 468,
    "Tool Name": "get_reference_info_by_drug_name",
    "Description": "Retrieves reference information for a specific drug from the FDA database by its name. Note: The underlying data is sparse, and many drugs may not have this information, resulting in a null return.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 469,
    "Tool Name": "get_drug_names_by_residue_warning",
    "Description": "Retrieves drug names based on a search query within their residue warning information. Note: The underlying data is sparse, and many search terms may result in a null return.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 470,
    "Tool Name": "get_residue_warning_by_drug_name",
    "Description": "Retrieves the residue warning for a specific drug from the FDA database by its name. Note: The underlying data is sparse, and many drugs may not have this information, resulting in a null return.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 471,
    "Tool Name": "get_drug_names_by_risk",
    "Description": "Retrieves drug names based on a search query within their risk information, especially regarding pregnancy or breastfeeding.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 472,
    "Tool Name": "get_risk_info_by_drug_name",
    "Description": "Retrieves risk information (especially regarding pregnancy or breastfeeding) for a specific drug from the FDA database by its name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 473,
    "Tool Name": "get_drug_names_by_route",
    "Description": "Retrieves drug names based on a search query within their route of administration.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 474,
    "Tool Name": "get_route_info_by_drug_name",
    "Description": "Retrieve the route of administration information based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 475,
    "Tool Name": "get_drug_names_by_safe_handling_warning",
    "Description": "Retrieves drug names based on a search query within their safe handling warning information. Note: The underlying data is sparse, and many search terms may result in a null return.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 476,
    "Tool Name": "get_safe_handling_warnings_by_drug_name",
    "Description": "Retrieves safe handling warnings for a specific drug from the FDA database by its name. Note: The underlying data is sparse, and many drugs may not have this information, resulting in a null return.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 477,
    "Tool Name": "get_drug_name_by_set_id",
    "Description": "Retrieve the drug name based on the Set ID of the labeling.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 478,
    "Tool Name": "get_drug_names_by_spl_indexing_data_elements",
    "Description": "Retrieve drug names based on Structured Product Labeling (SPL) indexing data elements.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 479,
    "Tool Name": "get_spl_indexing_data_elements_by_drug_name",
    "Description": "Retrieve Structured Product Labeling (SPL) indexing data elements based on drug names.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 480,
    "Tool Name": "get_drug_names_by_medication_guide",
    "Description": "Retrieve drug names based on the presence of specific information in the medication guide.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 481,
    "Tool Name": "get_medication_guide_info_by_drug_name",
    "Description": "Retrieve medication guide information based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 482,
    "Tool Name": "get_drug_name_from_patient_package_insert",
    "Description": "Retrieve the drug name based on the information provided in the patient package insert.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 483,
    "Tool Name": "get_patient_package_insert_from_drug_name",
    "Description": "Retrieve the patient package insert information based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 484,
    "Tool Name": "get_drug_names_by_ingredient",
    "Description": "Retrieve drug names based on a specific ingredient present in the drug product.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 485,
    "Tool Name": "get_ingredients_by_drug_name",
    "Description": "Retrieve a list of drug ingredients based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 486,
    "Tool Name": "get_spl_unclassified_section_by_drug_name",
    "Description": "Retrieve the SPL unclassified section information (Content not yet clearly categorized) based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 487,
    "Tool Name": "get_drug_name_by_stop_use_info",
    "Description": "Retrieve the drug name based on the stop use information provided.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 488,
    "Tool Name": "get_stop_use_info_by_drug_name",
    "Description": "Retrieve stop use information based on the drug name provided.\n\nQuery example: {\"drug_name\": \"aspirin\", \"limit\": 1, \"skip\": 0}\n\nReturns:\n    A dictionary containing the following fields:\n        meta:Pagination metadata,containing following fields:\n            skip:Current number of skipped records\n            limit:Number of records returned this time\n            total:Total number of records matched\n        results:Actual data array,containing the following fields:\n            openfda.brand_name:Drug brand name\n            openfda.generic_name:Drug generic name\n            stop_use:If you experience any of these symptoms, stop using this medicine.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 489,
    "Tool Name": "get_drug_name_by_storage_and_handling_info",
    "Description": "Retrieve the drug name based on storage and handling information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 490,
    "Tool Name": "get_storage_and_handling_info_by_drug_name",
    "Description": "Retrieve storage and handling information based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 491,
    "Tool Name": "get_drug_names_by_safety_summary",
    "Description": "Retrieve drug names based on the summary of safety and effectiveness information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 492,
    "Tool Name": "get_safety_summary_by_drug_name",
    "Description": "Retrieve a summary of safety and effectiveness information based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 493,
    "Tool Name": "get_drug_names_by_teratogenic_effects",
    "Description": "Retrieve drug names based on specific teratogenic effects categories.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 494,
    "Tool Name": "get_teratogenic_effects_by_drug_name",
    "Description": "Retrieve teratogenic effects information based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 495,
    "Tool Name": "get_drug_names_by_population_use",
    "Description": "Retrieve drug names based on their use in specific populations, such as pregnant women, nursing mothers, pediatric patients, and geriatric patients.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 496,
    "Tool Name": "get_population_use_info_by_drug_name",
    "Description": "Retrieve information about the use of a drug in specific populations based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 497,
    "Tool Name": "get_user_safety_warning_by_drug_names",
    "Description": "Retrieve specific user safety warnings based on drug names.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 498,
    "Tool Name": "get_drug_names_by_user_safety_warning",
    "Description": "Retrieve drug names that have specific user safety warnings.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 499,
    "Tool Name": "get_drug_name_by_warnings",
    "Description": "Retrieve the drug names based on specific warning information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 500,
    "Tool Name": "get_warnings_by_drug_name",
    "Description": "Retrieve warning information based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 501,
    "Tool Name": "get_warnings_and_cautions_by_drug_name",
    "Description": "Retrieve warnings and cautions information for a specific drug based on its name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 502,
    "Tool Name": "get_drug_names_by_warnings_and_cautions",
    "Description": "Retrieve drug names based on specific warnings and cautions information.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 503,
    "Tool Name": "get_when_using_info",
    "Description": "Retrieve information about side effects and substances or activities to avoid while using a specific drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 504,
    "Tool Name": "get_brand_name_generic_name",
    "Description": "Retrieve the brand name and generic name from generic name or brand name of a drug.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 505,
    "Tool Name": "get_do_not_use_info_by_drug_name",
    "Description": "Retrieve information about all contraindications for use based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 506,
    "Tool Name": "get_purpose_info_by_drug_name",
    "Description": "Retrieve about the drug product¡¯s indications for use based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 507,
    "Tool Name": "get_drug_generic_name",
    "Description": "Get the drug¡¯s generic name based on the drug's generic or brand name.",
    "category": "Databases",
    "Server Name": "Origene-FDADrug"
  },
  {
    "IDX": 508,
    "Tool Name": "get_associated_targets_by_disease_efoId",
    "Description": "Find targets associated with a specific disease or phenotype based on efoId.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 509,
    "Tool Name": "get_associated_diseases_phenotypes_by_target_ensemblID",
    "Description": "Find diseases or phenotypes associated with a specific target using ensemblId.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 510,
    "Tool Name": "target_disease_evidence",
    "Description": "Explore evidence that supports a specific target-disease association.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 511,
    "Tool Name": "get_drug_warnings_by_chemblId",
    "Description": "Retrieve warnings for a specific drug using ChEMBL ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 512,
    "Tool Name": "get_drug_mechanisms_of_action_by_chemblId",
    "Description": "Retrieve the mechanisms of action associated with a specific drug using chemblId.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 513,
    "Tool Name": "get_associated_drugs_by_disease_efoId",
    "Description": "Retrieve known drugs associated with a specific disease by disease efoId.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 514,
    "Tool Name": "get_similar_entities_by_disease_efoId",
    "Description": "Retrieve similar entities for a given disease efoId using a model trained with PubMed.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 515,
    "Tool Name": "get_similar_entities_by_drug_chemblId",
    "Description": "Retrieve similar entities for a given drug chemblId using a model trained with PubMed.",
    "category": "Model Services",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 516,
    "Tool Name": "get_similar_entities_by_target_ensemblID",
    "Description": "Retrieve similar entities for a given target ensemblID using a model trained with PubMed.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 517,
    "Tool Name": "get_associated_phenotypes_by_disease_efoId",
    "Description": "Find HPO phenotypes associated with the specified disease efoId.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 518,
    "Tool Name": "get_drug_withdrawn_blackbox_status_by_chemblId",
    "Description": "Find withdrawn and black-box warning statuses for a specific drug by chemblId.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 519,
    "Tool Name": "search_category_counts_by_query_string",
    "Description": "Get the count of entries in each entity category (disease, target, drug) based on a query string.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 520,
    "Tool Name": "get_disease_id_description_by_name",
    "Description": "Retrieve the efoId and additional details of a disease based on its name.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 521,
    "Tool Name": "get_drug_id_description_by_name",
    "Description": "Fetch the drug chemblId and description based on the drug name.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 522,
    "Tool Name": "get_drug_indications_by_chemblId",
    "Description": "Fetch indications (treatable phenotypes/diseases) for a given drug chemblId.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 523,
    "Tool Name": "get_target_gene_ontology_by_ensemblID",
    "Description": "Retrieve Gene Ontology annotations for a specific target by Ensembl ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 524,
    "Tool Name": "get_target_homologues_by_ensemblID",
    "Description": "Fetch homologues for a specific target by Ensembl ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 525,
    "Tool Name": "get_target_safety_profile_by_ensemblID",
    "Description": "Retrieve known target safety liabilities for a specific target Ensembl ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 526,
    "Tool Name": "get_biological_mouse_models_by_ensemblID",
    "Description": "Retrieve biological mouse models, including allelic compositions and genetic backgrounds, for a specific target.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 527,
    "Tool Name": "get_target_genomic_location_by_ensemblID",
    "Description": "Retrieve genomic location data for a specific target, including chromosome, start, end, and strand.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 528,
    "Tool Name": "get_target_subcellular_locations_by_ensemblID",
    "Description": "Retrieve information about subcellular locations for a specific target Ensembl ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 529,
    "Tool Name": "get_target_synonyms_by_ensemblID",
    "Description": "Retrieve synonyms for a specified target, including alternative names and symbols, using the given Ensembl ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 530,
    "Tool Name": "get_target_tractability_by_ensemblID",
    "Description": "Retrieve tractability assessments, including modality and values, for a specific target Ensembl ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 531,
    "Tool Name": "get_target_classes_by_ensemblID",
    "Description": "Retrieve the target classes associated with a specific target Ensembl ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 532,
    "Tool Name": "get_target_enabling_packages_by_ensemblID",
    "Description": "Retrieve the Target Enabling Packages (TEP) associated with a specific target Ensembl ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 533,
    "Tool Name": "get_target_interactions_by_ensemblID",
    "Description": "Retrieve interaction data for a specific target Ensembl ID, including interaction partners and evidence.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 534,
    "Tool Name": "get_disease_ancestors_parents_by_efoId",
    "Description": "Retrieve the disease ancestors and parents in the ontology using the disease EFO ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 535,
    "Tool Name": "get_disease_descendants_children_by_efoId",
    "Description": "Retrieve the disease descendants and children in the ontology using the disease EFO ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 536,
    "Tool Name": "get_disease_locations_by_efoId",
    "Description": "Retrieve the disease's direct location and indirect location disease terms and IDs using the disease EFO ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 537,
    "Tool Name": "get_disease_synonyms_by_efoId",
    "Description": "Retrieve disease synonyms by its EFO ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 538,
    "Tool Name": "get_disease_description_by_efoId",
    "Description": "Retrieve disease description, name, database cross-references, obsolete terms, and whether it's a therapeutic area, all using the specified EFO ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 539,
    "Tool Name": "get_disease_therapeutic_areas_by_efoId",
    "Description": "Retrieve the therapeutic areas associated with a specific disease EFO ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 540,
    "Tool Name": "get_drug_adverse_events_by_chemblId",
    "Description": "Retrieve significant adverse events reported for a specific drug ChEMBL ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 541,
    "Tool Name": "get_known_drugs_by_drug_chemblId",
    "Description": "Get a list of known drugs and associated information using the specified ChEMBL ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 542,
    "Tool Name": "get_parent_child_molecules_by_drug_chembl_ID",
    "Description": "Get parent and child molecules of specified drug ChEMBL ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 543,
    "Tool Name": "get_approved_indications_by_drug_chemblId",
    "Description": "Retrieve detailed information about multiple drugs using a list of ChEMBL IDs.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 544,
    "Tool Name": "get_drug_description_by_chemblId",
    "Description": "Get drug name, year of first approval, type, cross references, and max clinical trial phase based on specified chemblId.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 545,
    "Tool Name": "get_drug_synonyms_by_chemblId",
    "Description": "Retrieve the synonyms associated with a specific drug ChEMBL ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 546,
    "Tool Name": "get_drug_trade_names_by_chemblId",
    "Description": "Retrieve the trade names associated with a specific drug ChEMBL ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 547,
    "Tool Name": "get_drug_approval_status_by_chemblId",
    "Description": "Retrieve the approval status of a specific drug ChEMBL ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 548,
    "Tool Name": "get_chemical_probes_by_target_ensemblID",
    "Description": "Retrieve chemical probes associated with a specific target using its Ensembl ID.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 549,
    "Tool Name": "drug_pharmacogenomics_data",
    "Description": "Retrieve pharmacogenomics data for a specific drug, including evidence levels and genotype annotations.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 550,
    "Tool Name": "get_associated_drugs_by_target_ensemblID",
    "Description": "Get known drugs associated with a specific target Ensembl ID, including clinical trial phase and mechanism of action of the drugs.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 551,
    "Tool Name": "get_associated_diseases_by_drug_chemblId",
    "Description": "Retrieve the list of diseases associated with a specific drug ChEMBL ID based on clinical trial data or post-marketed drugs.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 552,
    "Tool Name": "get_associated_targets_by_drug_chemblId",
    "Description": "Retrieve the list of targets linked to a specific drug ChEMBL ID based on its mechanism of action.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 553,
    "Tool Name": "multi_entity_search_by_query_string",
    "Description": "Perform a multi-entity search based on a query string, filtering by entity names and pagination settings.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 554,
    "Tool Name": "get_gene_ontology_terms_by_goID",
    "Description": "Retrieve Gene Ontology terms based on a list of GO IDs.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 555,
    "Tool Name": "get_target_constraint_info_by_ensemblID",
    "Description": "Retrieve genetic constraint information for a specific target Ensembl ID, including expected and observed values, and scores.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 556,
    "Tool Name": "get_publications_by_disease_efoId",
    "Description": "Retrieve publications related to a disease EFO ID, including PubMed IDs and publication dates.",
    "category": "Literature Search",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 557,
    "Tool Name": "get_publications_by_target_ensemblID",
    "Description": "Retrieve publications related to a target Ensembl ID, including PubMed IDs and publication dates.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 558,
    "Tool Name": "get_publications_by_drug_chemblId",
    "Description": "Retrieve publications related to a drug ChEMBL ID, including PubMed IDs and publication dates.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 559,
    "Tool Name": "get_target_id_description_by_name",
    "Description": "Get the Ensembl ID and description based on the target name.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 560,
    "Tool Name": "get_general_info_by_disease_name",
    "Description": "Get disease EFO ID and description by disease name from OpenTargets.\nDescription information will include disease name, EFO ID, disease targets, related drugs and related disease phenotypes.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 561,
    "Tool Name": "get_target_ensembl_id",
    "Description": "Get target Ensembl ID by target name.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 562,
    "Tool Name": "get_disease_efo_id",
    "Description": "Get disease EFO ID by disease name.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 563,
    "Tool Name": "get_drug_chembl_id_by_name",
    "Description": "Find drug ChEMBL ID by drug name.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 564,
    "Tool Name": "get_associated_targets_by_disease_name",
    "Description": "Find targets associated with a specific disease or phenotype based on its name.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 565,
    "Tool Name": "get_associated_diseases_phenotypes_by_target_name",
    "Description": "Find diseases or phenotypes associated with a specific target.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 566,
    "Tool Name": "get_target_disease_evidence_by_name",
    "Description": "Explore evidence that supports a specific target-disease association. Input is disease name and target name.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 567,
    "Tool Name": "get_drug_warnings_by_name",
    "Description": "Retrieve warnings for a specific drug.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 568,
    "Tool Name": "get_drug_mechanisms_of_action_by_name",
    "Description": "Retrieve the mechanisms of action associated with a specific drug.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 569,
    "Tool Name": "get_associated_drugs_by_disease_name",
    "Description": "Retrieve known drugs associated with a specific disease by disease name.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 570,
    "Tool Name": "get_similar_entities_by_disease_name",
    "Description": "Retrieve similar entities for a given disease using a model trained with PubMed.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 571,
    "Tool Name": "get_similar_entities_by_drug_name",
    "Description": "Retrieve similar entities for a given drug using a model trained with PubMed.",
    "category": "Model Services",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 572,
    "Tool Name": "get_similar_entities_by_target_name",
    "Description": "Retrieve similar entities for a given target using a model trained with PubMed.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 573,
    "Tool Name": "get_associated_phenotypes_by_disease_name",
    "Description": "Find HPO phenotypes asosciated with the specified disease.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 574,
    "Tool Name": "get_drug_indications_by_name",
    "Description": "Fetch indications (treatable phenotypes/diseases) for a given drug.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 575,
    "Tool Name": "get_target_gene_ontology_by_name",
    "Description": "Retrieve Gene Ontology annotations for a specific target.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 576,
    "Tool Name": "get_target_homologues_by_name",
    "Description": "Fetch homologues for a specific target.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 577,
    "Tool Name": "get_target_safety_profile_by_name",
    "Description": "Retrieve known target safety liabilities for a specific target.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 578,
    "Tool Name": "get_biological_mouse_models_by_target_name",
    "Description": "Retrieve biological mouse models, including allelic compositions and genetic backgrounds, for a specific target.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 579,
    "Tool Name": "get_target_genomic_location_by_name",
    "Description": "Retrieve genomic location data for a specific target, including chromosome, start, end, and strand.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 580,
    "Tool Name": "get_target_subcellular_locations_by_name",
    "Description": "Retrieve information about subcellular locations for a specific target.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 581,
    "Tool Name": "get_target_synonyms_by_name",
    "Description": "Retrieve synonyms for specified target, including alternative names and symbols.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 582,
    "Tool Name": "get_target_tractability_by_name",
    "Description": "Retrieve tractability assessments, including modality and values.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 583,
    "Tool Name": "get_target_classes_by_name",
    "Description": "Retrieve the target classes associated with a specific target.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 584,
    "Tool Name": "get_target_enabling_packages_by_name",
    "Description": "Retrieve the Target Enabling Packages (TEP) associated with a specific target.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 585,
    "Tool Name": "get_target_interactions_by_name",
    "Description": "Retrieve interaction data for a specific target, including interaction partners and evidence.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 586,
    "Tool Name": "get_disease_ancestors_parents_by_name",
    "Description": "Retrieve the ancestors and parents of a specific disease.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 587,
    "Tool Name": "get_disease_descendants_children_by_name",
    "Description": "Retrieve the descendants and children of a specific disease.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 588,
    "Tool Name": "get_disease_locations_by_name",
    "Description": "Retrieve the locations of a specific disease.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 589,
    "Tool Name": "get_disease_synonyms_by_name",
    "Description": "Retrieve synonyms for a specific disease.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 590,
    "Tool Name": "get_disease_description_by_name",
    "Description": "Retrieve the description of a specific disease.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 591,
    "Tool Name": "get_disease_therapeutic_areas_by_name",
    "Description": "Retrieve the therapeutic areas associated with a specific disease.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 592,
    "Tool Name": "get_chemical_probes_by_target_name",
    "Description": "Retrieve chemical probes associated with a specific target.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 593,
    "Tool Name": "get_associated_drugs_by_target_name",
    "Description": "Get known drugs associated with a specific target, including clinical trial phase and mechanism of action of the drugs.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 594,
    "Tool Name": "get_associated_diseases_by_drug_name",
    "Description": "Retrieve the list of diseases associated with a specific drug based on clinical trial data or post-marketed drugs.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 595,
    "Tool Name": "get_associated_targets_by_drug_name",
    "Description": "Retrieve the list of targets linked to a specific drug based on its mechanism of action.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 596,
    "Tool Name": "get_target_constraint_info_by_name",
    "Description": "Retrieve genetic constraint information for a specific target, including expected and observed values, and scores.",
    "category": "Databases",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 597,
    "Tool Name": "get_publications_by_disease_name",
    "Description": "Retrieve publications related to a disease name, including PubMed IDs and publication dates.",
    "category": "Literature Search",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 598,
    "Tool Name": "get_publications_by_target_name",
    "Description": "Retrieve publications related to a target, including PubMed IDs and publication dates.",
    "category": "Literature Search",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 599,
    "Tool Name": "get_publications_by_drug_name",
    "Description": "Retrieve publications related to a drug, including PubMed IDs and publication dates.",
    "category": "Literature Search",
    "Server Name": "Origene-OpenTargets"
  },
  {
    "IDX": 600,
    "Tool Name": "get_joint_associated_diseases_by_HPO_ID_list",
    "Description": "Retrieve diseases associated with a list of phenotypes or symptoms by a list of HPO IDs.",
    "category": "Databases",
    "Server Name": "Origene-Monarch"
  },
  {
    "IDX": 601,
    "Tool Name": "get_phenotype_by_HPO_ID",
    "Description": "Retrieve a phenotype or symptom by its HPO ID.",
    "category": "Databases",
    "Server Name": "Origene-Monarch"
  },
  {
    "IDX": 602,
    "Tool Name": "get_HPO_ID_by_phenotype",
    "Description": "Retrieve one or more HPO ID of a phenotype or symptom.",
    "category": "Databases",
    "Server Name": "Origene-Monarch"
  },
  {
    "IDX": 603,
    "Tool Name": "interproscan_analyze",
    "Description": "Analyze protein sequence using InterProScan to identify functional domains, families, and GO terms.\n\nUse this tool when you need to:\n- Identify protein domains and families (Pfam, Gene3D, SUPERFAMILY, etc.)\n- Get Gene Ontology (GO) annotations\n- Understand protein function and structure\n\nReturns: domains found, GO terms, and InterPro annotations.\nNote: Analysis may take 2-10 minutes depending on sequence length.",
    "category": "Databases",
    "Server Name": "BioInfo-Tools"
  },
  {
    "IDX": 604,
    "Tool Name": "blast_search",
    "Description": "Search for similar protein sequences in UniProt Swiss-Prot database using BLAST.\n\nUse this tool when you need to:\n- Find homologous proteins\n- Identify the protein by sequence similarity\n- Get organism and gene name information\n- Check if a sequence matches known proteins\n\nReturns: matching proteins with identity%, coverage%, E-value, organism, and alignment details.",
    "category": "Databases",
    "Server Name": "BioInfo-Tools"
  },
  {
    "IDX": 605,
    "Tool Name": "analyze_protein",
    "Description": "Comprehensive protein analysis combining InterProScan (domain analysis) and BLAST (similarity search).\n\nUse this tool when you need a complete analysis of an unknown protein sequence, including:\n- Domain and family identification\n- GO term annotations\n- Similar proteins in UniProt\n- Organism and gene information\n\nThis runs both analyses in parallel for efficiency.\nNote: May take 2-10 minutes to complete.",
    "category": "Computational Tools",
    "Server Name": "BioInfo-Tools"
  },
  {
    "IDX": 606,
    "Tool Name": "exec_code",
    "Description": "Complete code to perform a certain calculation\nParameters:\n    code_snippet: An executable code\nReturns:\n    str: Execution information\n\nExample of calculating the concentration of the original library:\n```python\nimport csv\nimport math\nimport os\n\n# CONFIG: adjust these values as needed\nINPUT_CSV = \"measurements.csv\"            # input CSV produced by measure(...)\nOUTPUT_CSV = \"measurements_with_conc.csv\" # output CSV to write\nSLOPE = -3.32                             # standard curve slope (Cq vs log10(conc))\nINTERCEPT = 38.5                          # standard curve intercept\nLOG_BASE = 10                             # base of logarithm used in curve\nDEFAULT_DILUTION = 20.0                   # fallback dilution factor if missing in CSV\nPARAMETER_COLUMN = \"parameter\"\nVALUE_COLUMN = \"measured_value\"\nDILUTION_COLUMN = \"dilution_factor\"\nCQ_INDICATORS = (\"cq\", \"ct\")             # strings to detect Cq rows (case-insensitive)\n\ndef _parse_num(s):\n    if s is None:\n        return None\n    s = str(s).strip()\n    if s == \"\":\n        return None\n    try:\n        return float(s)\n    except ValueError:\n        # extract numeric substring like \"19.3 pM\"\n        num = ''.join(ch for ch in s if (ch.isdigit() or ch in \".-eE\"))\n        return float(num) if num else None\n\n# Read input CSV\nwith open(INPUT_CSV, newline=\"\", encoding=\"utf-8\") as fh:\n    reader = csv.DictReader(fh)\n    rows = list(reader)\n    fieldnames = reader.fieldnames or []\n\n# Process rows: add two new columns\nfor r in rows:\n    r[\"calculated_concentration\"] = \"\"\n    r[\"corrected_concentration\"] = \"\"\n    param = (r.get(PARAMETER_COLUMN) or \"\").lower()\n    if any(k in param for k in CQ_INDICATORS):\n        cq = _parse_num(r.get(VALUE_COLUMN))\n        if cq is None:\n            r[\"calculated_concentration\"] = \"\"\n            r[\"corrected_concentration\"] = \"\"\n        else:\n            if abs(SLOPE) < 1e-12:\n                calc = float(\"nan\")\n            else:\n                exponent = (cq - INTERCEPT) / SLOPE\n                try:\n                    calc = math.pow(LOG_BASE, exponent)\n                except OverflowError:\n                    calc = float(\"inf\")\n            # dilution\n            df = _parse_num(r.get(DILUTION_COLUMN))\n            if df is None:\n                df = DEFAULT_DILUTION\n            corrected = calc * df if (calc is not None and not (isinstance(calc, float) and math.isnan(calc))) else \"\"\n            r[\"calculated_concentration\"] = \"\" if calc in (None, \"\") else (\"{:.6g}\".format(calc) if isinstance(calc, (int, float)) else str(calc))\n            r[\"corrected_concentration\"] = \"\" if corrected in (None, \"\") else (\"{:.6g}\".format(corrected) if isinstance(corrected, (int, float)) else str(corrected))\n\n# Write output CSV (preserve original fields, then append new ones)\nout_dir = os.path.dirname(OUTPUT_CSV) or \".\"\nos.makedirs(out_dir, exist_ok=True)\nextra_cols = [\"calculated_concentration\", \"corrected_concentration\"]\nout_fieldnames = list(fieldnames) + [c for c in extra_cols if c not in fieldnames]\n\nwith open(OUTPUT_CSV, \"w\", newline=\"\", encoding=\"utf-8\") as fh:\n    writer = csv.DictWriter(fh, fieldnames=out_fieldnames)\n    writer.writeheader()\n    for r in rows:\n        writer.writerow({k: r.get(k, \"\") for k in out_fieldnames})\n\nprint(f\"Processed {len(rows)} rows. Output written to: {OUTPUT_CSV}\")\n```",
    "category": "Computational Tools",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 607,
    "Tool Name": "store",
    "Description": "Manually adjust storage condition for a container.\nParameters:\n    container (Container or str): Container to adjust.\n    condition (str): Storage condition name.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 608,
    "Tool Name": "spin",
    "Description": "Apply acceleration (spin) to a container.\nParameters:\n    container (Container or str): Container to spin.\n    acceleration (str): e.g., \"1000:g\"\n    duration (QuantityString or str): Duration, e.g., \"20:minute\"\n    flow_direction (str, optional): \"inward\" or \"outward\"\n    spin_direction (List[str], optional): List like [\"cw\", \"ccw\"]\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 609,
    "Tool Name": "incubate",
    "Description": "Incubate a container at a location for a duration.\nParameters:\n    container (Container or str): Container to incubate.\n    where (str): Location (e.g., \"warm_37\").\n    duration (QuantityString or str): Time to incubate.\n    shaking (bool): Whether to shake.\n    co2 (float): CO2 percent.\n    uncovered (bool): If incubating uncovered.\n    target_temperature (QuantityString or str, optional): Target device temperature.\n    shaking_params (dict, optional): Shaking parameters.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 610,
    "Tool Name": "seal",
    "Description": "Seal a container.\nParameters:\n    container (Container or str): Container to seal.\n    type (str, optional): Seal type (e.g., \"ultra-clear\").\n    mode (str, optional): \"thermal\" or \"adhesive\".\n    temperature (QuantityString or str, optional): Temperature for thermal seal.\n    duration (QuantityString or str, optional): Duration to press seal.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 611,
    "Tool Name": "cover",
    "Description": "Place a lid on a container.\nParameters:\n    container (Container or str): Container to cover.\n    lid (str, optional): Lid type.\n    retrieve_lid (bool, optional): Whether to retrieve lid from storage.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 612,
    "Tool Name": "wash",
    "Description": "Perform washing of a target.\nParameters:\n    target (Container or str): Washing target.\n    solution (str): Solution name.\n    volume (QuantityString or str, optional): Volume used per repetition.\n    repetitions (int): Number of repetitions.\n    duration (QuantityString or str, optional): Duration of each wash.\n    temperature (QuantityString or str, optional): Temperature.\n    method (str, optional): Method e.g., \"gentle\".\n    speed (QuantityString or str, optional): Speed (e.g., \"1000:rpm\").\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 613,
    "Tool Name": "dry",
    "Description": "Dry an item by a method.\nParameters:\n    item (str): Item to dry.\n    method (str): Drying method (e.g., \"heat\").\n    duration (QuantityString or str, optional): Duration.\n    temperature (QuantityString or str, optional): Temperature.\n    location (str, optional): Location/device.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 614,
    "Tool Name": "close",
    "Description": "Close an item (e.g., cap or seal).\nParameters:\n    item (Container or str): Item to close.\n    method (str, optional): Method like \"screw_cap\".\n    security (str, optional): Security level.\n    cap_type (str, optional): Cap type.\n    seal_type (str, optional): Seal type.\n    temperature (QuantityString or str, optional): Temperature.\n    ensure_proper_closure (bool): Ensure proper closure.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 615,
    "Tool Name": "check",
    "Description": "Perform a general validation/check operation.\nParameters:\n    item (str, optional): Item being checked.\n    parameter (str, optional): Parameter name.\n    expected_value (Any, optional): Expected value.\n    actual_value (Any, optional): Actual value.\n    condition (str, optional): Condition type.\n    method (str, optional): Method used.\n    tolerance (str, optional): Tolerance with units.\n    criteria (str/dict/list, optional): Criteria standard.\n    samples (list, optional): List of sample containers.\n    equipment (str, optional): Equipment used.\nReturns:\n    str: Execution information",
    "category": "Computational Tools",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 616,
    "Tool Name": "filter",
    "Description": "Filter a sample through a filter type.\nParameters:\n    sample (str): Sample name or id.\n    filter_type (str): e.g., \"syringe_filter\".\n    pore_size (str, optional): Pore size with units e.g., \"0.22:micrometer\".\n    volume (QuantityString or str, optional): Volume to filter.\n    method (str, optional): Filtration method.\n    container (Container or str, optional): Destination container.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 617,
    "Tool Name": "record",
    "Description": "Record data of a specified type.\nParameters:\n    data_type (str): Type of measurement (e.g., \"absorbance\").\n    value (str): Value with units.\n    sample (Container or str, optional): Sample container.\n    location (str, optional): Location like \"A1\".\n    method (str, optional): Measurement method.\n    details (str/dict, optional): Additional details.\n    instrument (str, optional): Instrument id.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 618,
    "Tool Name": "stand",
    "Description": "Let an item stand for some duration.\nParameters:\n    item (Container or str): Item to stand.\n    duration (QuantityString or str): Duration.\n    temperature (QuantityString or str, optional): Temperature.\n    location (str, optional): Location.\n    position (str, optional): Orientation.\n    purpose (str, optional): Purpose like \"equilibration\".\n    stand_type (str, optional): Type of stand.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 619,
    "Tool Name": "thaw",
    "Description": "Thaw materials at a specific temperature.\nParameters:\n    materials (str or List[str]): Names of materials.\n    temperature (QuantityString or str): Thaw temperature.\n    duration (QuantityString or str, optional): Duration.\n    method (str, optional): Thaw method.\n    container (Container or str, optional): Container used.\n    location (str, optional): Location.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 620,
    "Tool Name": "aspirate",
    "Description": "Aspirate liquid from a source.\nParameters:\n    source (Container or str): Source container.\n    volume (QuantityString or str, optional): Volume to aspirate.\n    liquid (str, optional): Liquid description.\n    tool (str, optional): Tool name.\n    technique (str, optional): Technique like \"reverse_pipetting\".\n    speed (str, optional): Speed e.g., \"slow\".\n    angle (int, optional): Angle in degrees.\n    note (str, optional): Notes.\n    destination (Container or str, optional): Destination container.\n    completeness (str, optional): \"complete\" or similar.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 621,
    "Tool Name": "thermocycle",
    "Description": "Thermocycle a container with groups of steps.\nParameters:\n    container (Container or str): Container to thermocycle.\n    groups (List[dict]): List of group dicts with cycles and steps.\n    volume (QuantityString or str, optional): Volume per well.\n    dataref (str, optional): Dataref name for qPCR.\n    dyes (dict, optional): Dye mapping to wells.\n    melting_* (QuantityString, optional): Melting curve parameters.\n    lid_temperature (QuantityString or str, optional): Lid temp.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 622,
    "Tool Name": "wait",
    "Description": "Wait for a specified duration.\nParameters:\n    duration (QuantityString or str): Duration to wait.\n    reason (str, optional): Reason for waiting.\n    purpose (str, optional): Purpose.\n    condition (str, optional): Environmental condition.\n    activity (str, optional): Activity during wait.\n    sample (Container or str, optional): Sample being waited on.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 623,
    "Tool Name": "heat",
    "Description": "Heat a target to a temperature for a duration.\nParameters:\n    target (Container or str): Target to heat.\n    temperature (QuantityString or str): Target temperature.\n    duration (QuantityString or str): Heating duration.\n    method (str, optional): Method.\n    equipment (str, optional): Equipment used.\n    purpose (str, optional): Purpose.\n    volume (QuantityString or str, optional): Volume involved.\n    power (float or str, optional): Power setting.\n    agitation_speed (int, optional): Agitation speed.\n    lid_temperature (QuantityString or str, optional): Lid temperature.\n    lid_status (str, optional): \"open\" or \"closed\".\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 624,
    "Tool Name": "vortex",
    "Description": "Vortex an item.\nParameters:\n    item (str or Container): Item to vortex.\n    duration (QuantityString or str): Duration.\n    speed (str, optional): Speed setting.\n    intensity (str, optional): Intensity.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 625,
    "Tool Name": "insert",
    "Description": "Insert an item into a target.\nParameters:\n    item (str): Item to insert.\n    target (str): Target location.\n    position (str, optional): Position within target.\n    method (str, optional): Method used.\n    orientation (str, optional): Orientation of insertion.\n    alignment (str, optional): Alignment.\n    precautions (str, optional): Precautions.\nReturns:\n    str: Execution information",
    "category": "Computational Tools",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 626,
    "Tool Name": "separate",
    "Description": "Separate phases/components from a target.\nParameters:\n    target (str): Separation target.\n    method (str): Separation method.\n    duration (QuantityString or str, optional): Duration.\n    phases (List[str], optional): Named phases.\n    equipment (str, optional): Equipment used.\n    conditions (dict, optional): Additional conditions.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 627,
    "Tool Name": "lyse",
    "Description": "Lyse a sample using a method and optional buffer.\nParameters:\n    sample (Container or str): Sample to lyse.\n    method (str): Lysis method.\n    buffer (str, optional): Lysis buffer.\n    duration (QuantityString or str, optional): Duration.\n    temperature (QuantityString or str, optional): Temperature.\n    volume (QuantityString or str, optional): Volume.\n    reagents (dict, optional): Additional reagents.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 628,
    "Tool Name": "add",
    "Description": "Add a component to a target container.\nParameters:\n    component (str): Component name.\n    volume (QuantityString or str): Volume to add.\n    target (Container or str): Destination container.\n    concentration (str, optional): Concentration string.\n    weight (QuantityString or str): Weight to add\n    source (Container or str, optional): Source container.\n    method (str, optional): Addition method.\n    mixing_required (bool): Whether to mix after adding.\n    avoid_wall (bool): Avoid container wall.\n    use_new_tip (bool): Use new tip for operation.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 629,
    "Tool Name": "inoculate",
    "Description": "Inoculate medium with inoculum.\nParameters:\n    medium (Container or str): Medium container.\n    inoculum (Container or str): Inoculum source.\n    volume (QuantityString or str, optional): Volume to inoculate.\n    method (str, optional): Method like \"pipetting\".\n    conditions (str, optional): Cultivation conditions.\n    antibiotic (str, optional): Antibiotic name.\n    antibiotic_concentration (str, optional): Antibiotic concentration.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 630,
    "Tool Name": "inactivate",
    "Description": "Inactivate a target using specified method and parameters.\nParameters:\n    target (Container or str): Target to inactivate.\n    method (str): Inactivation method.\n    duration (QuantityString or str, optional): Duration.\n    temperature (QuantityString or str, optional): Temperature.\n    volume (QuantityString or str, optional): Volume.\n    agent (str, optional): Inactivation agent.\n    concentration (str, optional): Agent concentration.\n    container (Container or str, optional): Container used.\n    sample_type (str, optional): Type of sample.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 631,
    "Tool Name": "anneal",
    "Description": "Anneal material at a temperature for a duration with optional cycles.\nParameters:\n    material (Container or str): Material to anneal.\n    temperature (QuantityString or str): Annealing temperature.\n    duration (QuantityString or str): Duration.\n    template (Container or str, optional): Template material.\n    cycles (int, optional): Number of cycles.\n    program_type (str, optional): Program type like \"gradient\".\n    temperature_decrease (str, optional): Decrease rate per cycle.\n    method (str, optional): Annealing method.\n    equipment (str, optional): Equipment used.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 632,
    "Tool Name": "transfer",
    "Description": "Transfer material from source to destination.\nParameters:\n    source (Container or str): Source location.\n    destination (Container or str): Destination location.\n    volume (QuantityString or str, optional): Volume to transfer.\n    method (str, optional): Transfer method.\n    tool (str, optional): Tool used.\n    material (str, optional): Material description.\n    use_filter_tip (bool): Use filter tip flag.\n    avoid_contact (bool): Avoid contact flag.\n    temperature (QuantityString or str, optional): Temperature during transfer.\n    notes (str, optional): Additional notes.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 633,
    "Tool Name": "resuspend",
    "Description": "Resuspend a material/pellet in buffer.\nParameters:\n    material (str): Material name.\n    buffer (str): Buffer name.\n    volume (QuantityString or str): Volume to resuspend.\n    method (str, optional): Resuspension method.\n    container (Container or str, optional): Container where resuspension occurs.\n    duration (QuantityString or str, optional): Duration.\n    speed (QuantityString or str, optional): Speed setting.\n    mix_times (int, optional): Number of mixing repetitions.\n    concentration (str, optional): Target concentration.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 634,
    "Tool Name": "open",
    "Description": "Open an item.\nParameters:\n    item (Container or str): Item to open.\n    method (str, optional): Opening method.\nReturns:\n    str: Execution information",
    "category": "Computational Tools",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 635,
    "Tool Name": "mix",
    "Description": "Mix components in a target container.\nParameters:\n    target (Container or str, optional): Container to mix.\n    method (str): Mixing method.\n    duration (QuantityString or str, optional): Duration.\n    components (list, optional): Components list.\n    volume (QuantityString or str, optional): Volume.\n    repetitions (int, optional): Number of repetitions.\n    intensity (str, optional): Intensity level.\n    ratio (str, optional): Mixing ratio.\n    container (Container or str, optional): Alternative container.\n    temperature (QuantityString or str, optional): Temperature.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 636,
    "Tool Name": "extract",
    "Description": "Extract target material from a sample.\nParameters:\n    sample (Container or str): Sample source.\n    method (str): Extraction method.\n    volume (QuantityString or str, optional): Volume used.\n    buffer (str, optional): Buffer used.\n    kit (str, optional): Kit name.\n    solvent (str, optional): Solvent used.\n    target (str, optional): Target molecule.\n    protocol (str, optional): Protocol name.\n    duration (QuantityString or str, optional): Duration.\n    temperature (QuantityString or str, optional): Temperature.\n    replicates (int, optional): Number of replicates.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 637,
    "Tool Name": "grind",
    "Description": "Grind material with specified method and parameters.\nParameters:\n    material (str): Material to grind.\n    method (str): Grinding method.\n    duration (QuantityString or str, optional): Duration.\n    tool (str, optional): Tool name.\n    buffer_type (str, optional): Buffer type.\n    buffer_volume (QuantityString or str, optional): Buffer volume.\n    temperature (QuantityString or str, optional): Temperature.\n    speed (str, optional): Speed setting.\n    equipment (str, optional): Equipment name.\n    grinding_medium (str, optional): Medium used for grinding.\n    container (Container or str, optional): Container used.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 638,
    "Tool Name": "prepare",
    "Description": "Prepare reagents or mixtures.\nParameters:\n    item (str): Name of item to prepare.\n    volume (QuantityString or str, optional): Volume.\n    concentration (str, optional): Concentration string.\n    components (dict, optional): Components mapping.\n    container (Container or str, optional): Container to hold prepared item.\n    temperature (QuantityString or str, optional): Preparation temperature.\n    quantity (int, optional): Number of aliquots.\n    storage_condition (str, optional): Storage condition.\n    purpose (str, optional): Intended use.\n    method (str, optional): Preparation method.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 639,
    "Tool Name": "freeze",
    "Description": "Freeze an item at given temperature.\nParameters:\n    item (Container or str): Item to freeze.\n    temperature (QuantityString or str): Freeze temperature.\n    duration (QuantityString or str, optional): Duration of freezing.\n    method (str, optional): Freezing method.\n    container (Container or str, optional): Storage container.\n    location (str, optional): Location or position.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 640,
    "Tool Name": "pipetting",
    "Description": "Perform pipetting from source to destination.\nParameters:\n    volume (QuantityString or str): Volume to pipette.\n    source (Container or str): Source location.\n    destination (Container or str): Destination location.\n    liquid (str, optional): Liquid type.\n    mix_times (int, optional): Number of mixes.\n    tip_change (bool): Whether to change tips.\n    use_new_tip (bool): Use new tip flag.\n    avoid_bubbles (bool): Avoid bubbles flag.\n    technique (str, optional): Technique name.\n    repetitions (int, optional): Repetitions.\n    pipette_type (str, optional): Pipette model.\n    tip_type (str, optional): Tip type.\n    action (str, optional): Action type.\n    sample_type (str, optional): Sample type.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 641,
    "Tool Name": "extend",
    "Description": "Perform extension (e.g., polymerase extension).\nParameters:\n    material (Container or str, optional): Material/sample.\n    temperature (QuantityString or str, optional): Extension temp.\n    duration (QuantityString or str, optional): Duration.\n    enzyme (str, optional): Enzyme name.\n    template (Container or str, optional): Template container.\n    cycles (int, optional): Number of cycles.\n    duration_per_kb (QuantityString or str, optional): Time per kb.\n    target_length (float, optional): Target length in kb.\n    method (str, optional): Method name.\n    equipment (str, optional): Equipment used.\n    data_collection (bool): Whether to collect data.\n    optimization_range (str, optional): Temp optimization range.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 642,
    "Tool Name": "invert",
    "Description": "Invert a container repeatedly to mix.\nParameters:\n    container (Container or str): Container to invert.\n    repetitions (int): Number of inversions.\n    purpose (str): Purpose like \"mix\".\n    method (str, optional): Inversion method.\n    duration (QuantityString or str, optional): Duration.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 643,
    "Tool Name": "place",
    "Description": "Place an item at a location.\nParameters:\n    item (str or Container): Item to place.\n    location (str): Target location.\n    duration (QuantityString or str, optional): Duration to remain.\n    temperature (QuantityString or str, optional): Temperature condition.\n    orientation (str, optional): Orientation.\n    equipment (str, optional): Equipment used.\n    purpose (str, optional): Purpose of placement.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 644,
    "Tool Name": "fill",
    "Description": "Fill a container with a substance.\nParameters:\n    container (Container or str): Container to fill.\n    substance (str): Substance name.\n    volume (QuantityString or str): Volume to fill.\n    method (str, optional): Filling method.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 645,
    "Tool Name": "tap",
    "Description": "Tap or flick a container to mix or settle.\nParameters:\n    container (Container or str): Container to tap.\n    intensity (str, optional): \"gentle\"/\"medium\"/\"hard\".\n    duration (QuantityString or str, optional): Duration.\n    purpose (str, optional): Purpose of tapping.\n    times (int, optional): Number of taps.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 646,
    "Tool Name": "centrifuge",
    "Description": "Centrifuge a container at a speed for a duration.\nParameters:\n    container (Container or str): Container to centrifuge.\n    speed (QuantityString or str): Speed e.g., \"3000:rpm\".\n    duration (QuantityString or str): Duration.\n    temperature (QuantityString or str, optional): Temperature.\n    equipment (str, optional): Equipment identifier.\n    purpose (str, optional): Purpose.\n    notes (str, optional): Notes.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 647,
    "Tool Name": "pellet",
    "Description": "Pellet material by specified method (e.g., centrifugation).\nParameters:\n    sample (Container or str, optional): Sample to pellet.\n    method (str, optional): Pelleting method.\n    speed (QuantityString or str, optional): Speed.\n    duration (QuantityString or str, optional): Duration.\n    container (Container or str, optional): Container used.\n    temperature (QuantityString or str, optional): Temperature.\n    visible (bool, optional): Whether pellet is visible.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 648,
    "Tool Name": "aliquot",
    "Description": "Aliquot a source into destination containers.\nParameters:\n    source (Container or str): Source container.\n    destination (Container or str): Destination container.\n    volume (QuantityString or str): Volume per aliquot.\n    number_of_aliquots (int): Number of aliquots.\n    material (str, optional): Material description.\n    storage_temperature (str, optional): Storage temperature.\n    label (str, optional): Label for aliquots.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 649,
    "Tool Name": "eliminate",
    "Description": "Eliminate a target (e.g., contaminant) using a method.\nParameters:\n    target (str): Target to eliminate.\n    method (str): Elimination method.\n    source (str, optional): Source container.\n    duration (QuantityString or str, optional): Duration.\n    reagents (List[str], optional): Reagents used.\n    volume (QuantityString or str, optional): Volume.\n    equipment (List[str], optional): Equipment list.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 650,
    "Tool Name": "measure",
    "Description": "Measure a parameter on a sample.\nParameters:\n    sample (Container or str or List): Sample to measure.\n    parameter (str): Parameter name e.g., \"absorbance\".\n    instrument (str, optional): Instrument id.\n    method (str, optional): Method used.\n    value (str, optional): Measured value with units.\n    equipment (str, optional): Equipment used.\n    container (Container or str, optional): Container reference.\n    output_path (str, optional): csv path of measurement\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 651,
    "Tool Name": "balance",
    "Description": "Balance items for equipment (e.g., centrifuge).\nParameters:\n    items (Container or List[Container]): Items to balance.\n    equipment (str, optional): Equipment used.\n    method (str, optional): Balancing method.\n    target_weight (str, optional): Target weight with units.\n    balance_item (Container or str, optional): Reference balance item.\n    duration (str, optional): Duration.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 652,
    "Tool Name": "shake",
    "Description": "Shake a target container.\nParameters:\n    target (Container or str): Target to shake.\n    duration (QuantityString or str, optional): Duration.\n    speed (QuantityString or str, optional): Speed e.g., \"300:rpm\".\n    temperature (QuantityString or str, optional): Temperature.\n    intensity (str, optional): \"gentle\"/\"vigorous\".\n    method (str, optional): Method like \"orbital\".\n    equipment (str, optional): Equipment id.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 653,
    "Tool Name": "decant",
    "Description": "Decant liquid from source to destination.\nParameters:\n    source (Container or str): Source container.\n    destination (Container or str, optional): Destination container.\n    volume (QuantityString or str, optional): Volume to decant.\n    liquid (str, optional): Liquid description e.g., \"supernatant\".\n    method (str, optional): Decanting method.\n    precautions (str, optional): Precautions to take.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 654,
    "Tool Name": "purify",
    "Description": "Purify a sample using a method/kit.\nParameters:\n    sample (Container or str): Sample to purify.\n    method (str): Purification method.\n    kit (str, optional): Kit name.\n    volume (QuantityString or str, optional): Volume.\n    reagents (str/list/dict, optional): Reagents used.\n    target (str, optional): Target molecule.\n    temperature (QuantityString or str, optional): Temperature.\n    time (QuantityString or str, optional): Time for step.\n    requirement (str, optional): Purity requirement.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 655,
    "Tool Name": "sequence",
    "Description": "Perform sequencing run with specified parameters.\nParameters:\n    kit (str, optional): Sequencing kit.\n    instrument (str, optional): Instrument name.\n    protocol (str, optional): Sequencing protocol.\n    flow_cell (str, optional): Flow cell type.\n    format (str, optional): Read format.\n    facility (str, optional): Facility name.\n    product (str, optional): Product to sequence.\n    reagent_cartridge (str, optional): Reagent cartridge type.\n    spike_ratio (float, optional): Spike-in ratio.\n    core (str, optional): Core facility name.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 656,
    "Tool Name": "digest",
    "Description": "Digest a sample with an enzyme under specified conditions.\nParameters:\n    sample (Container or str): Sample to digest.\n    enzyme (str): Enzyme name.\n    buffer (str, optional): Buffer used.\n    temperature (QuantityString or str, optional): Temperature.\n    duration (QuantityString or str, optional): Duration.\n    volume (QuantityString or str, optional): Volume.\n    concentration (str, optional): Concentration.\n    enzyme_units (float, optional): Units of enzyme.\n    samples (List[Container], optional): Batch samples.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 657,
    "Tool Name": "elute",
    "Description": "Elute material from a sample using a buffer.\nParameters:\n    sample (Container or str): Source sample.\n    buffer (str): Elution buffer.\n    volume (QuantityString or str): Elution volume.\n    temperature (QuantityString or str, optional): Temperature.\n    duration (QuantityString or str, optional): Duration.\n    mixing_cycles (int, optional): Mixing cycles count.\n    container (Container or str, optional): Column/container source.\n    target_tube (Container or str, optional): Target tube for eluate.\n    speed (QuantityString or str, optional): Speed used.\n    method (str, optional): Elution method.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 658,
    "Tool Name": "stain",
    "Description": "Stain a target with an agent.\nParameters:\n    target (Container or str): Target to stain.\n    agent (str): Staining agent.\n    duration (QuantityString or str, optional): Duration.\n    volume (QuantityString or str, optional): Volume of stain.\n    concentration (str, optional): Concentration of stain.\n    method (str, optional): Staining method.\n    temperature (QuantityString or str, optional): Temperature.\n    light_sensitive (bool): Protect from light if True.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 659,
    "Tool Name": "discard",
    "Description": "Discard items or materials safely.\nParameters:\n    item (str, optional): Item name.\n    method (str, optional): Discard method.\n    container (Container or str, optional): Container from which to discard.\n    volume (QuantityString or str, optional): Volume.\n    waste_type (str, optional): Waste category.\n    safety (str, optional): Safety instructions.\n    caution (bool): Whether caution is required.\n    decontamination (bool): Whether decontamination is required.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 660,
    "Tool Name": "dilute",
    "Description": "Dilute a sample by a factor.\nParameters:\n    sample (Container or str): Sample to dilute.\n    dilution_factor (float): Factor of dilution.\n    diluent (str, optional): Diluent name.\n    volume (QuantityString or str, optional): Volume parameter.\n    target_concentration (str, optional): Target concentration.\n    final_volume (QuantityString or str, optional): Final volume.\n    source_volume (QuantityString or str, optional): Source volume used.\n    aliquots (int, optional): Number of aliquots.\n    mixing_instructions (str, optional): Mixing instructions.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 661,
    "Tool Name": "dissolve",
    "Description": "Dissolve material in a solvent under specified conditions.\nParameters:\n    material (str): Material to dissolve.\n    solvent (str): Solvent to use.\n    volume (QuantityString or str, optional): Volume of solvent.\n    concentration (str, optional): Target concentration.\n    method (str, optional): Dissolution method.\n    temperature (QuantityString or str, optional): Temperature.\n    duration (QuantityString or str, optional): Total duration.\n    shaking_time (QuantityString or str, optional): Shaking time.\n    mixing_times (int, optional): Number of mixing cycles.\n    until_condition (str, optional): Condition to stop (e.g., \"fully dissolved\").\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 662,
    "Tool Name": "software_analysis",
    "Description": "Analyze real-time PCR run results using a specified analysis software.\nParameters:\n    software_name (str): Name of the analysis software (e.g., \"Sequence Detector\").\n    version (str, optional): Software version (e.g., \"v2.3\").\n    instrument (str, optional): Instrument used to generate the data (e.g., \"7500 Real-Time PCR-System\").\n    analysis_steps (List[str], optional): Ordered human-readable steps the software will perform (e.g., [\"import data\", \"apply baseline correction\", \"set threshold\", \"call Cq\"]).\n    channel_mapping (Dict[str, str], optional): Mapping of fluorophore channels to target/probe names (e.g., {\"FAM\": \"DnAprTm-v\", \"VIC\": \"DnAprTM-b\", \"LIZ\": \"Xeno\"}).\n    threshold (str, optional): Numeric threshold value for Cq calling.\n    baseline_correction (bool, optional): Whether automatic baseline correction is applied.\n    baseline_cycles (str, optional): Baseline cycles definition (e.g., \"cycles 3-15\" or \"automatic\").\n    replicate_handling (str, optional): How replicates are treated (e.g., \"average Cq\", \"report individual wells\", \"flag discordant\").\n    internal_controls (List[str], optional): Names of internal controls to check (e.g., [\"Xeno? Internal Positive Control\"]).\n    standard_curve_used (bool, optional): Whether an external standard curve was run for absolute quantification.\n    quantification_method (str, optional): \"absolute\" or \"relative\" or other description of quantification approach.\n    output_format (str, optional): Desired export format(s) (e.g., \"CSV\", \"XLSX\", \"JSON\").\n    notes (str, optional): Free-text field for any additional settings (e.g., threshold rationale, software-specific flags).\nReturns:\n    str: Execution information",
    "category": "Computational Tools",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 663,
    "Tool Name": "run_gel_and_image",
    "Description": "Abstract instruction to run an agarose gel and acquire an image using an appropriate DNA stain and imaging system.\n\nThis function captures the parameters relevant to electrophoresis and imaging only. It focuses on the single\nstatement: \"Run the gel and image using an appropriate DNA staining/imaging system. This scheme uses gels\nstained with ethidium bromide and images are acquired on a BioRad GelDoc XR system.\" Low-level laboratory\nactions (e.g., preparing the gel, loading samples, turn on/off instruments, waste disposal) are intentionally\nout of scope for this instruction and must be performed according to local laboratory procedures and safety rules.\n\nParameters:\n    samples (Container | str | list, optional): Identifier(s) for the samples to be loaded and imaged.\n    gel_percentage (float, optional): Agarose percentage of the gel (e.g., 1.0 for 1% agarose). Default: 1.0.\n    gel_buffer (str, optional): Buffer used for gel casting and electrophoresis (e.g., \"TAE\" or \"TBE\"). Default: \"TAE\".\n    gel_volume (QuantityString, optional): Informational field for gel volume or cast dimensions.\n    well_load_volume (QuantityString, optional): Total volume loaded per well (sample + loading dye). Default: \"16 ¦Ìl\".\n    loading_buffer (str, optional): Type of loading dye. Default: \"6¡Á DNA loading buffer\".\n    loading_buffer_volume_per_sample (QuantityString, optional): Volume of loading buffer added to each sample. Default: \"4 ¦Ìl\".\n    positive_control_volume_final (QuantityString, optional): Recommended load volume for a 1 ng positive control plasmid to avoid gel distortion during imaging. Default: \"8 ¦Ìl\".\n    ladder (str, optional): Identifier or description of the DNA ladder/marker used.\n    electrophoresis_voltage (QuantityString, optional): Voltage applied during electrophoresis (e.g., \"100 V\") ¡ª for record-keeping.\n    run_time (QuantityString, optional): Electrophoresis run time (e.g., \"45 min\") ¡ª for record-keeping.\n    stain (str, optional): DNA stain used for visualization. Default: \"Ethidium Bromide\".\n    imaging_system (str, optional): Imaging system used to capture gel images. Default: \"BioRad GelDoc XR\".\n    imaging_filter (str, optional): Imaging filter/channel appropriate for the stain (e.g., \"EtBr/UV\").\n    exposure_time (QuantityString, optional): Suggested or recorded exposure time for imaging.\n    image_output_path (str, optional): Path or filename to save the acquired gel image.\n    notes (str, optional): Free-text field for deviations, safety notes, or system-specific adjustments.\nReturns:\n    str: Execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-OP"
  },
  {
    "IDX": 664,
    "Tool Name": "protocol_generation",
    "Description": "Generation experiment protocol given user prompt.\nParams£º\n    user_prompt (str): requirements for the experiment\n\nReturns£º\n    protocol (str): generated protocol text",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-Plan"
  },
  {
    "IDX": 665,
    "Tool Name": "generate_executable_json",
    "Description": "Generate executable JSON file from protocol text.\nParams:\n    protocol (str): protocol text\n\nReturns:\n    json (str): executable JSON string",
    "category": "Computational Tools",
    "Server Name": "Thoth-Plan"
  },
  {
    "IDX": 666,
    "Tool Name": "extract_protocol_from_pdf",
    "Description": "Extract experimental protocol from a provided PDF\nParams:\n    pdf_url (str): PDF file URL\n\nReturns:\n    protocol (str): protocol text",
    "category": "Computational Tools",
    "Server Name": "Thoth-Plan"
  },
  {
    "IDX": 667,
    "Tool Name": "execute_json",
    "Description": "Execute JSON with PCR operation server\nParams:\n    json (str): executable JSON\n    server_url (str): PCR operation server URL; default to Thoth-OP server URL\n\nReturns:\n    exec_info (str): execution information",
    "category": "Wet-lab Operations",
    "Server Name": "Thoth-Plan"
  },
  {
    "IDX": 668,
    "Tool Name": "calculate_geometric_term",
    "Description": "Calculate the geometric term ¡Ì(¦Ða).",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 669,
    "Tool Name": "calculate_stress_intensity_factor",
    "Description": "Calculate the stress intensity factor KI.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 670,
    "Tool Name": "determine_fracture",
    "Description": "Determine if fracture occurs.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 671,
    "Tool Name": "calculate_mobility_ratio",
    "Description": "Calculate the ratio of the new mobility relative to the original mobility.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 672,
    "Tool Name": "calculate_cell_volume",
    "Description": "Calculate the unit cell volume of a crystal with lattice constant a.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 673,
    "Tool Name": "calculate_density",
    "Description": "Calculate density as mass divided by volume.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 674,
    "Tool Name": "calculate_proportionality_constant",
    "Description": "Calculate the proportionality constant K from known diameter, length, and weight.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 675,
    "Tool Name": "calculate_ring_moment_of_inertia",
    "Description": "Calculate the moment of inertia of a ring.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 676,
    "Tool Name": "calculate_spoke_moment_of_inertia",
    "Description": "Calculate the moment of inertia of a spoke (rod) about one end.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 677,
    "Tool Name": "calculate_total_spokes_inertia",
    "Description": "Calculate total moment of inertia of all spokes.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 678,
    "Tool Name": "calculate_weight",
    "Description": "Calculate the weight of concrete.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 679,
    "Tool Name": "convert_crack_length_to_meters",
    "Description": "Convert crack length to meters based on the specified unit.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 680,
    "Tool Name": "calculate_critical_stress_griffith_irwin",
    "Description": "Calculate critical stress using Griffith-Irwin equation.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 681,
    "Tool Name": "validate_strain_input",
    "Description": "Validate that the axial strain is non-zero.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 682,
    "Tool Name": "calculate_material_density",
    "Description": "Compute the material density based on specific gravity.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 683,
    "Tool Name": "calculate_volume_in_cubic_meters",
    "Description": "Calculate the volume in cubic meters given mass and density.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 684,
    "Tool Name": "calculate_angular_frequency",
    "Description": "Calculate the angular frequency of a harmonic oscillator.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 685,
    "Tool Name": "calculate_density_value",
    "Description": "Calculate density in g/cm3.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 686,
    "Tool Name": "convert_GPa_to_MPa",
    "Description": "Convert Young's modulus from GPa to MPa.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 687,
    "Tool Name": "determine_exponent_by_dimension",
    "Description": "Calculate the exponent used in the Mott VRH model based on system dimension.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 688,
    "Tool Name": "calculate_magnetic_anisotropy_energy",
    "Description": "Calculate the magnetic anisotropy energy.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 689,
    "Tool Name": "convert_gpa_to_pa",
    "Description": "Convert elastic modulus from GPa to Pa.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 690,
    "Tool Name": "calculate_stiffness_contribution",
    "Description": "Calculate the stiffness contribution (E * A) for a material.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 691,
    "Tool Name": "aggregate_stiffness",
    "Description": "Sum individual stiffness contributions to get total stiffness.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 692,
    "Tool Name": "calculate_stress_in_Pa",
    "Description": "Calculate tensile stress in Pascals.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 693,
    "Tool Name": "convert_Pa_to_MPa",
    "Description": "Convert stress from Pascals to Megapascals.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 694,
    "Tool Name": "calculate_volume_from_mass",
    "Description": "Calculate the volume of the oil droplet.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 695,
    "Tool Name": "calculate_packing_factor",
    "Description": "Calculate atomic packing factor (APF).",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 696,
    "Tool Name": "calculate_density_difference_percentage",
    "Description": "Calculate percentage difference between two packing factors.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 697,
    "Tool Name": "validate_stress_unit",
    "Description": "Validate that the stress unit is within the allowed units.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 698,
    "Tool Name": "calculate_elastic_modulus_rule_of_mixtures",
    "Description": "Calculate elastic modulus of alloy using rule of mixtures.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 699,
    "Tool Name": "compute_moment_of_inertia",
    "Description": "Calculate the moment of inertia of a hollow sphere about its diameter.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 700,
    "Tool Name": "calculate_hall_petch_strength",
    "Description": "Calculate yield strength using Hall-Petch equation.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 701,
    "Tool Name": "calculate_strain",
    "Description": "Calculate the strain (dimensionless) from elongation and original length.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 702,
    "Tool Name": "calculate_energy_density",
    "Description": "Calculate the elastic energy density (J/m3) from Young's modulus and strain.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 703,
    "Tool Name": "calculate_total_energy",
    "Description": "Calculate the total elastic energy stored in the material.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 704,
    "Tool Name": "convert_stress_to_pa",
    "Description": "Convert stress to Pascals (Pa).",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 705,
    "Tool Name": "convert_pa_to_output_unit",
    "Description": "Convert a stress value from Pascals to the specified output unit.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 706,
    "Tool Name": "compute_compression_ratio",
    "Description": "Compute the compression ratio.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 707,
    "Tool Name": "calculate_mass_from_volume_density",
    "Description": "Calculate mass from volume and density.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 708,
    "Tool Name": "calculate_weight_N",
    "Description": "Calculate weight in Newtons.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 709,
    "Tool Name": "calculate_buoyancy_force_N",
    "Description": "Calculate buoyant force in Newtons.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 710,
    "Tool Name": "calculate_thickness_nm",
    "Description": "Calculate total thickness in nanometers.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 711,
    "Tool Name": "compute_transverse_strain_from_poisson",
    "Description": "Calculate the transverse strain using Poisson's ratio and longitudinal strain.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 712,
    "Tool Name": "calculate_length_change",
    "Description": "Calculate the length change due to thermal expansion.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 713,
    "Tool Name": "calculate_final_length",
    "Description": "Calculate the final length after thermal contraction.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 714,
    "Tool Name": "calculate_relative_length_change",
    "Description": "Calculate the relative change in length as a percentage.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 715,
    "Tool Name": "calculate_minimum_thickness_for_waterproofing",
    "Description": "Calculate the minimum coating thickness needed for waterproofing (¦È ¡Ý 90¡ã).",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 716,
    "Tool Name": "calculate_atomic_radius_bcc",
    "Description": "Calculate the atomic radius for BCC structure.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 717,
    "Tool Name": "calculate_density_bcc",
    "Description": "Calculate the density of a BCC structured material.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 718,
    "Tool Name": "calculate_stress",
    "Description": "Calculate stress according to Hooke's law.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 719,
    "Tool Name": "check_debonding",
    "Description": "Determine whether interface debonding will occur.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 720,
    "Tool Name": "calculate_safety_factor",
    "Description": "Calculate safety factor.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 721,
    "Tool Name": "convert_mpa_to_pa",
    "Description": "Convert stress from MPa to Pascals.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 722,
    "Tool Name": "calculate_plastic_section_modulus",
    "Description": "Calculate the plastic section modulus Z for a rectangular cross-section.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 723,
    "Tool Name": "calculate_plastic_moment",
    "Description": "Calculate the plastic moment Mp in kN¡¤m.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 724,
    "Tool Name": "compute_density_value",
    "Description": "Compute the density from mass and volume.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 725,
    "Tool Name": "convert_GPa_to_Pa",
    "Description": "Convert shear modulus from GPa to Pa.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 726,
    "Tool Name": "calculate_solid_cylinder_inertia",
    "Description": "Calculate the moment of inertia of a hollow cylinder about its central axis.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 727,
    "Tool Name": "calculate_final_roughness",
    "Description": "Calculate the final surface roughness after reduction.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 728,
    "Tool Name": "calculate_hardness_increase_percent",
    "Description": "Calculate the percentage increase in hardness.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 729,
    "Tool Name": "calculate_impact_toughness_in_j_per_m2",
    "Description": "Calculate impact toughness in J/m^2.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 730,
    "Tool Name": "convert_j_per_m2_to_kj_per_m2",
    "Description": "Convert impact toughness from J/m^2 to kJ/m^2.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 731,
    "Tool Name": "calculate_solids_volume",
    "Description": "Calculate volume occupied by suspended solids (m3/m2).",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 732,
    "Tool Name": "calculate_stress_amplitude_value",
    "Description": "Calculate stress amplitude.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 733,
    "Tool Name": "calculate_mean_stress",
    "Description": "Calculate mean stress.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 734,
    "Tool Name": "calculate_stress_ratio",
    "Description": "Calculate stress ratio R = ¦Ò_min / ¦Ò_max.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 735,
    "Tool Name": "compute_numerator",
    "Description": "Compute the numerator for Jurin's law: 2 * sigma * cos(theta).",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 736,
    "Tool Name": "sort_stresses_descending",
    "Description": "Sorts a list of three stresses in descending order.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 737,
    "Tool Name": "calculate_max_shear_stress",
    "Description": "Calculates the maximum shear stress from principal stresses.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 738,
    "Tool Name": "calculate_plastic_zone_radius",
    "Description": "Calculate the plastic zone radius r_y = K^2 / (2¦Ð¡¤¦Ò_ys^2).",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 739,
    "Tool Name": "calculate_ssa",
    "Description": "Calculate specific surface area (SSA) based on constant k and diameter.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 740,
    "Tool Name": "d_spacing_bcc",
    "Description": "Calculate the interplanar spacing for a BCC crystal.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 741,
    "Tool Name": "is_allowed_reflection_bcc",
    "Description": "Check if the (h,k,l) plane is allowed for diffraction in BCC structure.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 742,
    "Tool Name": "calculate_min_grain_size_from_corrosion_rate",
    "Description": "Calculate minimum grain size according to maximum corrosion rate and porosity.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 743,
    "Tool Name": "calculate_elongation",
    "Description": "Calculate the elongation of the specimen.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 744,
    "Tool Name": "check_allowable_stress",
    "Description": "Determine if interface stress exceeds allowable stress.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 745,
    "Tool Name": "find_zero_strain_conductivity",
    "Description": "Find the conductivity at zero strain.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 746,
    "Tool Name": "compute_characteristic_length",
    "Description": "Calculate the characteristic length as the ratio of diffusion coefficient to front velocity.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 747,
    "Tool Name": "convert_modulus_to_SI",
    "Description": "Convert Young's modulus to Pascals.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 748,
    "Tool Name": "calculate_total_mass",
    "Description": "Calculate the total mass of a material.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 749,
    "Tool Name": "calculate_shear_stress_Pa",
    "Description": "Calculate shear stress in pascals.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 750,
    "Tool Name": "assess_safety",
    "Description": "Assess safety based on shear stress and allowable stress.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 751,
    "Tool Name": "calculate_rectangular_inertia",
    "Description": "Calculate the moment of inertia for a rectangular cross-section.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 752,
    "Tool Name": "calculate_theoretical_packing_factor",
    "Description": "Return the theoretical packing factor for FCC.\n\nReturns:\n    float: Theoretical packing factor (~0.74).",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 753,
    "Tool Name": "calculate_stress_amplitude",
    "Description": "Calculate the stress amplitude from maximum and minimum stresses.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 754,
    "Tool Name": "calculate_fatigue_life_Nf",
    "Description": "Calculate the fatigue life (number of cycles) from stress ratio, fatigue strength coefficient, and exponent.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 755,
    "Tool Name": "determine_elastic_state",
    "Description": "Determine if material is in elastic deformation stage.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 756,
    "Tool Name": "calculate_effective_elastic_modulus",
    "Description": "Calculate effective elastic modulus of porous materials.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 757,
    "Tool Name": "set_material_parameters",
    "Description": "Assign default parameters for a given material type.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 758,
    "Tool Name": "calculate_kd_over_sqrt_d",
    "Description": "Calculate kd divided by sqrt(d).",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 759,
    "Tool Name": "calculate_yield_strength",
    "Description": "Calculate the yield strength based on sigma_0 and kd/¡Ìd.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 760,
    "Tool Name": "calculate_wire_weight",
    "Description": "Calculate the weight of the wire (gravitational force).",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 761,
    "Tool Name": "calculate_moment_of_inertia_solid_sphere",
    "Description": "Calculate the moment of inertia for a solid sphere.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 762,
    "Tool Name": "determine_stronger_material",
    "Description": "Determine which material has higher remaining strength.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 763,
    "Tool Name": "calculate_initial_volume",
    "Description": "Calculate the initial volume of a material.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 764,
    "Tool Name": "calculate_maximum_force",
    "Description": "Calculate the maximum compressive force.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 765,
    "Tool Name": "calculate_modulus_difference",
    "Description": "Calculate the difference between fiber and matrix elastic moduli.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 766,
    "Tool Name": "calculate_fiber_contribution",
    "Description": "Calculate the contribution of fibers to the composite elastic modulus.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 767,
    "Tool Name": "calculate_mass_from_density_volume",
    "Description": "Calculate mass from density and volume.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 768,
    "Tool Name": "calculate_polar_moment_of_inertia",
    "Description": "Calculate the polar moment of inertia J = (¦Ð/2) * r^4.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 769,
    "Tool Name": "calculate_allowable_shear_stress",
    "Description": "Calculate allowable shear stress, which is half of the yield strength.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 770,
    "Tool Name": "determine_plastic_deformation",
    "Description": "Determine if plastic deformation will occur.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 771,
    "Tool Name": "get_material_coefficient",
    "Description": "Retrieve the coefficient for the specified material.",
    "category": "Databases",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 772,
    "Tool Name": "calculate_material_volume",
    "Description": "Calculate the volume of a material based on its mass and density.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 773,
    "Tool Name": "calculate_surface_roughness",
    "Description": "Calculate surface roughness (Ra) based on thickness.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 774,
    "Tool Name": "calculate_interface_stress",
    "Description": "Calculate interface stress.",
    "category": "Computational Tools",
    "Server Name": "Materials_Mechanics_and_Fracture_Analysis"
  },
  {
    "IDX": 775,
    "Tool Name": "convert_resistance_kOhm_to_Ohm",
    "Description": "Convert resistance from kilo-ohms to ohms.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 776,
    "Tool Name": "calculate_voltage",
    "Description": "Calculate required voltage using Ohm's Law V = I ¡Á R.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 777,
    "Tool Name": "calculate_noise_current",
    "Description": "Calculate noise current.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 778,
    "Tool Name": "compute_capacitance_value",
    "Description": "Calculate the capacitance in Farads based on physical parameters.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 779,
    "Tool Name": "calculate_initial_charge",
    "Description": "Calculate the electric charge in a capacitor.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 780,
    "Tool Name": "calculate_new_capacitance",
    "Description": "Calculate the capacitance after inserting dielectric.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 781,
    "Tool Name": "calculate_voltage_from_charge_and_capacitance",
    "Description": "Calculate voltage from charge and capacitance.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 782,
    "Tool Name": "calculate_cell_potential",
    "Description": "Calculate the cell potential (V) from cathode and anode potentials.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 783,
    "Tool Name": "calculate_conductivity",
    "Description": "Compute the electrical conductivity using the Mott VRH model.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 784,
    "Tool Name": "calculate_electric_field",
    "Description": "Calculate the electric field strength between two plates.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 785,
    "Tool Name": "calculate_total_charge",
    "Description": "Calculate the total charge on the oil droplet.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 786,
    "Tool Name": "calculate_mass_from_electric_force",
    "Description": "Calculate the mass of the oil droplet based on electric force balance.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 787,
    "Tool Name": "calculate_charge",
    "Description": "Calculate the charge stored on a capacitor.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 788,
    "Tool Name": "calculate_energy_difference_eV",
    "Description": "Calculate the energy difference between two energy levels (in eV).",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 789,
    "Tool Name": "calculate_resistivity",
    "Description": "Calculate resistivity from conductivity.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 790,
    "Tool Name": "calculate_vacuum_capacitance",
    "Description": "Calculate the vacuum capacitance of a parallel plate capacitor.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 791,
    "Tool Name": "calculate_energy_stored",
    "Description": "Calculate the energy stored in a capacitor.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 792,
    "Tool Name": "calculate_potential_at_position",
    "Description": "Calculate the electric potential at a given position ratio between plates.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 793,
    "Tool Name": "calculate_electric_field_strength",
    "Description": "Calculate the electric field strength produced by a point charge.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 794,
    "Tool Name": "calculate_force_on_charge",
    "Description": "Calculate the force experienced by a charge in an electric field.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 795,
    "Tool Name": "calculate_parallel_resistance",
    "Description": "Calculate the equivalent resistance of resistors in parallel.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 796,
    "Tool Name": "calculate_output_voltage",
    "Description": "Calculate the maximum output voltage of a series-connected battery pack.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 797,
    "Tool Name": "calculate_duty_cycle",
    "Description": "Calculate the duty cycle for a boost converter.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 798,
    "Tool Name": "calculate_induced_emf_magnitude",
    "Description": "Calculate the magnitude of the induced EMF.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 799,
    "Tool Name": "calculate_voltage_across_resistor",
    "Description": "Calculate the voltage across a resistor using Ohm's Law.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 800,
    "Tool Name": "calculate_output_current",
    "Description": "Calculate the current needed for a specified output power given the battery voltage.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 801,
    "Tool Name": "calculate_total_current",
    "Description": "Calculate total current consumption including quiescent and output currents.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 802,
    "Tool Name": "calculate_battery_life_hours",
    "Description": "Calculate the battery life in hours.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 803,
    "Tool Name": "calculate_rc_time_constant",
    "Description": "Calculate the RC time constant (¦Ó) for an RC circuit.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 804,
    "Tool Name": "calculate_absolute_permittivity",
    "Description": "Calculate the absolute permittivity of a dielectric material.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 805,
    "Tool Name": "calculate_capacitance",
    "Description": "Calculate the capacitance of a parallel plate capacitor.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 806,
    "Tool Name": "calculate_stored_energy",
    "Description": "Calculate the energy stored in a capacitor.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 807,
    "Tool Name": "calculate_current_from_heat",
    "Description": "Calculate the current required to generate a specified heat energy in a resistor over a given time.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 808,
    "Tool Name": "calculate_vacuum_permeability",
    "Description": "Calculate the vacuum magnetic permeability (¦Ì?).\n\nReturns:\n    float: Vacuum permeability in T¡¤m/A.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 809,
    "Tool Name": "calculate_magnetic_flux_density",
    "Description": "Calculate the magnetic flux density B.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 810,
    "Tool Name": "calculate_vacuum_permittivity",
    "Description": "Return the vacuum permittivity (epsilon_0) in F/m.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 811,
    "Tool Name": "calculate_charge_coulombs",
    "Description": "Convert energy to charge (Coulombs).",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 812,
    "Tool Name": "convert_coulombs_to_mAh",
    "Description": "Convert Coulombs to milliamp-hours.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 813,
    "Tool Name": "calculate_time_constant",
    "Description": "Calculate the RC time constant.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 814,
    "Tool Name": "calculate_collector_current",
    "Description": "Calculate the collector current in a transistor circuit.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 815,
    "Tool Name": "calculate_emitter_current",
    "Description": "Calculate the emitter current in a transistor circuit.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 816,
    "Tool Name": "calculate_voltage_drop",
    "Description": "Calculate voltage drop across a resistor.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 817,
    "Tool Name": "calculate_power_supply_voltage",
    "Description": "Calculate the power supply voltage Vcc.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 818,
    "Tool Name": "calculate_voltage_formula",
    "Description": "Calculate Vcc using the analytical formula for verification.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 819,
    "Tool Name": "calculate_total_turns",
    "Description": "Calculate the total number of turns in a solenoid.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 820,
    "Tool Name": "calculate_turns_per_meter",
    "Description": "Calculate the number of turns per meter for a solenoid.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 821,
    "Tool Name": "calculate_conduction_electrons",
    "Description": "Calculate the number of conduction electrons.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 822,
    "Tool Name": "calculate_current",
    "Description": "Calculate the current in a circuit.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 823,
    "Tool Name": "calculate_terminal_voltage",
    "Description": "Calculate the terminal voltage of the battery.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 824,
    "Tool Name": "calculate_resistivity_at_temperature",
    "Description": "Calculate resistivity at a specific temperature.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 825,
    "Tool Name": "calculate_thickness_from_resistance",
    "Description": "Calculate the coating thickness in meters.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 826,
    "Tool Name": "calculate_new_short_circuit_current",
    "Description": "Calculate new short-circuit current.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 827,
    "Tool Name": "calculate_new_max_power",
    "Description": "Calculate maximum power output (W).",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 828,
    "Tool Name": "calculate_effective_voltage",
    "Description": "Calculate the effective voltage after voltage drop.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 829,
    "Tool Name": "calculate_voltage_drop_percentage",
    "Description": "Calculate the percentage of voltage drop.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 830,
    "Tool Name": "calculate_average_current",
    "Description": "Calculate the average current over a cycle.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 831,
    "Tool Name": "calculate_average_voltage_drop",
    "Description": "Calculate the average voltage drop over a cycle.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 832,
    "Tool Name": "calculate_velocity",
    "Description": "Calculate the velocity of a charged particle in a uniform magnetic field.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 833,
    "Tool Name": "calculate_conductivity_ratio",
    "Description": "Calculate the conductivity ratio assuming proportionality to mobility.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 834,
    "Tool Name": "calculate_new_conductivity",
    "Description": "Calculate the new conductivity if original conductivity is known.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 835,
    "Tool Name": "calculate_power_consumption",
    "Description": "Calculate the power consumption of a resistor.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 836,
    "Tool Name": "assess_resistor_sufficiency",
    "Description": "Assess whether the resistor's rated power is sufficient and provide recommendation.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 837,
    "Tool Name": "determine_particle_charge",
    "Description": "Determine the electric charge of a particle based on its type.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 838,
    "Tool Name": "calculate_resistance_change_percentage",
    "Description": "Calculate percentage change in resistance.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 839,
    "Tool Name": "calculate_lsb_voltage",
    "Description": "Calculate voltage per LSB for an ADC.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 840,
    "Tool Name": "calculate_inherent_noise_voltage",
    "Description": "Convert inherent noise in LSB to voltage value.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 841,
    "Tool Name": "find_zero_field_jc",
    "Description": "Find the critical current density Jc0 at zero magnetic field.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 842,
    "Tool Name": "find_zero_jc_field",
    "Description": "Find the magnetic field strength H0 at zero critical current density.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 843,
    "Tool Name": "determine_flux_quanta_integer",
    "Description": "Round the exact flux quantum number to an integer using specified method.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 844,
    "Tool Name": "calculate_new_mobility",
    "Description": "Calculate the new mobility after decrease.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 845,
    "Tool Name": "calculate_absolute_decrease",
    "Description": "Calculate the absolute decrease between two values.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 846,
    "Tool Name": "calculate_minimum_field",
    "Description": "Calculate the minimum magnetic field strength needed for magnetization.",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 847,
    "Tool Name": "verify_duty_cycle_range",
    "Description": "Verify that the duty cycle is within the valid range [0, 1].",
    "category": "Computational Tools",
    "Server Name": "Electrical_Engineering_and_Circuit_Calculations"
  },
  {
    "IDX": 848,
    "Tool Name": "calculate_heat_released",
    "Description": "Calculate heat released during combustion.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 849,
    "Tool Name": "calculate_net_buoyancy_coefficient",
    "Description": "Calculate the net buoyancy coefficient from the density ratio.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 850,
    "Tool Name": "calculate_acceleration",
    "Description": "Calculate the upward acceleration of the hot air balloon.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 851,
    "Tool Name": "calculate_ratio",
    "Description": "Calculate the ratio of target phonon mean free path to reference value.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 852,
    "Tool Name": "calculate_new_speed",
    "Description": "Calculate the new welding speed after applying a speed factor.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 853,
    "Tool Name": "calculate_new_power",
    "Description": "Calculate the new laser power to maintain welding depth after increasing speed.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 854,
    "Tool Name": "calculate_theoretical_depth",
    "Description": "Calculate the theoretical welding depth.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 855,
    "Tool Name": "convert_celsius_to_kelvin",
    "Description": "Convert temperature from Celsius to Kelvin.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 856,
    "Tool Name": "calculate_thermal_energy",
    "Description": "Calculate the thermal energy stored during a phase change.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 857,
    "Tool Name": "celsius_to_kelvin",
    "Description": "Convert Celsius temperature to Kelvin.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 858,
    "Tool Name": "get_water_vapor_pressure",
    "Description": "Get water vapor pressure at a specific temperature.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 859,
    "Tool Name": "calculate_thermal_resistance",
    "Description": "Calculate the thermal resistance R = d / k.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 860,
    "Tool Name": "calculate_temperature_change",
    "Description": "Calculate the temperature difference.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 861,
    "Tool Name": "explain_relationship_between_alpha_and_m",
    "Description": "Explain the relationship between thermal expansion coefficient (¦Á) and thermal strain rate sensitivity index (m).",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 862,
    "Tool Name": "aggregate_results",
    "Description": "Aggregate intermediate values and results into a dictionary.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 863,
    "Tool Name": "calculate_heat_absorbed",
    "Description": "Calculate heat absorbed.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 864,
    "Tool Name": "calculate_heat_energy",
    "Description": "Calculate the required heat energy.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 865,
    "Tool Name": "calculate_specific_heat_capacity",
    "Description": "Calculate the specific heat capacity (c) using heat energy, mass, and temperature change.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 866,
    "Tool Name": "get_degeneracy",
    "Description": "Return degeneracy array for each energy level.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 867,
    "Tool Name": "calculate_energy_ratio",
    "Description": "Calculate the ratio of new energy to initial energy.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 868,
    "Tool Name": "calculate_thermal_conductivity_from_ratio",
    "Description": "Calculate thermal conductivity of the target material from a known reference and ratio.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 869,
    "Tool Name": "generate_analytical_expression",
    "Description": "Generate a string for analytical temperature distribution expression.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 870,
    "Tool Name": "calculate_temperature_difference",
    "Description": "Calculate the temperature difference.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 871,
    "Tool Name": "calculate_cooling_time",
    "Description": "Calculate the cooling time.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 872,
    "Tool Name": "compute_temperature_change",
    "Description": "Calculate the temperature change.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 873,
    "Tool Name": "compute_temperature_gradient",
    "Description": "Calculate the temperature gradient (¡ãC/m) along a material.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 874,
    "Tool Name": "calculate_heat_flux",
    "Description": "Calculate the heat flux (W/m^2) using Fourier's law.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 875,
    "Tool Name": "calculate_total_heat_conduction_rate",
    "Description": "Calculate the total heat conduction rate (W) through a material.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 876,
    "Tool Name": "calculate_kinetic_energy",
    "Description": "Calculate kinetic energy from energy density and volume.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 877,
    "Tool Name": "convert_temperature_to_celsius",
    "Description": "Convert normalized temperature to Celsius.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 878,
    "Tool Name": "calculate_temperature_points",
    "Description": "Calculate temperature points T1 and T2 based on T_c and their ratios.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 879,
    "Tool Name": "calculate_energy_released",
    "Description": "Calculate the energy released during phase change.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 880,
    "Tool Name": "calculate_potential_energy",
    "Description": "Calculate the gravitational potential energy.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 881,
    "Tool Name": "calculate_absorbed_heat_energy",
    "Description": "Calculate the heat energy absorbed.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 882,
    "Tool Name": "calculate_temperature_increase",
    "Description": "Calculate the temperature increase.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 883,
    "Tool Name": "calculate_thermal_expansion_difference",
    "Description": "Calculate the difference between the thermal expansion coefficients of two materials.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 884,
    "Tool Name": "convert_to_kelvin",
    "Description": "Convert temperature to Kelvin.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 885,
    "Tool Name": "calculate_carnot_efficiency_value",
    "Description": "Calculate Carnot efficiency.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 886,
    "Tool Name": "calculate_heat_flux_value",
    "Description": "Calculate the heat flux using Fourier's law.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 887,
    "Tool Name": "calculate_heat_flow_rate",
    "Description": "Calculate the total heat flow rate.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 888,
    "Tool Name": "calculate_heat_resistance_sum",
    "Description": "Calculate the total thermal resistance by summing individual resistances.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 889,
    "Tool Name": "include_optional_resistance",
    "Description": "Add optional resistance to total if specified.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 890,
    "Tool Name": "calculate_heat_required",
    "Description": "Calculate the heat energy required for heating.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 891,
    "Tool Name": "calculate_heat_opposite",
    "Description": "Calculate the opposite (negation) of a heat value.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 892,
    "Tool Name": "calculate_melted_mass_from_energy",
    "Description": "Calculate the mass of ice melted from the given energy.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 893,
    "Tool Name": "adjust_energy_for_efficiency",
    "Description": "Adjust energy to account for transfer efficiency.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 894,
    "Tool Name": "calculate_ideal_gas_work",
    "Description": "Calculate the work done during an isothermal process for an ideal gas.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 895,
    "Tool Name": "calculate_internal_energy_change",
    "Description": "Calculate the change in internal energy for an ideal gas during an isothermal process.\n\nReturns:\n    float: Change in internal energy, which is zero for an ideal gas in an isothermal process.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 896,
    "Tool Name": "calculate_heat_transfer",
    "Description": "Calculate the heat transfer based on the first law of thermodynamics.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 897,
    "Tool Name": "convert_energy_MeV_to_J",
    "Description": "Convert energy from MeV to Joules (J).",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 898,
    "Tool Name": "calculate_total_heat",
    "Description": "Calculate total heat based on enthalpy and moles.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 899,
    "Tool Name": "calculate_heat_per_gram",
    "Description": "Calculate heat per gram of adsorbent.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 900,
    "Tool Name": "compute_temperature_difference",
    "Description": "Calculate the temperature difference across a flat plate.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 901,
    "Tool Name": "compute_boltzmann_constant",
    "Description": "Return the Boltzmann constant (J/K).",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 902,
    "Tool Name": "calculate_heat_per_volume",
    "Description": "Calculate heat generation rate per unit volume.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 903,
    "Tool Name": "calculate_ventilation_heat_coefficient",
    "Description": "Calculate heat loss coefficient due to ventilation.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 904,
    "Tool Name": "calculate_final_temperature",
    "Description": "Calculate the final temperature in Kelvin.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 905,
    "Tool Name": "calculate_time_per_mm",
    "Description": "Calculate time (seconds) to weld per millimeter.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 906,
    "Tool Name": "calculate_heat_input",
    "Description": "Calculate heat input in Joules per millimeter.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 907,
    "Tool Name": "compute_clearing_temperature",
    "Description": "Calculate the clearing temperature of a liquid crystal polymer.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 908,
    "Tool Name": "calculate_dynamic_viscosity",
    "Description": "Calculate dynamic viscosity ¦Ì = ¦Í * ¦Ñ.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 909,
    "Tool Name": "calculate_velocity_gradient_at_wall",
    "Description": "Calculate the velocity gradient at y = -h for the flow profile.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 910,
    "Tool Name": "calculate_drag_coefficient",
    "Description": "Calculate the drag coefficient (Cd) based on roughness.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 911,
    "Tool Name": "calculate_drag_force",
    "Description": "Calculate the drag force based on fluid properties and drag coefficient.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 912,
    "Tool Name": "calculate_total_influent_volume",
    "Description": "Calculate total influent volume (m3/m2) over a given period.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 913,
    "Tool Name": "calculate_filter_pore_volume",
    "Description": "Calculate the pore volume of the filter (m3/m2).",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 914,
    "Tool Name": "calculate_fraction_occupied",
    "Description": "Calculate the fraction of filter pore volume occupied by solids.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 915,
    "Tool Name": "calculate_traditional_infiltration",
    "Description": "Calculate infiltration through traditional pavement.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 916,
    "Tool Name": "calculate_permeable_infiltration",
    "Description": "Calculate infiltration through permeable pavement.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 917,
    "Tool Name": "calculate_additional_infiltration",
    "Description": "Calculate the additional infiltration volume.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 918,
    "Tool Name": "calculate_water_mass",
    "Description": "Calculate water mass in grams.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 919,
    "Tool Name": "compute_denominator",
    "Description": "Compute the denominator for Jurin's law: rho * g * d.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 920,
    "Tool Name": "calculate_capillary_height",
    "Description": "Calculate the capillary rise height.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 921,
    "Tool Name": "calculate_liquid_mass",
    "Description": "Calculate the mass of the liquid.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 922,
    "Tool Name": "compute_density_ratio",
    "Description": "Calculate the ratio of two densities.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 923,
    "Tool Name": "calculate_permeability",
    "Description": "Calculate gas permeability.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 924,
    "Tool Name": "calculate_minimum_release_height",
    "Description": "Calculate the minimum release height for a rolling solid sphere on a circular track.",
    "category": "Computational Tools",
    "Server Name": "Thermal _Fluid_Dynamics"
  },
  {
    "IDX": 925,
    "Tool Name": "calculate_incident_photon_rate",
    "Description": "Calculate incident photon rate.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 926,
    "Tool Name": "calculate_power_increase_factor",
    "Description": "Calculate the multiplicative effect of laser power increase on scattered photon rate.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 927,
    "Tool Name": "calculate_reflectance_fresnel",
    "Description": "Calculate reflectance at normal incidence using Fresnel equation.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 928,
    "Tool Name": "compute_wavenumber",
    "Description": "Calculate the wavenumber for a hydrogen spectral line.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 929,
    "Tool Name": "compute_wave_number_ratio",
    "Description": "Calculate the wave number ratio k'/k in a dielectric material.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 930,
    "Tool Name": "calculate_wave_number_in_medium",
    "Description": "Calculate the wave number in a medium.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 931,
    "Tool Name": "calculate_transmission_coefficient",
    "Description": "Calculate the transmission coefficient at the interface.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 932,
    "Tool Name": "calculate_min_wavelength",
    "Description": "Calculate the minimum wavelength based on film thickness.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 933,
    "Tool Name": "calculate_max_frequency",
    "Description": "Calculate the maximum frequency for wave propagation.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 934,
    "Tool Name": "calculate_frequency_range",
    "Description": "Calculate the frequency range of electromagnetic waves.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 935,
    "Tool Name": "calculate_absorption_ratio",
    "Description": "Calculate the absorption ratio based on absorption coefficient and thickness.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 936,
    "Tool Name": "convert_wavelength_nm_to_meters",
    "Description": "Convert wavelength from nanometers to meters.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 937,
    "Tool Name": "generate_wave_functions",
    "Description": "Generate wave function expressions for each quantum number.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 938,
    "Tool Name": "calculate_photon_energy_J",
    "Description": "Calculate photon energy in joules.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 939,
    "Tool Name": "determine_sufficiency",
    "Description": "Determine if photon energy is sufficient to excite electron.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 940,
    "Tool Name": "calculate_numerical_aperture",
    "Description": "Calculate the numerical aperture (NA) of an optical fiber.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 941,
    "Tool Name": "calculate_emission_rate_with_enhancement",
    "Description": "Calculate the emission rate W given free space rate and enhancement.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 942,
    "Tool Name": "calculate_field_change_rate",
    "Description": "Calculate the rate of change of magnetic field.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 943,
    "Tool Name": "electron_wavelength",
    "Description": "Calculate the electron wavelength considering relativistic effects.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 944,
    "Tool Name": "bragg_angle",
    "Description": "Calculate the Bragg angle in radians.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 945,
    "Tool Name": "calculate_intensity",
    "Description": "Calculate laser intensity.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 946,
    "Tool Name": "calculate_radiation_pressure",
    "Description": "Calculate radiation pressure.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 947,
    "Tool Name": "calculate_photon_energy_eV",
    "Description": "Calculate photon energy in eV from wavelength.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 948,
    "Tool Name": "calculate_total_power",
    "Description": "Calculate total power in the spot area.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 949,
    "Tool Name": "calculate_photon_flux",
    "Description": "Calculate the photon flux per second.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 950,
    "Tool Name": "compute_minimum_thickness",
    "Description": "Calculate the minimum film thickness for minimal normal incidence reflection.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 951,
    "Tool Name": "calculate_irradiance_ratio",
    "Description": "Calculate the ratio of new irradiance to standard irradiance.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 952,
    "Tool Name": "calculate_new_mpp_current",
    "Description": "Calculate new max power point current.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 953,
    "Tool Name": "compute_signal_difference",
    "Description": "Compute the difference between light and dark signal means.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 954,
    "Tool Name": "convert_wavelength_to_meters",
    "Description": "Convert wavelength to meters.",
    "category": "Computational Tools",
    "Server Name": "Optics_and_Electromagnetics"
  },
  {
    "IDX": 955,
    "Tool Name": "calculate_volume_ml",
    "Description": "Calculate volume in milliliters from mass and density.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 956,
    "Tool Name": "calculate_molar_mass",
    "Description": "Calculate the molar mass of a compound.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 957,
    "Tool Name": "convert_mass_to_moles",
    "Description": "Convert mass to moles.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 958,
    "Tool Name": "calculate_cell_mass",
    "Description": "Calculate the mass of the unit cell based on density and volume.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 959,
    "Tool Name": "calculate_atoms_per_cell",
    "Description": "Calculate number of atoms per unit cell.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 960,
    "Tool Name": "calculate_volume_in_liters",
    "Description": "Calculate volume in liters.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 961,
    "Tool Name": "calculate_molarity",
    "Description": "Calculate molarity of the solution.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 962,
    "Tool Name": "calculate_pH_from_pOH",
    "Description": "Calculate pH from pOH.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 963,
    "Tool Name": "calculate_molar_mass_ratio",
    "Description": "Calculate the square root of the ratio of two molar masses.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 964,
    "Tool Name": "calculate_solute_mass",
    "Description": "Calculate the mass of solute in a saturated solution.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 965,
    "Tool Name": "calculate_solution_mass",
    "Description": "Calculate total mass of the saturated solution.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 966,
    "Tool Name": "calculate_mass_percent",
    "Description": "Calculate the mass percentage of the solute in the solution.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 967,
    "Tool Name": "calculate_reaction_rate_constant",
    "Description": "Calculate the reaction rate constant k from the half-life.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 968,
    "Tool Name": "calculate_initial_moles",
    "Description": "Calculate initial moles of gas using PV = nRT.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 969,
    "Tool Name": "calculate_molecules_from_moles",
    "Description": "Calculate the number of molecules from moles.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 970,
    "Tool Name": "calculate_remaining_molecules_half_life",
    "Description": "Calculate remaining molecules after reaction_time_hours based on half-life.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 971,
    "Tool Name": "calculate_hydrogen_energy_level",
    "Description": "Calculate the energy of an electron in a hydrogen atom at principal quantum number n.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 972,
    "Tool Name": "calculate_ionization_energy",
    "Description": "Calculate the ionization energy as the absolute value of the energy level.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 973,
    "Tool Name": "calculate_moles_hydrogen",
    "Description": "Calculate moles of hydrogen gas using ideal gas law.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 974,
    "Tool Name": "calculate_moles",
    "Description": "Calculate the amount (moles) of gas.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 975,
    "Tool Name": "calculate_partial_pressure",
    "Description": "Calculate the partial pressure of a gas (atm).",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 976,
    "Tool Name": "calculate_solution_volume",
    "Description": "Calculate solution volume in mL from mass and density.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 977,
    "Tool Name": "determine_electron_transfer_number",
    "Description": "Determine the number of electrons transferred in the reaction.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 978,
    "Tool Name": "calculate_gibbs_free_energy",
    "Description": "Calculate Gibbs free energy change (J/mol).",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 979,
    "Tool Name": "calculate_decay_constant",
    "Description": "Calculate the decay constant ¦Ë given the known final mass, initial mass, and elapsed time.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 980,
    "Tool Name": "calculate_half_life",
    "Description": "Calculate the half-life T?/? from decay constant ¦Ë.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 981,
    "Tool Name": "calculate_remaining_mass",
    "Description": "Calculate the remaining mass after a certain time using decay constant.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 982,
    "Tool Name": "calculate_content_mass",
    "Description": "Calculate the mass of the content substance.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 983,
    "Tool Name": "calculate_moles_from_mol_percent",
    "Description": "Calculate the molar amount of a component given total moles and mol%.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 984,
    "Tool Name": "get_molar_mass",
    "Description": "Retrieve the molar mass of a specified acid.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 985,
    "Tool Name": "convert_density_to_g_per_L",
    "Description": "Convert density from g/mL to g/L.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 986,
    "Tool Name": "calculate_acid_mass_in_solution",
    "Description": "Calculate the mass of acid in 1 liter of solution.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 987,
    "Tool Name": "calculate_mass",
    "Description": "Calculate mass of a component based on weight fraction.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 988,
    "Tool Name": "convert_weight_percentage_to_fraction",
    "Description": "Convert weight percentage to weight fraction.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 989,
    "Tool Name": "compute_final_concentration",
    "Description": "Calculate the final doping concentration after growth.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 990,
    "Tool Name": "get_molecular_mass_amu",
    "Description": "Return molecular mass in amu for a given molecule type.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 991,
    "Tool Name": "calculate_moles_deposited",
    "Description": "Calculate the moles of Ni deposited.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 992,
    "Tool Name": "calculate_mass_deposited",
    "Description": "Calculate the mass of Ni deposited in grams.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 993,
    "Tool Name": "calculate_remaining_moles",
    "Description": "Calculate remaining moles of Ni ions after deposition.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 994,
    "Tool Name": "calculate_final_moles",
    "Description": "Calculate the target moles of Ni ions after treatment.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 995,
    "Tool Name": "calculate_moles_to_adsorb",
    "Description": "Calculate the moles of Ni ions to be adsorbed.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 996,
    "Tool Name": "calculate_adsorbent_mass",
    "Description": "Calculate the mass of adsorbent in kilograms.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 997,
    "Tool Name": "calculate_substance_mass",
    "Description": "Calculate the mass of the solute in grams.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 998,
    "Tool Name": "calculate_pure_water_solubility",
    "Description": "Calculate Mg(OH)? solubility in pure water.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 999,
    "Tool Name": "calculate_naoh_solution_solubility",
    "Description": "Calculate Mg(OH)? solubility in NaOH solution.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1000,
    "Tool Name": "calculate_mass_in_grams",
    "Description": "Calculate mass in grams from volume and density.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1001,
    "Tool Name": "calculate_required_hydroxide_moles",
    "Description": "Calculate the required moles of calcium hydroxide based on acetic acid moles.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1002,
    "Tool Name": "convert_concentration_mgL_to_gm3",
    "Description": "Convert concentration from mg/L to g/m3.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1003,
    "Tool Name": "compute_decay_constant",
    "Description": "Calculate decay constant ¦Ë from half-life.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1004,
    "Tool Name": "calculate_time_from_activity_ratio",
    "Description": "Calculate time t from decay constant and activity ratio.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1005,
    "Tool Name": "calculate_total_moles",
    "Description": "Calculate the total moles of a substance.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1006,
    "Tool Name": "calculate_total_atoms",
    "Description": "Calculate total number of atoms from moles.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1007,
    "Tool Name": "calculate_nanoparticle_count",
    "Description": "Calculate the number of nanoparticles.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1008,
    "Tool Name": "calculate_nanoparticle_moles",
    "Description": "Calculate the total moles of nanoparticles.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1009,
    "Tool Name": "calculate_nanoparticle_concentration",
    "Description": "Calculate the molar concentration of nanoparticles.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1010,
    "Tool Name": "calculate_scaling_factor",
    "Description": "Calculate the scaling factor based on actual and reference PVC masses.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1011,
    "Tool Name": "calculate_additive_mass",
    "Description": "Calculate the additive mass needed based on the scaling factor.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1012,
    "Tool Name": "calculate_co2_moles_from_propane",
    "Description": "Calculate the moles of CO2 produced from propane.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1013,
    "Tool Name": "calculate_bond_order",
    "Description": "Calculate the bond order given bonding and antibonding electrons.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1014,
    "Tool Name": "estimate_bond_length",
    "Description": "Estimate the bond length for a given bond type.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1015,
    "Tool Name": "analyze_bond_length_vs_bond_order",
    "Description": "Generate an analysis statement relating bond order to bond length.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1016,
    "Tool Name": "compute_time_from_alpha",
    "Description": "Calculate the curing time needed to reach a target degree of curing.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1017,
    "Tool Name": "calculate_alpha_at_time",
    "Description": "Calculate the degree of curing at time t.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1018,
    "Tool Name": "calculate_corrosion_rate",
    "Description": "Calculate corrosion rate.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1019,
    "Tool Name": "calculate_cells_per_mole",
    "Description": "Calculate the number of unit cells in one mole of the compound.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1020,
    "Tool Name": "calculate_total_volume_per_mole",
    "Description": "Calculate the total volume occupied by one mole of the compound.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1021,
    "Tool Name": "calculate_molar_mass_from_volume_density",
    "Description": "Calculate the molar mass of the compound.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1022,
    "Tool Name": "calculate_component_mass",
    "Description": "Calculate the mass of a component.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1023,
    "Tool Name": "compute_oxides_masses",
    "Description": "Calculate the mass of each oxide.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1024,
    "Tool Name": "calculate_component_masses",
    "Description": "Calculate zinc and copper masses from total mass and zinc percentage.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1025,
    "Tool Name": "calculate_atoms",
    "Description": "Calculate the number of atoms from moles.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1026,
    "Tool Name": "calculate_total_molecules",
    "Description": "Calculate total gas molecules based on volume and molecules per cc.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1027,
    "Tool Name": "calculate_atomic_radius",
    "Description": "Calculate the atomic radius in FCC structure.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1028,
    "Tool Name": "calculate_two_radii",
    "Description": "Calculate twice the atomic radius.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1029,
    "Tool Name": "calculate_crystallinity_percentage",
    "Description": "Calculate the crystallinity percentage based on density measurements.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1030,
    "Tool Name": "calculate_plating_time",
    "Description": "Calculate the plating time needed to reach the target thickness.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1031,
    "Tool Name": "get_element_magnetic_moment",
    "Description": "Get the magnetic moment value for an element.",
    "category": "Databases",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1032,
    "Tool Name": "calculate_proof",
    "Description": "Calculate alcohol proof from volume percentage.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1033,
    "Tool Name": "calculate_life_material_mass",
    "Description": "Calculate required LiFePO? cathode material mass.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1034,
    "Tool Name": "calculate_number_of_atoms",
    "Description": "Calculate the number of atoms from moles.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1035,
    "Tool Name": "determine_Z_number",
    "Description": "Determine the number of atoms per unit cell based on crystal structure.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1036,
    "Tool Name": "calculate_mass_of_sodium_benzoate",
    "Description": "Calculate the mass of sodium benzoate.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1037,
    "Tool Name": "compute_fukui_plus",
    "Description": "Calculate the nucleophilic Fukui function f?.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1038,
    "Tool Name": "compute_fukui_minus",
    "Description": "Calculate the electrophilic Fukui function f?.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1039,
    "Tool Name": "compute_dual_descriptor",
    "Description": "Calculate the Dual Descriptor f.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1040,
    "Tool Name": "explain_physical_meaning",
    "Description": "Return detailed explanation of the physical meaning of adsorption enthalpy (¦¤H_ads).",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1041,
    "Tool Name": "calculate_total_coordination_number",
    "Description": "Calculate total coordination number based on cyclic ligand denticity and additional ligands.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1042,
    "Tool Name": "determine_typical_lanthanide_range",
    "Description": "Return the typical coordination number range for lanthanide elements.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1043,
    "Tool Name": "assess_coordination_number_range",
    "Description": "Assess whether the total coordination number is within the typical lanthanide range.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1044,
    "Tool Name": "calculate_distance",
    "Description": "Calculate the sum of two ionic radii, representing the distance between ions.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1045,
    "Tool Name": "calculate_tolerance_factor",
    "Description": "Calculate the tolerance factor t.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1046,
    "Tool Name": "assess_stability",
    "Description": "Assess perovskite structure stability based on tolerance factor t.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1047,
    "Tool Name": "parse_electron_configuration",
    "Description": "Parse electron configuration string, extract number of d electrons.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1048,
    "Tool Name": "calculate_expected_electrons",
    "Description": "Calculate expected d electron count based on oxidation state, for verification.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1049,
    "Tool Name": "calculate_spin_S",
    "Description": "Calculate total spin angular momentum S from number of d electrons.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1050,
    "Tool Name": "calculate_standard_cell_potential",
    "Description": "Calculate the standard cell potential (E¡ã) from cathode and anode potentials.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1051,
    "Tool Name": "calculate_reaction_quotient",
    "Description": "Calculate the reaction quotient Q for the cell reaction.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1052,
    "Tool Name": "calculate_single_suppository_mass",
    "Description": "Calculate the mass of a single suppository.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1053,
    "Tool Name": "calculate_gas_constant",
    "Description": "Return gas constant R in J/(mol¡¤K).",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1054,
    "Tool Name": "calculate_faraday_constant",
    "Description": "Return Faraday constant F in C/mol.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1055,
    "Tool Name": "calculate_coefficient_RT_over_nF",
    "Description": "Calculate the coefficient (RT)/(nF) for the Nernst equation.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1056,
    "Tool Name": "generate_molecular_weight_range",
    "Description": "Generate a sequence of molecular weight values within a specified range.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1057,
    "Tool Name": "calculate_total_solids_mass",
    "Description": "Calculate total mass of suspended solids in grams.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1058,
    "Tool Name": "calculate_final_density",
    "Description": "Calculate the final density of a material.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1059,
    "Tool Name": "calculate_total_peg_mass",
    "Description": "Calculate total PEG mass needed for all suppositories.",
    "category": "Computational Tools",
    "Server Name": "Chemistry_and_Reaction_Calculations"
  },
  {
    "IDX": 1060,
    "Tool Name": "calculate_length_plus_width",
    "Description": "Calculate the sum of length and width given the perimeter.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1061,
    "Tool Name": "calculate_length_times_width",
    "Description": "Calculate the product of length and width given the volume and height.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1062,
    "Tool Name": "adjust_dimensions_order",
    "Description": "Ensure that length >= width in the dimensions tuple.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1063,
    "Tool Name": "select_final_dimensions",
    "Description": "Select the first dimension set as the final dimensions.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1064,
    "Tool Name": "round_to_decimal_places",
    "Description": "Round a value to a specified number of decimal places.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1065,
    "Tool Name": "calculate_area",
    "Description": "Calculate the area of a rectangle.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1066,
    "Tool Name": "calculate_theta_rad",
    "Description": "Calculate the angle ¦È in radians from shear strain.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1067,
    "Tool Name": "convert_rad_to_deg",
    "Description": "Convert radians to degrees.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1068,
    "Tool Name": "calculate_phi_deg",
    "Description": "Calculate ¦Õ in degrees from ¦È in degrees.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1069,
    "Tool Name": "calculate_cot_2phi",
    "Description": "Calculate cot(2¦Õ) as tan(¦È).",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1070,
    "Tool Name": "verify_cot_value",
    "Description": "Verify if cot(2¦Õ) is close to shear_strain / 2.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1071,
    "Tool Name": "calculate_volume",
    "Description": "Calculate the volume of a rectangular prism.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1072,
    "Tool Name": "calculate_volume_cm3",
    "Description": "Calculate volume in cubic centimeters.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1073,
    "Tool Name": "round_to_nearest_int",
    "Description": "Round a float to the nearest integer.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1074,
    "Tool Name": "calculate_segment_area",
    "Description": "Calculate the area of a curb segment.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1075,
    "Tool Name": "sum_areas",
    "Description": "Sum a list of area values.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1076,
    "Tool Name": "round_volume_value",
    "Description": "Round the volume in cubic feet to a specified number of decimal places.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1077,
    "Tool Name": "calculate_cross_section_area",
    "Description": "Calculate the cross-sectional area of a rectangular prism.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1078,
    "Tool Name": "calculate_denominator",
    "Description": "Return the denominator (plate separation d).",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1079,
    "Tool Name": "compute_length_change",
    "Description": "Calculate the change in length.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1080,
    "Tool Name": "round_to_significant_figures",
    "Description": "Round a number to a specified number of significant figures.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1081,
    "Tool Name": "symbolic_gamma_mgf_derivation",
    "Description": "Derive the gamma distribution's MGF symbolically using sympy.\n\nReturns:\n    sympy.Expr or str: The symbolic expression of the MGF or a message if derivation is not straightforward.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1082,
    "Tool Name": "round_value",
    "Description": "Round a value to a specified number of decimal places.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1083,
    "Tool Name": "calculate_diameter_from_radius",
    "Description": "Calculate the diameter of a sphere given its radius.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1084,
    "Tool Name": "calculate_atoms_total_volume",
    "Description": "Calculate total volume of all atoms in the unit cell.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1085,
    "Tool Name": "calculate_total_volume",
    "Description": "Calculate total volume from component volumes.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1086,
    "Tool Name": "calculate_volume_fraction",
    "Description": "Calculate volume fraction of a component.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1087,
    "Tool Name": "calculate_outer_volume",
    "Description": "Calculate the volume of a sphere given its radius.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1088,
    "Tool Name": "generate_sphere_mesh",
    "Description": "Generate mesh coordinates for a sphere surface.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1089,
    "Tool Name": "calculate_cross_sectional_area",
    "Description": "Calculate cross-sectional area.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1090,
    "Tool Name": "calculate_increase_factor",
    "Description": "Calculate the increase factor of K-points.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1091,
    "Tool Name": "calculate_percentage_decrease",
    "Description": "Calculate the percentage decrease based on absolute decrease and initial value.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1092,
    "Tool Name": "convert_degrees_to_radians",
    "Description": "Convert angle from degrees to radians.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1093,
    "Tool Name": "convert_percentage_to_ratio",
    "Description": "Convert a percentage to a decimal ratio.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1094,
    "Tool Name": "compute_ratio",
    "Description": "Compute the ratio of numerator to denominator.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1095,
    "Tool Name": "calculate_cylinder_volume",
    "Description": "Calculate the volume of a hollow cylinder.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1096,
    "Tool Name": "calculate_c_analytically",
    "Description": "Compute the constant c analytically so that the integral of f(x) = cx^2 over [0,1] equals 1.\n\nReturns:\n    float: The computed value of c.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1097,
    "Tool Name": "convert_percentage_to_decimal",
    "Description": "Convert a percentage value to decimal.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1098,
    "Tool Name": "calculate_cosine",
    "Description": "Calculate the cosine of an angle in radians.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1099,
    "Tool Name": "calculate_relative_change",
    "Description": "Calculate the ratio of initial to final values.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1100,
    "Tool Name": "calculate_sine",
    "Description": "Calculate the sine of an angle in radians.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1101,
    "Tool Name": "calculate_height_from_length_and_sine",
    "Description": "Calculate the height of the ramp.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1102,
    "Tool Name": "round_to_nearest_tenth",
    "Description": "Round a number to the nearest tenth.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1103,
    "Tool Name": "calculate_absolute_difference",
    "Description": "Calculate the absolute difference between two numbers.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1104,
    "Tool Name": "calculate_bags_needed",
    "Description": "Calculate the number of bags needed to cover an area.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1105,
    "Tool Name": "calculate_total_area",
    "Description": "Calculate the total area of a rectangular region.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1106,
    "Tool Name": "calculate_grass_area",
    "Description": "Calculate the area to be seeded with grass.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1107,
    "Tool Name": "calculate_rectangle_area",
    "Description": "Calculate the area of a rectangle.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1108,
    "Tool Name": "calculate_unit_cell_volume",
    "Description": "Calculate the volume of a unit cell.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1109,
    "Tool Name": "convert_to_percentage",
    "Description": "Convert a ratio to a percentage.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1110,
    "Tool Name": "calculate_percentage_increase",
    "Description": "Calculate the percentage increase from original_value to new_value.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1111,
    "Tool Name": "calculate_max_capacity",
    "Description": "Calculate maximum number of items that can fit into the storage space.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1112,
    "Tool Name": "calculate_max_storage_days",
    "Description": "Calculate maximum storage days based on volume and space utilization.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1113,
    "Tool Name": "calculate_face_diagonal",
    "Description": "Calculate the face diagonal length of the cubic cell.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1114,
    "Tool Name": "calculate_atom_center_distance",
    "Description": "Calculate the distance between neighboring atom centers.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1115,
    "Tool Name": "express_original_expression",
    "Description": "Return the original algebraic expression string for the difference in tank volumes.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1116,
    "Tool Name": "factorize_volume_difference",
    "Description": "Factorize the formula for the volume difference.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1117,
    "Tool Name": "validate_factorization",
    "Description": "Validate the correctness of the factorization.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1118,
    "Tool Name": "parameter_curve_z",
    "Description": "Calculate z-coordinate of the curve at parameter t.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1119,
    "Tool Name": "derivative_z",
    "Description": "Calculate dz/dt at parameter t.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1120,
    "Tool Name": "compute_analytical_probability",
    "Description": "Calculate the analytical probability that X < value.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1121,
    "Tool Name": "degrees_to_radians",
    "Description": "Convert angle from degrees to radians.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1122,
    "Tool Name": "calculate_sqrt_grain_size",
    "Description": "Calculate the square root of the grain size.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1123,
    "Tool Name": "calculate_area_of_square",
    "Description": "Calculate the area of a square.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1124,
    "Tool Name": "calculate_minimum_velocity_at_top",
    "Description": "Calculate the minimum velocity at the top of the track.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1125,
    "Tool Name": "format_volume_output",
    "Description": "Format the volume to a string with two decimal places.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1126,
    "Tool Name": "calculate_orbital_L",
    "Description": "Estimate total orbital angular momentum L, simplified.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1127,
    "Tool Name": "convert_set_to_sorted_list",
    "Description": "Convert a set to a sorted list.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1128,
    "Tool Name": "calculate_surface_area_cube",
    "Description": "Calculate the surface area of a cube-shaped warehouse.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1129,
    "Tool Name": "validate_non_zero_sum_squares",
    "Description": "Validate that the sum of squared Miller indices is not zero.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1130,
    "Tool Name": "compute_sum_of_squares",
    "Description": "Compute the sum of squares of Miller indices.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1131,
    "Tool Name": "symbolic_expectation",
    "Description": "Calculates the symbolic expectation E[X] of the distribution using sympy.\n\nReturns:\n    sympy.Expr: The symbolic expression for E[X].",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1132,
    "Tool Name": "symbolic_second_moment",
    "Description": "Calculates the symbolic second moment E[X^2] using sympy.\n\nReturns:\n    sympy.Expr: The symbolic expression for E[X^2].",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1133,
    "Tool Name": "calculate_initial_area",
    "Description": "Calculate the initial area of the window.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1134,
    "Tool Name": "convert_diameter_to_radius_m",
    "Description": "Convert diameter (millimeter) to radius (meter).",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1135,
    "Tool Name": "convert_diameter_to_radius",
    "Description": "Convert wire diameter to radius.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1136,
    "Tool Name": "calculate_length_contraction",
    "Description": "Calculate the length contraction due to thermal contraction.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1137,
    "Tool Name": "calculate_final_volume",
    "Description": "Calculate the final volume after volume reduction.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1138,
    "Tool Name": "calculate_manufacturing_overhead",
    "Description": "Calculate total manufacturing overhead costs.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1139,
    "Tool Name": "calculate_direct_manufacturing_costs",
    "Description": "Calculate direct manufacturing costs.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1140,
    "Tool Name": "calculate_total_manufacturing_costs",
    "Description": "Calculate total manufacturing costs.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1141,
    "Tool Name": "calculate_non_manufacturing_costs",
    "Description": "Calculate non-manufacturing costs.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1142,
    "Tool Name": "calculate_total_costs",
    "Description": "Calculate total costs.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1143,
    "Tool Name": "calculate_gross_profit",
    "Description": "Calculate gross profit.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1144,
    "Tool Name": "calculate_net_profit",
    "Description": "Calculate net profit.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1145,
    "Tool Name": "calculate_total_input_cost",
    "Description": "Calculate total input cost from labor, material, and overhead costs.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1146,
    "Tool Name": "calculate_productivity_ratio",
    "Description": "Calculate productivity ratio as output value divided by total input cost.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1147,
    "Tool Name": "calculate_explicit_costs",
    "Description": "Calculate explicit costs.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1148,
    "Tool Name": "calculate_opportunity_costs",
    "Description": "Calculate opportunity costs.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1149,
    "Tool Name": "calculate_accounting_profit",
    "Description": "Calculate accounting profit.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1150,
    "Tool Name": "calculate_economic_profit",
    "Description": "Calculate economic profit.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1151,
    "Tool Name": "calculate_total_units",
    "Description": "Calculate total units to account for.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1152,
    "Tool Name": "calculate_completed_units",
    "Description": "Calculate units completed and transferred out.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1153,
    "Tool Name": "calculate_indirect_materials",
    "Description": "Calculate the indirect materials cost.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1154,
    "Tool Name": "calculate_indirect_labor",
    "Description": "Calculate the indirect labor cost.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1155,
    "Tool Name": "calculate_total_direct_cost",
    "Description": "Calculate total direct costs.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1156,
    "Tool Name": "calculate_element_contribution",
    "Description": "Calculate individual element's contribution to total magnetic moment.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1157,
    "Tool Name": "calculate_total_materials_used",
    "Description": "Calculate total materials used by summing direct and indirect materials.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1158,
    "Tool Name": "calculate_ending_inventory",
    "Description": "Calculate ending inventory of raw materials.",
    "category": "Computational Tools",
    "Server Name": "Geometry_and_mathematical_calculations"
  },
  {
    "IDX": 1159,
    "Tool Name": "calculate_absolute_error",
    "Description": "Calculate the absolute error between measured value and true value.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1160,
    "Tool Name": "calculate_percentage_error",
    "Description": "Calculate the percentage error.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1161,
    "Tool Name": "analyze_error_reason",
    "Description": "Analyze possible causes of measurement error.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1162,
    "Tool Name": "convert_to_numpy_array",
    "Description": "Convert input data to a NumPy array.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1163,
    "Tool Name": "calculate_scientific_notation",
    "Description": "Convert a number to scientific notation with specified significant digits.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1164,
    "Tool Name": "format_scientific_notation",
    "Description": "Format mantissa and exponent into a scientific notation string.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1165,
    "Tool Name": "compute_absolute_error",
    "Description": "Calculate the absolute error between measured and true values.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1166,
    "Tool Name": "compute_absolute_error_magnitude",
    "Description": "Calculate the magnitude (absolute value) of an error.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1167,
    "Tool Name": "compute_relative_error_percent",
    "Description": "Calculate the relative error as a percentage.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1168,
    "Tool Name": "format_result",
    "Description": "Round a value to a specified number of decimal places.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1169,
    "Tool Name": "calculate_variance_from_Sxx",
    "Description": "Calculate sample variance from Sxx and sample size.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1170,
    "Tool Name": "format_output",
    "Description": "Format output pore volume and unit.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1171,
    "Tool Name": "calculate_percentage_change",
    "Description": "Calculate percentage change from relative change ratio.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1172,
    "Tool Name": "calculate_percent_error",
    "Description": "Calculate the percentage error based on absolute difference and reference value.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1173,
    "Tool Name": "calculate_mean_square",
    "Description": "Calculate the mean square (MS) from sum of squares (SS) and degrees of freedom (DF).",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1174,
    "Tool Name": "calculate_f_value",
    "Description": "Calculate the F-value for an ANOVA test.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1175,
    "Tool Name": "compute_f_critical_value",
    "Description": "Compute the critical F value for given significance level and degrees of freedom.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1176,
    "Tool Name": "determine_significance",
    "Description": "Determine if the F value indicates significance.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1177,
    "Tool Name": "calculate_max_value",
    "Description": "Calculate the maximum integer value for a given bit resolution.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1178,
    "Tool Name": "convert_humidity_to_percent",
    "Description": "Convert normalized humidity to percentage.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1179,
    "Tool Name": "calculate_relative_uncertainty",
    "Description": "Calculate the relative uncertainty.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1180,
    "Tool Name": "calculate_combined_relative_uncertainty",
    "Description": "Calculate the combined relative uncertainty using root-sum-square.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1181,
    "Tool Name": "calculate_relative_error",
    "Description": "Calculate the relative error as a ratio.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1182,
    "Tool Name": "check_data_length",
    "Description": "Check if two arrays have the same length.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1183,
    "Tool Name": "print_parameters",
    "Description": "Print the input parameters with formatting for clarity.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1184,
    "Tool Name": "calculate_additional_space_needed",
    "Description": "Calculate extra space needed if storage is insufficient.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1185,
    "Tool Name": "convert_to_scientific_notation",
    "Description": "Convert a ratio to scientific notation with mantissa rounded to one decimal.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1186,
    "Tool Name": "pdf",
    "Description": "Calculate the probability density function at point(s) x.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1187,
    "Tool Name": "ensure_array",
    "Description": "Convert input to a numpy array.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1188,
    "Tool Name": "generate_q_range",
    "Description": "Generate a list of integer q values from start to end inclusive.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1189,
    "Tool Name": "calculate_theoretical_variance",
    "Description": "Calculate the theoretical variance of a normal distribution.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1190,
    "Tool Name": "calculate_median",
    "Description": "Calculate the median of a list of numerical data.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1191,
    "Tool Name": "calculate_mode",
    "Description": "Calculate the mode(s) of a list of numerical data.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1192,
    "Tool Name": "convert_to_set",
    "Description": "Convert an input iterable to a set.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1193,
    "Tool Name": "calculate_max_quantization_error",
    "Description": "Calculate the maximum quantization error, equal to half the LSB voltage.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1194,
    "Tool Name": "compute_variance_difference",
    "Description": "Compute the difference between light and dark signal variances.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1195,
    "Tool Name": "compute_standard_deviation",
    "Description": "Calculate the standard deviation from the variance difference.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1196,
    "Tool Name": "validate_interval_vs_total",
    "Description": "Validate that interval individuals do not exceed total individuals.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1197,
    "Tool Name": "compute_relative_density",
    "Description": "Compute interval relative density.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1198,
    "Tool Name": "calculate_probability",
    "Description": "Calculate probability that device life exceeds threshold.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1199,
    "Tool Name": "generate_collision_time_meshgrid",
    "Description": "Create mesh grids for electron density and temperature ranges.",
    "category": "Computational Tools",
    "Server Name": "Data_processing_and_statistical_analysis"
  },
  {
    "IDX": 1200,
    "Tool Name": "convert_length_mm_to_m",
    "Description": "Convert length from millimeters to meters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1201,
    "Tool Name": "convert_ml_to_dl",
    "Description": "Convert volume from milliliters to deciliters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1202,
    "Tool Name": "convert_current_mA_to_A",
    "Description": "Convert current from milliamperes to amperes.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1203,
    "Tool Name": "convert_pm_to_cm",
    "Description": "Convert length from picometers to centimeters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1204,
    "Tool Name": "convert_capacitance_unit",
    "Description": "Convert capacitance from Farads to specified units.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1205,
    "Tool Name": "convert_kg_to_ton",
    "Description": "Convert weight from kilograms to tons.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1206,
    "Tool Name": "convert_mass_kg_to_g",
    "Description": "Convert mass from kilograms to grams.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1207,
    "Tool Name": "convert_pressure_mmHg_to_atm",
    "Description": "Convert pressure from mm mercury to atmospheres.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1208,
    "Tool Name": "convert_cubic_meters_to_cubic_feet",
    "Description": "Convert volume from cubic meters to cubic feet.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1209,
    "Tool Name": "convert_length_to_meters",
    "Description": "Convert length from kilometers to meters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1210,
    "Tool Name": "convert_cm_to_meters",
    "Description": "Convert length from centimeters to meters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1211,
    "Tool Name": "convert_mmHg_to_atm",
    "Description": "Convert pressure from mm Hg to atm.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1212,
    "Tool Name": "convert_mass_to_grams",
    "Description": "Convert mass to grams.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1213,
    "Tool Name": "convert_volume_to_mL",
    "Description": "Convert volume to milliliters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1214,
    "Tool Name": "convert_thickness_nm_to_m",
    "Description": "Convert thickness from nanometers to meters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1215,
    "Tool Name": "convert_thickness_to_meters",
    "Description": "Convert film thickness to meters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1216,
    "Tool Name": "convert_area_cm2_to_m2",
    "Description": "Convert area from square centimeters to square meters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1217,
    "Tool Name": "convert_mm2_to_m2",
    "Description": "Convert cross-sectional area from mm2 to m2.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1218,
    "Tool Name": "convert_mm_to_meters",
    "Description": "Convert length from millimeters to meters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1219,
    "Tool Name": "convert_kN_to_N",
    "Description": "Convert force from kilonewtons to newtons.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1220,
    "Tool Name": "convert_kg_to_g",
    "Description": "Convert mass from kilograms to grams.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1221,
    "Tool Name": "convert_grain_size_to_meters",
    "Description": "Convert grain size from micrometers to meters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1222,
    "Tool Name": "convert_nm_to_m",
    "Description": "Convert length from nanometers to meters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1223,
    "Tool Name": "convert_J_to_eV",
    "Description": "Convert energy from Joules to electron volts.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1224,
    "Tool Name": "convert_area_to_square_meters",
    "Description": "Convert area to square meters based on unit system.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1225,
    "Tool Name": "convert_mass_g_to_kg",
    "Description": "Convert mass from grams to kilograms.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1226,
    "Tool Name": "convert_volume_L_to_m3",
    "Description": "Convert volume from liters to cubic meters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1227,
    "Tool Name": "convert_time_minutes_to_seconds",
    "Description": "Convert time from minutes to seconds.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1228,
    "Tool Name": "convert_nm_to_um",
    "Description": "Convert thickness from nanometers to micrometers.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1229,
    "Tool Name": "convert_force_pN_to_N",
    "Description": "Convert force from picoNewtons to Newtons.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1230,
    "Tool Name": "convert_distance_to_meters",
    "Description": "Convert distance to meters if input is in centimeters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1231,
    "Tool Name": "convert_length_to_mm",
    "Description": "Convert length change from meters to millimeters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1232,
    "Tool Name": "convert_velocity_kmh_to_ms",
    "Description": "Convert velocity from kilometers per hour to meters per second.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1233,
    "Tool Name": "convert_amu_to_kg",
    "Description": "Convert atomic mass units to kilograms.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1234,
    "Tool Name": "convert_liters_to_milliliters",
    "Description": "Convert volume from liters to milliliters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1235,
    "Tool Name": "convert_grams_to_kilograms",
    "Description": "Convert mass from grams to kilograms.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1236,
    "Tool Name": "convert_grams_to_pounds",
    "Description": "Convert mass from grams to pounds.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1237,
    "Tool Name": "convert_work_function_to_joules",
    "Description": "Convert work function to Joules.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1238,
    "Tool Name": "format_density_unit",
    "Description": "Format the density unit string.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1239,
    "Tool Name": "convert_torr_to_atm",
    "Description": "Convert pressure from Torr to atm.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1240,
    "Tool Name": "convert_energy_to_joules",
    "Description": "Convert energy to Joules.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1241,
    "Tool Name": "convert_dimensions_to_meters",
    "Description": "Convert dimension to meters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1242,
    "Tool Name": "convert_hours_to_days",
    "Description": "Convert time in hours to days.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1243,
    "Tool Name": "convert_radius_to_meters",
    "Description": "Convert radius from centimeters to meters if necessary.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1244,
    "Tool Name": "convert_inches_to_centimeters",
    "Description": "Convert length from inches to centimeters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1245,
    "Tool Name": "convert_meters_to_centimeters",
    "Description": "Convert length from meters to centimeters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1246,
    "Tool Name": "convert_length_m_to_mm",
    "Description": "Convert length from meters to millimeters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1247,
    "Tool Name": "convert_ml_to_m3",
    "Description": "Convert volume from milliliters to cubic meters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1248,
    "Tool Name": "convert_mH_to_H",
    "Description": "Convert inductance from millihenries to henries.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1249,
    "Tool Name": "convert_h2_production_rate_to_molecules_per_s",
    "Description": "Convert H2 production rate from mmol/cm2/h to molecules/cm2/s.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1250,
    "Tool Name": "convert_length_units",
    "Description": "Convert length from meters to specified units.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1251,
    "Tool Name": "convert_length_to_SI",
    "Description": "Convert length to meters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1252,
    "Tool Name": "convert_area_to_SI",
    "Description": "Convert area to square meters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1253,
    "Tool Name": "convert_length_from_SI",
    "Description": "Convert length from meters to target units.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1254,
    "Tool Name": "convert_work_function_to_J",
    "Description": "Convert work function from electron volts to joules.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1255,
    "Tool Name": "convert_force_kN_to_N",
    "Description": "Convert force from kilonewtons to newtons.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1256,
    "Tool Name": "convert_liters_to_cubic_centimeters",
    "Description": "Convert volume from liters to cubic centimeters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1257,
    "Tool Name": "convert_thickness_mm_to_um",
    "Description": "Convert thickness from millimeters to micrometers.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1258,
    "Tool Name": "convert_modulus_to_mpa",
    "Description": "Convert elastic modulus from GPa to MPa.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1259,
    "Tool Name": "convert_energy_units",
    "Description": "Convert energy units based on a conversion factor.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1260,
    "Tool Name": "convert_speed_kmh_to_ms",
    "Description": "Convert speed from km/h to m/s.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1261,
    "Tool Name": "convert_thickness_km_to_m",
    "Description": "Convert thickness from kilometers to meters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1262,
    "Tool Name": "convert_concentration_cm3_to_m3",
    "Description": "Convert carrier concentration from cm^-3 to m^-3.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1263,
    "Tool Name": "convert_mobility_cm2Vs_to_m2Vs",
    "Description": "Convert mobility from cm^2/(V¡¤s) to m^2/(V¡¤s).",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1264,
    "Tool Name": "convert_volume_to_cm3",
    "Description": "Convert volume to cubic centimeters based on the input unit.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1265,
    "Tool Name": "convert_density_to_SI",
    "Description": "Convert density from g/cm3 to g/m3.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1266,
    "Tool Name": "convert_volume_liters_to_cm3",
    "Description": "Convert volume from liters to cubic centimeters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1267,
    "Tool Name": "convert_meters_to_micrometers",
    "Description": "Convert length from meters to micrometers.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1268,
    "Tool Name": "convert_gallons_to_liters",
    "Description": "Convert gallons to liters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1269,
    "Tool Name": "convert_milliliters_to_cubic_centimeters",
    "Description": "Convert milliliters to cubic centimeters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1270,
    "Tool Name": "convert_activity_uCi_to_Bq",
    "Description": "Convert microcuries (¦ÌCi) to becquerel (Bq).",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1271,
    "Tool Name": "calculate_max_power_kw",
    "Description": "Convert maximum power from watts to kilowatts.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1272,
    "Tool Name": "convert_ev_to_joules",
    "Description": "Convert energy from electron volts to joules.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1273,
    "Tool Name": "convert_radius_m_to_cm",
    "Description": "Convert radius from meters to centimeters.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1274,
    "Tool Name": "convert_distance_mm_to_m",
    "Description": "Convert distance from mm to m.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1275,
    "Tool Name": "convert_pressure_kpa_to_pa",
    "Description": "Convert pressure from kPa to Pa.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1276,
    "Tool Name": "convert_bulk_modulus_kpa_to_pa",
    "Description": "Convert bulk modulus from kPa to Pa.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1277,
    "Tool Name": "convert_stress_pa_to_mpa",
    "Description": "Convert stress from pascal (Pa) to megapascal (MPa).",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1278,
    "Tool Name": "convert_speed_to_mm_per_s",
    "Description": "Convert welding speed from m/min to mm/s.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1279,
    "Tool Name": "convert_cm_to_mm",
    "Description": "Convert length from centimeters to millimeters, rounded to 2 decimal places.",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1280,
    "Tool Name": "get_physical_constants",
    "Description": "Return a dictionary of fundamental physical constants.\n\nReturns:\n    dict: A dictionary with keys 'h' (Planck's constant) and 'c' (speed of light).",
    "category": "Computational Tools",
    "Server Name": "Physical_Quantities_Conversion"
  },
  {
    "IDX": 1281,
    "Tool Name": "ChemicalStructureAnalyzer",
    "Description": "Complete structure analysis from compound name (SMILES + properties).",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1282,
    "Tool Name": "MolecularWeightCalculator",
    "Description": "Calculate molecular weight from compound name.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1283,
    "Tool Name": "FunctionalGroupAnalyzer",
    "Description": "Identify all functional groups in a molecule.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1284,
    "Tool Name": "MolecularDescriptorCalculator",
    "Description": "Calculate comprehensive molecular descriptors.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1285,
    "Tool Name": "LipinskiRuleChecker",
    "Description": "Check Lipinski's Rule of Five for drug-likeness.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1286,
    "Tool Name": "MolecularSimilarityComparator",
    "Description": "Compare two molecules for similarity using Tanimoto coefficient.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1287,
    "Tool Name": "FingerprintGenerator",
    "Description": "Generate multiple types of molecular fingerprints.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1288,
    "Tool Name": "StereochemistryAnalyzer",
    "Description": "Analyze stereochemical properties.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1289,
    "Tool Name": "RingSystemAnalyzer",
    "Description": "Analyze ring systems in molecules.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1290,
    "Tool Name": "Conformation3DAnalyzer",
    "Description": "Analyze 3D conformational properties.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1291,
    "Tool Name": "ElectronicPropertyCalculator",
    "Description": "Calculate electronic properties of molecules.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1292,
    "Tool Name": "MolecularShapeDescriptor",
    "Description": "Calculate molecular shape descriptors.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1293,
    "Tool Name": "ChiralityAnalyzer",
    "Description": "Comprehensive chirality analysis.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1294,
    "Tool Name": "BondAnalyzer",
    "Description": "Analyze bond properties in molecules.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1295,
    "Tool Name": "HeteroatomAnalyzer",
    "Description": "Analyze heteroatom content.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1296,
    "Tool Name": "FormalChargeCalculator",
    "Description": "Calculate formal charges.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1297,
    "Tool Name": "FragmentationAnalyzer",
    "Description": "Analyze molecular fragments.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1298,
    "Tool Name": "TopologicalIndexCalculator",
    "Description": "Calculate topological indices.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1299,
    "Tool Name": "KappaIndicesCalculator",
    "Description": "Calculate molecular kappa shape indices.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1300,
    "Tool Name": "AromaticityAnalyzer",
    "Description": "Analyze aromatic systems.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1301,
    "Tool Name": "StructureFormatConverter",
    "Description": "Convert between molecular structure formats.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1302,
    "Tool Name": "NameToAllFormats",
    "Description": "Convert compound name to all structure formats.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1303,
    "Tool Name": "InChIKeyResolver",
    "Description": "Resolve InChIKey to other formats.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1304,
    "Tool Name": "SELFIESConverter",
    "Description": "Bidirectional SELFIES conversion.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1305,
    "Tool Name": "CASNumberLookup",
    "Description": "Lookup CAS number for a compound.",
    "category": "Databases",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1306,
    "Tool Name": "StructureValidator",
    "Description": "Validate and standardize molecular structures.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1307,
    "Tool Name": "TautomerGenerator",
    "Description": "Generate tautomers of a molecule.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1308,
    "Tool Name": "InChIValidator",
    "Description": "Validate InChI and InChIKey strings.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1309,
    "Tool Name": "MoleculeStandardizer",
    "Description": "Standardize molecular structures.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1310,
    "Tool Name": "SubstructureSearcher",
    "Description": "Search for substructures in molecules.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1311,
    "Tool Name": "ProteinPropertyCalculator",
    "Description": "Calculate comprehensive protein properties.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1312,
    "Tool Name": "ProteinStabilityAnalyzer",
    "Description": "Analyze protein stability indicators.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1313,
    "Tool Name": "SequenceAlignmentAnalyzer",
    "Description": "Perform sequence alignment analysis.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1314,
    "Tool Name": "ProteinMotifFinder",
    "Description": "Find common protein motifs.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1315,
    "Tool Name": "ProteinSolubilityPredictor",
    "Description": "Predict protein solubility.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1316,
    "Tool Name": "AntibodyAnalyzer",
    "Description": "Analyze antibody sequence features.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1317,
    "Tool Name": "ProteinLocalizationPredictor",
    "Description": "Predict subcellular localization.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1318,
    "Tool Name": "AminoAcidCompositionAnalyzer",
    "Description": "Analyze amino acid composition.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1319,
    "Tool Name": "ProteinInteractionPredictor",
    "Description": "Predict protein-protein interaction potential.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1320,
    "Tool Name": "DrugTargetInteractionPredictor",
    "Description": "Predict drug-target interaction.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1321,
    "Tool Name": "DNASequenceAnalyzer",
    "Description": "Comprehensive DNA sequence analysis.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1322,
    "Tool Name": "GeneticCodeTranslator",
    "Description": "Translate DNA to protein.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1323,
    "Tool Name": "DNAComplementFinder",
    "Description": "Find DNA complement and reverse complement.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1324,
    "Tool Name": "PalindromeFinder",
    "Description": "Find palindromic sequences.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1325,
    "Tool Name": "CodonOptimizer",
    "Description": "Optimize codon usage.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1326,
    "Tool Name": "DNARNACodonOptimizer",
    "Description": "Optimize DNA/RNA codons.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1327,
    "Tool Name": "PCRPrimerDesigner",
    "Description": "Design PCR primers.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1328,
    "Tool Name": "RestrictionSiteAnalyzer",
    "Description": "Analyze restriction enzyme sites.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1329,
    "Tool Name": "CircularDNAAnalyzer",
    "Description": "Analyze circular DNA.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1330,
    "Tool Name": "RandomDNAGenerator",
    "Description": "Generate random DNA sequence.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1331,
    "Tool Name": "PeptidePropertyCalculator",
    "Description": "Calculate peptide properties.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1332,
    "Tool Name": "AlanineScanningDesigner",
    "Description": "Design alanine scanning library.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1333,
    "Tool Name": "TruncationLibraryDesigner",
    "Description": "Design truncation library.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1334,
    "Tool Name": "OverlapPeptideDesigner",
    "Description": "Design overlapping peptide library.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1335,
    "Tool Name": "PositionalScanningDesigner",
    "Description": "Design positional scanning library.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1336,
    "Tool Name": "ProteaseDigestionAnalyzer",
    "Description": "Analyze protease digestion patterns.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1337,
    "Tool Name": "DegenerateCodonCalculator",
    "Description": "Calculate degenerate codons.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1338,
    "Tool Name": "OligonucleotideCalculator",
    "Description": "Calculate oligonucleotide properties.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1339,
    "Tool Name": "PeptideToSMILESConverter",
    "Description": "Convert peptide sequence to SMILES.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1340,
    "Tool Name": "ComputeAffinityCalculator",
    "Description": "Compute binding affinity.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1341,
    "Tool Name": "ComprehensiveADMETPredictor",
    "Description": "Predict ADMET properties for compounds.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1342,
    "Tool Name": "DrugLikenessAnalyzer",
    "Description": "Analyze drug-likeness properties.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1343,
    "Tool Name": "CompoundToADMET",
    "Description": "From compound name to ADMET prediction.",
    "category": "Model Services",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1344,
    "Tool Name": "BBBPenetrancePredictor",
    "Description": "Predict blood-brain barrier penetrance.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1345,
    "Tool Name": "BioavailabilityPredictor",
    "Description": "Predict oral bioavailability.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1346,
    "Tool Name": "PubChemCompoundLookup",
    "Description": "Complete PubChem compound information.",
    "category": "Databases",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1347,
    "Tool Name": "CompoundSimilaritySearch",
    "Description": "Search for similar compounds.",
    "category": "Databases",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1348,
    "Tool Name": "CompoundPropertiesRetriever",
    "Description": "Retrieve compound properties from PubChem.",
    "category": "Databases",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1349,
    "Tool Name": "CompoundNameResolver",
    "Description": "Resolve compound name to all identifiers.",
    "category": "Databases",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1350,
    "Tool Name": "CIDToPropertiesConverter",
    "Description": "Convert CID to comprehensive properties.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1351,
    "Tool Name": "MolecularOptimizer",
    "Description": "Optimize molecule for drug-likeness.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1352,
    "Tool Name": "LeadOptimizationAnalyzer",
    "Description": "Analyze lead compound for optimization.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1353,
    "Tool Name": "ScaffoldAnalyzer",
    "Description": "Analyze molecular scaffold.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1354,
    "Tool Name": "RandomMoleculeGenerator",
    "Description": "Generate random molecules.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1355,
    "Tool Name": "TautomerEnumerator",
    "Description": "Enumerate all tautomers.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1356,
    "Tool Name": "ProteinLigandInteractionPredictor",
    "Description": "Predict protein-ligand interactions.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1357,
    "Tool Name": "UniProtProteinAnalyzer",
    "Description": "Analyze protein from UniProt.",
    "category": "Databases",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1358,
    "Tool Name": "TargetProteinProfiler",
    "Description": "Complete target protein profile.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1359,
    "Tool Name": "DrugTargetValidator",
    "Description": "Validate drug-target pair.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1360,
    "Tool Name": "SmallMoleculeAffinityCalculator",
    "Description": "Calculate small molecule similarity.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1361,
    "Tool Name": "MOFStructureAnalyzer",
    "Description": "Analyze MOF structure.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1362,
    "Tool Name": "MaterialDensityCalculator",
    "Description": "Calculate material density.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1363,
    "Tool Name": "MaterialSymmetryAnalyzer",
    "Description": "Analyze material symmetry.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1364,
    "Tool Name": "MaterialCompositionAnalyzer",
    "Description": "Analyze element composition.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1365,
    "Tool Name": "MaterialStructureInfoRetriever",
    "Description": "Get structure information.",
    "category": "Databases",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1366,
    "Tool Name": "MOFToCompoundConverter",
    "Description": "Convert MOF to SMILES.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1367,
    "Tool Name": "CompoundToMaterialPrice",
    "Description": "Get material price from compound.",
    "category": "Databases",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1368,
    "Tool Name": "MaterialLatticeAnalyzer",
    "Description": "Analyze MOF lattice parameters.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1369,
    "Tool Name": "ComprehensiveMaterialAnalyzer",
    "Description": "Complete material analysis.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1370,
    "Tool Name": "SMILESToCASConverter",
    "Description": "Convert SMILES to CAS number.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1371,
    "Tool Name": "SmallMoleculeToProteinInteraction",
    "Description": "Analyze small molecule-protein interaction.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1372,
    "Tool Name": "PeptideDrugDesigner",
    "Description": "Design peptide-based drugs.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1373,
    "Tool Name": "BioactiveCompoundScreener",
    "Description": "Screen bioactive compounds.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1374,
    "Tool Name": "ProteinSMILESConverter",
    "Description": "Convert peptide to SMILES and analyze.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1375,
    "Tool Name": "DrugPeptideOptimizer",
    "Description": "Optimize peptide for drug properties.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1376,
    "Tool Name": "CompoundToADMETProfile",
    "Description": "Complete compound ADMET profile.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1377,
    "Tool Name": "DrugCandidateScreener",
    "Description": "Screen drug candidates.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1378,
    "Tool Name": "LeadCompoundProfiler",
    "Description": "Profile lead compound.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1379,
    "Tool Name": "CompoundDrugLikenessScorer",
    "Description": "Score compound drug-likeness.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1380,
    "Tool Name": "DrugSimilarityAnalyzer",
    "Description": "Analyze drug similarity.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1381,
    "Tool Name": "ProteinTargetDrugDesigner",
    "Description": "Design drugs for protein target.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1382,
    "Tool Name": "AntibodyDrugConjugateDesigner",
    "Description": "Design antibody-drug conjugates.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1383,
    "Tool Name": "BiologicDrugAnalyzer",
    "Description": "Analyze biologic drugs.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1384,
    "Tool Name": "TargetSequenceValidator",
    "Description": "Validate target sequence for druggability.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1385,
    "Tool Name": "ProteinDrugInteractionProfiler",
    "Description": "Profile protein-drug interactions.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1386,
    "Tool Name": "MOFCompoundAnalyzer",
    "Description": "Analyze MOF as compound.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1387,
    "Tool Name": "MaterialCompoundConverter",
    "Description": "Convert material to compound format.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1388,
    "Tool Name": "MOFPropertiesCalculator",
    "Description": "Calculate MOF properties.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1389,
    "Tool Name": "CrystalStructureAnalyzer",
    "Description": "Analyze crystal structure.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1390,
    "Tool Name": "MaterialDrugDeliveryAnalyzer",
    "Description": "Analyze material for drug delivery.",
    "category": "Computational Tools",
    "Server Name": "InternAgent"
  },
  {
    "IDX": 1391,
    "Tool Name": "ComputeProtPara",
    "Description": "Compute various physical and chemical parameters for a given protein sequence using Expasy ProtParam API.\n    Parameters include molecular weight, theoretical pI, amino acid composition, atomic composition, \n    extinction coefficient, half-life, instability index, aliphatic index, and GRAVY.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1392,
    "Tool Name": "ComputeProtScale",
    "Description": "Predict the hydrophilicity of a protein sequence.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1393,
    "Tool Name": "ComputeExtinctionCoefficient",
    "Description": "This tool compute the molar extinction coefficient and protein concentration of the protein, \n    and also provides information such as the protein isoelectric point. \n\n    Returns:\n        str: The Markdown content with molar extinction coefficient.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1394,
    "Tool Name": "ComputePiMw",
    "Description": "Compute the theoretical isoelectric point (pI) and molecular weight (mW) of a protein sequence.\n    The input should be a protein sequence.\n\n    Returns:\n        str: A string containing the calculated pI and mW.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1395,
    "Tool Name": "CipherOptimizer",
    "Description": "Codon Optimization Tool: Used to optimize codons for expression of recombinant genes in mainstream hosts. The parameters optimized include up to a dozen key parameters for both transcription and translation processes.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1396,
    "Tool Name": "CalculatorPeptideProperty",
    "Description": "Peptide Property Calculator: Calculate molecular weight, extinction coefficient, net peptide charge, peptide isoelectric point, and average hydrophobicity (GRAVY) of peptide properties.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1397,
    "Tool Name": "CalculatorOligonucleotide",
    "Description": "Oligonucleotide (primer) Calculator: The annealing temperature (Tm), molecular weight (MW), extinction coefficient (OD/¦Ìmol, ¦Ìg/OD) of the oligonucleotides were calculated.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1398,
    "Tool Name": "ProteinCodonOptimization",
    "Description": "This tool optimize the the expression of recombinant gene condons of the protein,\n\n    Input:\n        protein: protein sequence\n\n    Returns:\n        str: The Markdown content with new sequence.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1399,
    "Tool Name": "DNARNACodonOptimization",
    "Description": "This tool optimize the the expression of recombinant gene condons of the DNA/RNA,\n\n    Input:\n        protein: DNA or RNA sequence. The sequence length must be a multiple of 3¡£\n\n    Returns:\n        str: The Markdown content with new sequence.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1400,
    "Tool Name": "ComputeHydrophilicity",
    "Description": "This tool compute the hydrophilicity of the protein,",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1401,
    "Tool Name": "ComputeAnnealingTemperature",
    "Description": "This tool compute the annealing temperature of an oligonucleotide.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1402,
    "Tool Name": "ConvertingPeptide2SMILES",
    "Description": "This tool translate the polypeptide sequence to SMILES.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1403,
    "Tool Name": "ProteinIsoelectricPointCalculator",
    "Description": "This tool compute the isoelectric point of the protein or peptide.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1404,
    "Tool Name": "ComputeAffinity",
    "Description": "This tool compute affinity based on the molar Gibbs free energy",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1405,
    "Tool Name": "PeptideWeightCalculator",
    "Description": "This tool compute the Average molecular weight of the polypeptide.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1406,
    "Tool Name": "PeptideFormulaCalculator",
    "Description": "This tool compute the chemical formula of the polypeptide.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1407,
    "Tool Name": "DegenerateCodonCalculatorbyAminoAcid",
    "Description": "This tool calculates the optimal degenerate codons that encode one or more input amino acids. It can be applied to library construction.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1408,
    "Tool Name": "OverlapPeptideLibraryDesign",
    "Description": "This tool design overlapping peptide library.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1409,
    "Tool Name": "AlanineScanningLibraryDesign",
    "Description": "This tool design peptide library.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1410,
    "Tool Name": "TruncationLibraryDesign",
    "Description": "This tool design truncation peptide library.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1411,
    "Tool Name": "PositionalScanningLibraryDesign",
    "Description": "This tool design positional scanning peptide library.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1412,
    "Tool Name": "ProteaseDigestion",
    "Description": "This tool can simulate the hydrolytic behavior of protein-degrading enzymes. Its main purpose is to predict the hydrolysis outcomes of peptide substrates.\"",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1413,
    "Tool Name": "CDRLabelingAntibody",
    "Description": "This tool label the variable regions of antibodies with CDR and FR regions; Users need to choose a numbering system, and the numbering schemes include: imgt, chothia, kabat, martin; The definition scheme includes chothia, kabat, imgt, and contact. Note: The Kabat numbering scheme is not compatible with the Contact definition scheme. To number the amino acid sequence of antibodies, please visit the antibody sequence numbering tool\n    You should provide 3 input arguments.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1414,
    "Tool Name": "AntibodySequenceNumbering",
    "Description": "This tool number the amino acid sequence of the antibody; Identify the input sequence and distinguish between immunoglobulin (IG) and T cell receptor (TR); The numbering system includes: IMGT, Chothia, Kabat, Martin (extended version Chothia), and AHo; TR sequences can only be numbered using IMGT or AHo. To label the CDR and FR of antibody variable regions, please use the antibody variable region CDR labeling tool",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1415,
    "Tool Name": "CircularDNAAlignment",
    "Description": "This tool aligns circular DNA sequences. Comparing homologous sequences of different species is an important method for reconstructing the evolutionary history of genes in species and their genomes. The circular DNA sequence alignment tool can be applied to sequence alignment of circular nucleic acid molecules, including plasmids, mitochondrial DNA (mtDNA), circular bacterial chromosomes, cccDNA (covalently closed circular DNA), chloroplast DNA (cpDNA), and other plastids.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1416,
    "Tool Name": "CompareSequenceByLogExpectation",
    "Description": "Multiple Sequence Comparison by Log Expectation is a tool used to compare protein or nucleic acid sequences. The average accuracy and speed are better than ClustalW2 or T-Coffee.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1417,
    "Tool Name": "DoubleSequenceGlobalAlignment",
    "Description": "This tool compares two sequences in global alignment style.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1418,
    "Tool Name": "ProteinSolubilityPredictor",
    "Description": "Predicts the solubility of a given protein sequence using a fine-tuned BioT5 model\n    for the protein solubility prediction task.\n    \n    Parameters:\n    - protein_sequence: The amino acid sequence of the protein.\n    \n    Returns:\n    - A Markdown formatted string containing the protein sequence and its solubility prediction.",
    "category": "Model Services",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1419,
    "Tool Name": "InherentDisorderedRegionsPredictor",
    "Description": "This tool predict the inherent disordered regions of the protein based on sequence.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1420,
    "Tool Name": "DoubleSequenceLocalAlignment",
    "Description": "This tool compare two sequences in local alignment style.Local alignment: Unlike global alignment, local alignment does not require alignment between two complete sequences, but instead uses certain local region fragments in each sequence for alignment. The demand for it lies in the discovery that although some protein sequences exhibit significant differences in overall sequence, they can independently perform the same function in certain local regions, and the sequences are quite conservative. At this point, it is obvious that relying on global alignment cannot obtain these locally similar sequences. Secondly, in eukaryotic genes, intron fragments exhibit great variability, while exon regions are relatively conserved. At this point, global alignment shows its limitations and cannot identify these local similarity sequences. Its representative is the Smith Waterman local alignment algorithm.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1421,
    "Tool Name": "ProteinMotifAnalysis",
    "Description": "This tool analyse the motif of the protein based on sequence. \"Motif (motif) refers to a conserved region in a DNA or protein sequence, or a small sequence pattern shared by a group of sequences. In biology, it is a data-driven mathematical statistical model. Functional prediction can be made based on protein sequence characteristics (such as protein motifs). Proteins with the same motif or domain can be classified into a large group called super families. Protein domains: are structural entities that typically represent a part of a protein's independent folding and movement functions. Therefore, proteins are often constructed from different combinations of these structural domains. At the motif level, the main emphasis is on the concept of structure rather than function, while domain emphasizes functional units, so it is mostly referred to as domain by function. If a protein has a Ca+2 binding domain, it means that the main function of a certain domain of the protein is to bind Ca+2, and the domain must have a Ca+2 binding motif (E-F hand motif) that provides Ca+2 binding.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1422,
    "Tool Name": "SequenceSimilarityCalculator",
    "Description": "Sequence similarity calculation takes a set of aligned sequences (FASTA or GCG format) as input to calculate their similarity",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1423,
    "Tool Name": "ORFFind",
    "Description": "The ORF search tool can help you find open reading frames in DNA sequences, and the returned results include the start and end positions of the ORF as well as the translation results of the open reading frames.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1424,
    "Tool Name": "TranslateDNAtoAminoAcidSequence",
    "Description": "This tool translate DNA sequence to protein(amino acid) sequence.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1425,
    "Tool Name": "RepeatDNASequenceSearch",
    "Description": "This tool search repeat DNA sequence in DNA sequence.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1426,
    "Tool Name": "RepeatProteinSequenceSearch",
    "Description": "This tool search repeat protein sequence in DNA sequence.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1427,
    "Tool Name": "PalindromicSequencesFinder",
    "Description": "This tool searches for palindrome sequences in the sequence and enters the length range of nucleic acid sequences and palindrome sequences in the text box below. The 'U' in the nucleic acid sequence will be replaced with 'T'",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1428,
    "Tool Name": "CalculateAminoAcidbyDegenerateCodon",
    "Description": "This tool calculate amino acid by degenerate codon. Input a degenerate codon (3nt) and calculate the amino acid encoded by that codon",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1429,
    "Tool Name": "ProteinNuclearLocalizationSequencePrediction",
    "Description": "This tool predict nuclear localization sequence of protein based on sequence. Nuclear localization sequence or signal - a structural domain of a protein, usually a short amino acid sequence that can interact with the nuclear carrier to transport the protein into the nucleus. NLS has no special requirements for the proteins it connects to and is not cleaved after nuclear input.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1430,
    "Tool Name": "SmallMoleculeSimilarityCalculation",
    "Description": "This tool calculate the similarity of small molecules based on SMILES.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1431,
    "Tool Name": "DNAMolecularWeightCalculator",
    "Description": "This tool calculate the molecular weight of DNA based on sequence.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1432,
    "Tool Name": "CpGIslandPrediction",
    "Description": "The CpG island prediction tool can predict potential CpG islands using the Gardiner Garden and Frommer (1987) method The calculation method is to use a 200bp window, with each shift of 1 bp. The CpG island is defined as an Obs/Exp value greater than 0.6 and a GC content greater than 50%. The calculation method for the number of CpG dimers in a window is to multiply the base number of 'C' in the window by the base number of 'G' in the window and then divide by the window length. CpG islands are typically found in the 5 'region of vertebrate genes.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1433,
    "Tool Name": "PCRPrimerProperties",
    "Description": "This tool calculate the properties of PCR primer based on sequence.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1434,
    "Tool Name": "AminoAcidStatistics",
    "Description": "This tool count the number of amino acids in the protein. Protein statistics summary: Based on the input protein sequence, count the number of each amino acid and calculate the proportion, which can quickly compare different sequences.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1435,
    "Tool Name": "SummaryEnzymeCleavageSites",
    "Description": "The enzyme digestion site summary tool counts the number and location of commonly used restriction endonucleating recognition sites in DNA sequences.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1436,
    "Tool Name": "RandomDNAGeneration",
    "Description": "This tool generate random DNA sequence.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1437,
    "Tool Name": "GenerateMoleculeDescription",
    "Description": "Given a molecule SELFIES, generates its English description using a pre-trained T5 model.\n    \n    Parameters:\n    - selfies: A string representing the molecule in SELFIES format.\n    \n    Returns:\n    - A Markdown string containing the molecule SELFIES and its generated English description.",
    "category": "Model Services",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1438,
    "Tool Name": "TexToMoleculeSELFIES",
    "Description": "Given a molecule description in English, generates its SELFIES and SMILES representation using a pre-trained T5 model.\n    \n    Parameters:\n    - description: A string containing the molecule description in English.\n    \n    Returns:\n    - A Markdown string containing the molecule description, SELFIES, and SMILES representation.\n    \n    Raises:\n    - ValueError: If the description input is not a valid string.",
    "category": "Model Services",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1439,
    "Tool Name": "PredictDrugTargetInteraction",
    "Description": "Predicts whether a given molecule (SELFIES format) and a protein sequence can interact with each other\n    and returns the result in Markdown format with a brief explanation.\n    \n    Parameters:\n    - input_data: A string containing the SELFIES of the molecule and the amino acid sequence of the protein,\n                  separated by a dot. Extra quotes around the input are handled.\n    \n    Returns:\n    - A Markdown formatted string containing the prediction and its explanation.",
    "category": "Model Services",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1440,
    "Tool Name": "PredictHumanProteinInteraction",
    "Description": "Predicts whether two given protein sequences can interact with each other\n    using a fine-tuned BioT5 model for the protein-protein interaction task with human dataset.\n    \n    Parameters:\n    - protein_sequence_1: The amino acid sequence of the first protein.\n    - protein_sequence_2: The amino acid sequence of the second protein.\n    \n    Returns:\n    - A Markdown formatted string containing the two protein sequences and the interaction prediction.",
    "category": "Model Services",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1441,
    "Tool Name": "PredictYeastProteinInteraction",
    "Description": "Predicts whether two given yeast protein sequences can interact with each other\n    using a fine-tuned BioT5 model for the protein-protein interaction task with yeast dataset.\n    \n    Parameters:\n    - input_data: A string containing two amino acid sequences of yeast proteins separated by a dot.\n                  Extra quotes around the input are handled.\n    \n    Returns:\n    - A Markdown formatted string containing the two yeast protein sequences and the interaction prediction.",
    "category": "Model Services",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1442,
    "Tool Name": "PredictProteinSolubility",
    "Description": "Predicts the solubility of a given protein sequence using a fine-tuned BioT5 model\n    for the protein solubility prediction task.\n    \n    Parameters:\n    - protein_sequence: The amino acid sequence of the protein.\n    \n    Returns:\n    - A Markdown formatted string containing the protein sequence and its solubility prediction.",
    "category": "Model Services",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1443,
    "Tool Name": "PredictProteinBinaryLocalization",
    "Description": "Predicts the binary localization of a given protein sequence using a fine-tuned BioT5 model\n    for the protein binary localization prediction task.\n    \n    Parameters:\n    - protein_sequence: The amino acid sequence of the protein.\n    \n    Returns:\n    - A Markdown formatted string containing the protein sequence and its binary localization prediction.",
    "category": "Model Services",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1444,
    "Tool Name": "CalculateMolecularWeight",
    "Description": "Calculates the molecular weight of a structure in a PDB file using MDAnalysis.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1445,
    "Tool Name": "GetAminoAcidFrequency",
    "Description": "Calculates the frequency of each amino acid in a protein sequence.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1446,
    "Tool Name": "GetReverseComplement",
    "Description": "Generates the reverse complement of a DNA sequence.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1447,
    "Tool Name": "CalculateHydrophobicityAndPolarity",
    "Description": "Calculates the hydrophobicity and polarity of a protein sequence.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Bio"
  },
  {
    "IDX": 1448,
    "Tool Name": "CASToPrice",
    "Description": "Fetches average price for multiple chemical substances identified by their CAS numbers.\n    This tool is usually a must-have if someone asks directly about the price of mof materials. And it is often used as a third tool.\n    Note: If someone ask the price of MOF materials, you need to convert the MOF material to Smiles first, then convert Smiles into CAS, and finally query the price through CAS.",
    "category": "Databases",
    "Server Name": "SciToolAgent-Mat"
  },
  {
    "IDX": 1449,
    "Tool Name": "MOFToSMILES",
    "Description": "Convert multiple MOF materials into SMILES representations and return the results in Markdown table format.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Mat"
  },
  {
    "IDX": 1450,
    "Tool Name": "SMILESToCAS",
    "Description": "Query a SMILES and return their CAS numbers in string.",
    "category": "Databases",
    "Server Name": "SciToolAgent-Mat"
  },
  {
    "IDX": 1451,
    "Tool Name": "MofLattice",
    "Description": "Obtain lattice structure information from the provided MOF cif file name.\n    Note: Please directly pass the cif file name, such as 'HKUST-1', as Action Input. The function will construct the file path based on a predefined directory.\n    The function returns the information in a Markdown table format.",
    "category": "Databases",
    "Server Name": "SciToolAgent-Mat"
  },
  {
    "IDX": 1452,
    "Tool Name": "GetStructureInfo",
    "Description": "Reads a structure file and returns basic information about the structure.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Mat"
  },
  {
    "IDX": 1453,
    "Tool Name": "CalculateDensity",
    "Description": "Calculates the density of a structure from a file.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Mat"
  },
  {
    "IDX": 1454,
    "Tool Name": "GetElementComposition",
    "Description": "Returns the elemental composition of a structure from a file.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Mat"
  },
  {
    "IDX": 1455,
    "Tool Name": "CalculateSymmetry",
    "Description": "Calculates the symmetry of a structure from a file.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Mat"
  },
  {
    "IDX": 1456,
    "Tool Name": "NameToSMILES",
    "Description": "Query a molecule name and return its SMILES string in Markdown format.",
    "category": "Databases",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1457,
    "Tool Name": "NameToCas",
    "Description": "Query a molecule name and return its CAS number in Markdown format.",
    "category": "Databases",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1458,
    "Tool Name": "FuncGroups",
    "Description": "Identify and list the functional groups in a molecule given its SMILES string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1459,
    "Tool Name": "SMILESToWeight",
    "Description": "Calculate the molecular weight of a molecule given its SMILES string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1460,
    "Tool Name": "MolSimilarity",
    "Description": "Calculate the Tanimoto similarity between two molecules given their SMILES strings.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1461,
    "Tool Name": "SMILESToInChI",
    "Description": "Convert a SMILES string to an InChI string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1462,
    "Tool Name": "InChIKeyToSMILES",
    "Description": "Convert an InChIKey string to a SMILES string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1463,
    "Tool Name": "InChIKeyToInChI",
    "Description": "Convert an InChIKey string to InChI.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1464,
    "Tool Name": "InChIKeyToMOL",
    "Description": "Convert an InChI string to a MOL string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1465,
    "Tool Name": "IsValidInChIKey",
    "Description": "Check if an InChIKey string is valid.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1466,
    "Tool Name": "InChIToSMILES",
    "Description": "Convert an InChI string to a SMILES string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1467,
    "Tool Name": "InChIToInChIKey",
    "Description": "Convert an InChI string to an InChIKey string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1468,
    "Tool Name": "InChIToCSID",
    "Description": "Convert InChI to ChemSpider ID.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1469,
    "Tool Name": "SMILEStoSELFIES",
    "Description": "Translates a SMILES string into its corresponding SELFIES string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1470,
    "Tool Name": "SELFIEStoSMILES",
    "Description": "Translates a SELFIES string into its corresponding SMILES string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1471,
    "Tool Name": "RandomMoelcule",
    "Description": "Generates a random molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1472,
    "Tool Name": "Length_SELFIES",
    "Description": "Computes the length of a SELFIES string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1473,
    "Tool Name": "Split_SELFIES",
    "Description": "Splits a SELFIES string into its individual tokens.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1474,
    "Tool Name": "GetAtomPairFingerprintAsBitVect",
    "Description": "Generate the atom pair fingerprint of a molecule as a SparseBitVect. \n    This fingerprint represents the presence of atom pairs, not just their counts.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1475,
    "Tool Name": "AssignPattyTypes",
    "Description": "Assign Patty types to the atoms of a molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1476,
    "Tool Name": "TestMolecule",
    "Description": "Perform a series of tests on a molecule, including sanitization, removal of hydrogens,\n    and canonicalization check. This function helps in validating the molecule's structure \n    and consistency.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1477,
    "Tool Name": "ShowMol",
    "Description": "Generate a molecule image from its SMILES representation and embed it directly in Markdown.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1478,
    "Tool Name": "TypeAtomsInMolecule",
    "Description": "Assigns EState types to each atom in a molecule based on its SMILES representation.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1479,
    "Tool Name": "CalculateEstateIndices",
    "Description": "Calculate EState indices for each atom in a molecule based on its SMILES representation.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1480,
    "Tool Name": "CalculateEstateVsa",
    "Description": "Calculate EState VSA indices for a molecule based on its SMILES representation.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1481,
    "Tool Name": "GenerateEstateFingerprint",
    "Description": "Generate the EState fingerprint for a molecule based on its SMILES representation.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1482,
    "Tool Name": "CalculateShapeSimilarity",
    "Description": "Calculate shape similarity scores using USRCAT for a list of molecules defined by their SMILES.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1483,
    "Tool Name": "CalculatePmi",
    "Description": "Calculate the normalized principal moments of inertia (NPR1 and NPR2) for a molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1484,
    "Tool Name": "CalculateDistanceMatrix",
    "Description": "Calculate the distance matrix for a list of molecules based on their fingerprints.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1485,
    "Tool Name": "ClusterMolecules",
    "Description": "Clusters molecules based on their fingerprints and returns the clustering results in Markdown format.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1486,
    "Tool Name": "ProcessFingerprintMol",
    "Description": "Process the molecular fingerprint generated by FingerprintMol function.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1487,
    "Tool Name": "FingerprintsFromSmiles",
    "Description": "Generate fingerprints for a list of SMILES strings.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1488,
    "Tool Name": "GetRdkFingerprintFromSmiles",
    "Description": "Generate an RDKit fingerprint from a SMILES string using default parameters.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1489,
    "Tool Name": "GenerateFraggleFragments",
    "Description": "Generate all possible Fraggle fragmentations for a molecule represented by a SMILES string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1490,
    "Tool Name": "CheckValidRingCut",
    "Description": "Check if the molecule represented by a SMILES string is a valid ring cut.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1491,
    "Tool Name": "BuildAtomPairFpFromSmiles",
    "Description": "Generate an Atom Pair Fingerprint from a SMILES string and display the results in a readable format.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1492,
    "Tool Name": "BuildTorsionsFpFromSmiles",
    "Description": "Generate a Torsions Fingerprint from a SMILES string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1493,
    "Tool Name": "BuildRdkitFpFromSmiles",
    "Description": "Generate an RDKit fingerprint from a SMILES string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1494,
    "Tool Name": "BuildPharm2DFpFromSmiles",
    "Description": "Generate a Pharm2D fingerprint from a SMILES string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1495,
    "Tool Name": "BuildMorganFpFromSmiles",
    "Description": "Generate a Morgan fingerprint from a SMILES string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1496,
    "Tool Name": "BuildAvalonFpFromSmiles",
    "Description": "Generate an Avalon fingerprint from a SMILES string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1497,
    "Tool Name": "ConvertSmilesToInchi",
    "Description": "Converts a SMILES string to its corresponding InChI string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1498,
    "Tool Name": "GenerateMolKeyFromSmiles",
    "Description": "Generates a molecular key for a given molecule represented by a SMILES string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1499,
    "Tool Name": "GetStereoCodeFromSmiles",
    "Description": "Generates the stereo code for a given molecule represented by a SMILES string.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1500,
    "Tool Name": "DetermineBondOrders",
    "Description": "The tool is used to determine the bond orders between atoms in a molecule based on their atomic coordinates. It assigns the connectivity information to the molecule by disregarding pre-existing bonds. This function is useful for inferring the chemical bonds in a molecule when the bond information is not already available or needs to be updated based on the 3D structure of the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1501,
    "Tool Name": "DetermineBonds",
    "Description": "The tool is used to determine the bond orders between atoms in a molecule based on their atomic coordinates. It assigns the connectivity information to the molecule by disregarding pre-existing bonds. This function is useful for inferring the chemical bonds in a molecule when the bond information is not already available or needs to be updated based on the 3D structure of the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1502,
    "Tool Name": "GetPatternFingerprint",
    "Description": "This tool is used to generate a pattern fingerprint for a molecule. The pattern fingerprint is a bit vector that encodes the presence or absence of particular substructures in the molecule. The substructures are defined by SMARTS patterns. The SMARTS patterns are converted to molecular fingerprints and then combined to generate the pattern fingerprint.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1503,
    "Tool Name": "IsSubstructure",
    "Description": "This tool is used to check if a molecule(target) is a substructure of another molecule(template). It returns true if the molecule is a substructure of the other molecule and false otherwise. The substructure search is performed by matching the SMARTS pattern of the query molecule to the target molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1504,
    "Tool Name": "GetTemplateMolecule",
    "Description": "This tool is used to get the template molecule from a TautomerQuery object.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1505,
    "Tool Name": "GetTautomers",
    "Description": "This tool obtains all possible tautomers of a TautomerQuery object. Tautomers are molecules that have the same atomic composition but differ in the connectivity of atoms. Retrieving all possible tautomers can help in understanding and analyzing changes in chemical reactions and molecular conformations.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1506,
    "Tool Name": "GetModifiedAtoms",
    "Description": "This tool is used to get the modified atoms of a TautomerQuery object. Modified atoms are atoms that have changed their connectivity in the tautomerization process.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1507,
    "Tool Name": "GetModifiedBonds",
    "Description": "This tool is used to get the modified bonds of a TautomerQuery object. Modified bonds are bonds that have changed their connectivity in the tautomerization process.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1508,
    "Tool Name": "GetSubstructMatches",
    "Description": "This tool is to search for substructures in a given target molecule that match the tautomer query.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1509,
    "Tool Name": "CanSerialize",
    "Description": "This tool is used to check if a TautomerQuery object can be serialized.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1510,
    "Tool Name": "AssignCIPLabels",
    "Description": "This tool is used to assign CIP labels to the atoms in a molecule. The CIP labels are used to describe the stereochemistry of the molecule. The labels are assigned based on the 3D structure of the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1511,
    "Tool Name": "Enumerate",
    "Description": "The rdkit.Chem.rdMolEnumerator.Enumerate function is used to perform enumeration on a given molecule and returns a MolBundle object containing multiple molecules generated during the enumeration process.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1512,
    "Tool Name": "Deprotect",
    "Description": "The rdkit.Chem.rdDeprotect.Deprotect function removes protecting groups from a molecule, returning the deprotected version.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1513,
    "Tool Name": "CondenseAbbreviationSubstanceGroups",
    "Description": "This tool finds and replaces abbreviation substance groups in a molecule, resulting in a compressed version of the molecule where the abbreviations are expanded.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1514,
    "Tool Name": "SlnToSmiles",
    "Description": "This tool is used to convert a SLN string to a SMILES string. Input SMILES directly without any other characters.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1515,
    "Tool Name": "CreateShingling",
    "Description": "This tool is used to create a shingling for a molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1516,
    "Tool Name": "EncodeMolecule",
    "Description": "This tool creates an MHFP vector from a molecule  using MHFP encoder, capturing structural information of the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1517,
    "Tool Name": "EncodeSECFP",
    "Description": "This tool creates an SECFP vector from a molecule using SECFP encoder, capturing structural information of the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1518,
    "Tool Name": "GetBCUT",
    "Description": "This tool computes the 2D BCUT descriptors for a given molecule, representing mass, Gasteiger charge, Crippen logP, and Crippen MR values.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1519,
    "Tool Name": "GetAutocorrelation2D",
    "Description": "This tool computes the 2D autocorrelation descriptors for a given molecule, capturing the spatial arrangement of atoms in the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1520,
    "Tool Name": "GetAutocorrelation3D",
    "Description": "This tool computes the 3D autocorrelation descriptors for a given molecule, capturing the spatial arrangement of atoms in the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1521,
    "Tool Name": "GetAsphericity",
    "Description": "This tool calculates the asphericity descriptor for a molecule, which measures how much the molecule deviates from a perfectly spherical shape.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1522,
    "Tool Name": "GetChi0n",
    "Description": "This tool calculates the chi^0 (chi-zero) cluster index, which represents a topological descriptor related to molecular branching.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1523,
    "Tool Name": "GetChi0v",
    "Description": "This function calculates the Chi^0v (Chi-zero-v) valence molecular graph index for a molecule, which is used to describe the topology of the molecule. It returns a float value.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1524,
    "Tool Name": "GetChi1n",
    "Description": "This tool calculates the chi^1 (chi-one) cluster index, which represents a topological descriptor related to molecular branching.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1525,
    "Tool Name": "GetChi1v",
    "Description": "This function calculates the Chi^1v (Chi-one-v) valence molecular graph index for a molecule, which is used to describe the topology of the molecule. It returns a float value.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1526,
    "Tool Name": "GetChi2n",
    "Description": "This tool calculates the chi^2 (chi-two) cluster index, which represents a topological descriptor related to molecular branching.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1527,
    "Tool Name": "GetChi2v",
    "Description": "This function calculates the Chi^2v (Chi-two-v) valence molecular graph index for a molecule, which is used to describe the topology of the molecule. It returns a float value.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1528,
    "Tool Name": "GetChi3n",
    "Description": "This tool calculates the chi^3 (chi-three) cluster index, which represents a topological descriptor related to molecular branching.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1529,
    "Tool Name": "GetChi3v",
    "Description": "This function calculates the Chi^3v (Chi-three-v) valence molecular graph index for a molecule, which is used to describe the topology of the molecule. It returns a float value.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1530,
    "Tool Name": "GetChi4n",
    "Description": "This tool calculates the chi^4 (chi-four) cluster index, which represents a topological descriptor related to molecular branching.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1531,
    "Tool Name": "GetChi4v",
    "Description": "This function calculates the Chi^4v (Chi-four-v) valence molecular graph index for a molecule, which is used to describe the topology of the molecule. It returns a float value.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1532,
    "Tool Name": "GetCoulombMat",
    "Description": "This tool calculates the Coulomb matrix for a molecule, which represents the electrostatic interactions between atoms in the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1533,
    "Tool Name": "GetCrippenDescriptors",
    "Description": "This function calculates the Wildman-Crippen logP and MR (molecular refractivity) values for a given molecule in RDKit.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1534,
    "Tool Name": "GetEEMCharges",
    "Description": "This function computes the EEM (Electronegativity Equalization Method) atomic partial charges for a given molecule using its atomic properties.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1535,
    "Tool Name": "GetEccentricity",
    "Description": "This function calculates the eccentricity of a molecule, which is a measure of its shape.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1536,
    "Tool Name": "GetExactMolceularWeight",
    "Description": "This function calculates the exact molecular weight of a molecule, which is the sum of the atomic weights of all atoms in the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1537,
    "Tool Name": "GetFractionCSP3",
    "Description": "This function calculates the fraction of sp3-hybridized carbon atoms in a molecule, which is a measure of its shape.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1538,
    "Tool Name": "GetGETAWAY",
    "Description": "This function calculates the GETAWAY descriptors for a molecule, which capture the shape and size of the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1539,
    "Tool Name": "GetHallKierAlpha",
    "Description": "This function calculates the Hall-Kier alpha index for a molecule, which is a measure of its shape.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1540,
    "Tool Name": "GetInertialShapeFactor",
    "Description": "This function calculates the Inertial Shape Factor of a molecule, which is a measure of its shape. The Inertial Shape Factor ranges from 0 to 1, where values closer to 1 indicate a more spherical shape and values closer to 0 indicate a more linear shape.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1541,
    "Tool Name": "GetKappa1",
    "Description": "This function computes the Kappa1 (¦Ê1) value of a molecule, which is a topological descriptor representing its shape complexity or branching degree. The Kappa1 value is a floating-point number calculated based on the molecular graph topology.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1542,
    "Tool Name": "GetKappa2",
    "Description": "This function computes the Kappa2 (¦Ê2) value of a molecule, which is a topological descriptor representing its shape complexity or branching degree. The Kappa2 value is a floating-point number calculated based on the molecular graph topology.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1543,
    "Tool Name": "GetKappa3",
    "Description": "This function computes the Kappa3 (¦Ê3) value of a molecule, which is a topological descriptor representing its shape complexity or branching degree. The Kappa3 value is a floating-point number calculated based on the molecular graph topology.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1544,
    "Tool Name": "GetLabuteASA",
    "Description": "This function calculates the Labute accessible surface area (ASA) value for a molecule, which is a measure of the solvent-accessible surface area of the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1545,
    "Tool Name": "GetMolFormula",
    "Description": "This function calculates the molecular formula of a molecule, which is a string representing the number and type of atoms in the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1546,
    "Tool Name": "GetMORSE",
    "Description": "This tool calculates the Molecule Representation of Structures based on Electron diffraction (MORSE) descriptors for a given molecule. MORSE descriptors provide a representation of molecular structures based on electron diffraction concepts.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1547,
    "Tool Name": "GetNPR1",
    "Description": "This function calculates the NPR1 (Normalized Principal Moments Ratio) descriptor for a molecule, which serves as a descriptor for the distribution of charges within the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1548,
    "Tool Name": "GetNPR2",
    "Description": "This function calculates the NPR2 (Normalized Principal Moments Ratio) descriptor for a molecule, which serves as a descriptor for the distribution of charges within the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1549,
    "Tool Name": "GetAliphaticCarbocyclesNum",
    "Description": "This function calculates the number of aliphatic carbocycles in a molecule. Aliphatic carbocycles are cyclic structures that contain at least one non-aromatic bond.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1550,
    "Tool Name": "GetAliphaticHeterocyclesNum",
    "Description": "This function calculates the number of aliphatic heterocycles in a molecule. Aliphatic heterocycles are cyclic structures that contain at least one non-aromatic bond and at least one heteroatom (an atom other than carbon).",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1551,
    "Tool Name": "GetAliphaticRingsNum",
    "Description": "This tool calculates the number of aliphatic rings in a molecule. Aliphatic rings are ring structures that contain at least one non-aromatic bond.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1552,
    "Tool Name": "GetAmideBondsNum",
    "Description": "This function calculates the number of amide bonds in a molecule. Amide bonds are chemical bonds formed between a carbonyl group and an amino group.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1553,
    "Tool Name": "GetAromaticCarbocyclesNum",
    "Description": "This function calculates the number of aromatic carbocycles in a molecule. Aromatic carbocycles are cyclic structures composed entirely of carbon atoms with alternating single and double bonds (aromaticity) in at least one ring.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1554,
    "Tool Name": "GetAromaticHeterocyclesNum",
    "Description": "This function calculates the number of aromatic heterocycles in a molecule. Aromatic heterocycles are cyclic structures that contain at least one heteroatom (an atom other than carbon) and exhibit aromaticity. es.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1555,
    "Tool Name": "GetAromaticRingsNum",
    "Description": "This tool calculates the number of aromatic rings in a molecule. Aromatic rings are cyclic structures composed of alternating single and double bonds (aromaticity) and exhibit stability due to delocalization of electrons.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1556,
    "Tool Name": "GetAtomStereoCentersNum",
    "Description": "This function calculates the number of atom stereo centers in a molecule. Atom stereo centers are atoms that are chiral centers and are not part of a ring structure.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1557,
    "Tool Name": "GetAtomsNum",
    "Description": "This function calculates the number of atoms in a molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1558,
    "Tool Name": "GetBridgeheadAtomsNum",
    "Description": "This function calculates the number of bridgehead atoms in a molecule. Bridgehead atoms are atoms that are part of a bridged ring structure.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1559,
    "Tool Name": "GetHBANum",
    "Description": "This function calculates the number of hydrogen bond acceptors (HBA) in a molecule. Hydrogen bond acceptors are atoms capable of forming hydrogen bonds by accepting a hydrogen atom.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1560,
    "Tool Name": "GetHBDNum",
    "Description": "This function calculates the number of hydrogen bond donors (HBD) in a molecule. Hydrogen bond donors are atoms capable of forming hydrogen bonds by donating a hydrogen atom.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1561,
    "Tool Name": "GetHeavyAtomsNum",
    "Description": "This tool calculates the number of heavy atoms in a molecule. Heavy atoms are atoms other than hydrogen.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1562,
    "Tool Name": "GetHeteroatomsNum",
    "Description": "This tool calculates the number of heteroatoms in a molecule. Heteroatoms are atoms other than carbon and hydrogen.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1563,
    "Tool Name": "GetHeterocyclesNum",
    "Description": "This tool calculates the number of heterocycles in a molecule. Heterocycles are cyclic structures that contain at least one heteroatom (an atom other than carbon).",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1564,
    "Tool Name": "GetLipinskiHBANum",
    "Description": "This tool calculates the number of Lipinski hydrogen bond acceptors (HBA) in a molecule, which is a measure used in drug-likeness evaluation according to Lipinski's rule of five. Lipinski's rule suggests that molecules with no more than five hydrogen bond acceptors tend to have better oral bioavailability.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1565,
    "Tool Name": "GetLipinskiHBDNum",
    "Description": "This tool calculates the number of Lipinski hydrogen bond donors (HBD) in a molecule, which is a measure used in drug-likeness evaluation according to Lipinski's rule of five.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1566,
    "Tool Name": "GetRingsNum",
    "Description": "This tool calculates the number of rings in a molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1567,
    "Tool Name": "GetRotatableBondsNum",
    "Description": "This tool calculates the number of rotatable bonds in a molecule. Rotatable bonds are single bonds that are not part of a ring structure and are not terminal (i.e., not connected to a hydrogen atom).",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1568,
    "Tool Name": "GetSaturatedCarbocyclesNum",
    "Description": "This function calculates the number of saturated carbocycles in a molecule. Saturated carbocycles are cyclic structures composed entirely of carbon atoms with single bonds (no double bonds).",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1569,
    "Tool Name": "GetSaturatedHeterocyclesNum",
    "Description": "This function calculates the number of saturated heterocycles in a molecule. Saturated heterocycles are cyclic structures that contain at least one heteroatom (an atom other than carbon) and are composed entirely of single bonds (no double bonds).",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1570,
    "Tool Name": "GetSaturatedRingsNum",
    "Description": "This tool calculates the number of saturated rings in a molecule. Saturated rings are ring structures composed entirely of single bonds (no double bonds).",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1571,
    "Tool Name": "GetSpiroAtomsNum",
    "Description": "This function calculates the number of spiro atoms in a molecule. Spiro atoms are atoms that are part of a spiro ring structure, which consists of two rings that share a single atom.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1572,
    "Tool Name": "GetUnspecifiedAtomStereoCentersNum",
    "Description": "This tool calculates the number of unspecified atomic stereocenters in a molecule. Unspecified atomic stereocenters are atoms that have the potential to be stereocenters but lack explicit specification of their stereochemistry.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1573,
    "Tool Name": "GenerateRDKFingerprintsFromCSV",
    "Description": "Generate RDKfingerprints for the SMILES strings in a CSV file and save to a new CSV file.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1574,
    "Tool Name": "GenerateMorganfingerprintsFromCSV",
    "Description": "Generate morgan fingerprints for the SMILES strings in a CSV file and save to a new CSV file.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1575,
    "Tool Name": "GenerateElectricalDescriptorsFromCSV",
    "Description": "Generate electrical RDKit descriptors for the SMILES strings in a CSV file and save to a new CSV file.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1576,
    "Tool Name": "MLPClassifier",
    "Description": "General MLP classifier function that predicts based on processed feature files.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1577,
    "Tool Name": "AdaBoostClassifier",
    "Description": "General AdaBoost classifier function that predicts based on processed feature files.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1578,
    "Tool Name": "RandomForestClassifier",
    "Description": "General Random Forest classifier function that predicts based on processed feature files.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1579,
    "Tool Name": "AssignOxidationNumbers",
    "Description": "Adds the oxidation number/state to the atoms of a molecule as property OxidationNumber on each atom.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1580,
    "Tool Name": "CalculatePBF",
    "Description": "This tool calculates the PBF (plane of best fit) descriptor for a given molecule. PBF is a molecular descriptor that characterizes the flatness or planarity of a molecule. It is calculated based on the arrangement of atoms in 3D space.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1581,
    "Tool Name": "CalculatePMI1",
    "Description": "This tool calculates the first principal moment of inertia (PMI1) for a given molecule. PMI1 is a molecular descriptor used to characterize the shape and spatial distribution of atoms in a molecule. PMI1 measures the asymmetry or elongation of a molecule along its principal axis. It provides information about the molecule's overall shape and can be useful in various computational chemistry and drug design applications.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1582,
    "Tool Name": "CalculatePMI2",
    "Description": "This tool is designed to compute the PMI2 (Partial Molecular Information 2) value of a molecule, which serves as a descriptor indicating the shape and structure of the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1583,
    "Tool Name": "CalculatePMI3",
    "Description": "This tool is designed to compute the PMI3 (Partial Molecular Information 3) value of a molecule, which serves as a descriptor characterizing the shape and structure of the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1584,
    "Tool Name": "CalculatePhi",
    "Description": "This tool calculates the Phi (¦Õ) angle of a molecule, which is a torsional angle describing the rotation about a single bond.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1585,
    "Tool Name": "CalculateRDF",
    "Description": "This tool calculates the RDF (Radial Distribution Function) descriptor for a given molecule. RDF is a molecular descriptor that characterizes the distribution of atoms in 3D space. It is calculated based on the distances between pairs of atoms in a molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1586,
    "Tool Name": "CalculateRadiusOfGyration",
    "Description": "This tool is designed to compute the radius of gyration for a given molecule, providing insights into its overall shape and compactness.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1587,
    "Tool Name": "CalculateSpherocityIndex",
    "Description": "This function calculates the sphericity index for a given molecule. Sphericity index is a measure of how close the shape of a molecule resembles a perfect sphere.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1588,
    "Tool Name": "CalculateTPSA",
    "Description": "This tool calculates the TPSA (Topological Polar Surface Area) descriptor for a given molecule, which is a measure of the accessible polar surface area in a molecule. TPSA is a molecular descriptor that characterizes the polarity and hydrophilicity of a molecule. It is calculated based on the distribution of polar atoms and bonds in a molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1589,
    "Tool Name": "CalculateWHIM",
    "Description": "This tool calculates the WHIM (Weighted Holistic Invariant Molecular) descriptor for a given molecule. WHIM is a molecular descriptor that characterizes the 3D shape and electronic properties of a molecule. It is calculated based on the distribution of atomic properties and their spatial arrangement in a molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1590,
    "Tool Name": "CustomPropertyVSA",
    "Description": "This function computes a custom property for a given molecule using the Van der Waals Surface Area (VSA) method, based on user-defined parameters.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1591,
    "Tool Name": "GetAtomFeature",
    "Description": "This function computes a set of atom features for a given molecule, including atomic number, valence, and hybridization.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1592,
    "Tool Name": "GetAtomPairFingerprint",
    "Description": "This function computes the atom pair for a given molecule. The atom pair fingerprint is a molecular descriptor that characterizes the presence of pairs of atoms in a molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1593,
    "Tool Name": "GetConnectivityInvariants",
    "Description": "This tool computes connectivity invariants, similar to ECFP (Extended Connectivity Fingerprints), for a given molecule. These invariants serve as a fingerprint representation of the molecule's structural connectivity, aiding in tasks such as similarity comparison and molecular structure representation.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1594,
    "Tool Name": "GetFeatureInvariants",
    "Description": "This tool computes feature invariants, similar to FCFP (Feature Centroid Fingerprints), for a given molecule. These invariants provide a fingerprint representation of the molecule's features, aiding in tasks such as similarity comparison and molecular structure analysis.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1595,
    "Tool Name": "GetAtomPairCode",
    "Description": "This function computes atom pair code (hash) for each atom in a molecular. The atom pair code is a molecular descriptor that characterizes the presence of pairs of atoms in a molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1596,
    "Tool Name": "GetHybridization",
    "Description": "This function computes the hybridization of each atom in a molecule. Hybridization is a property of an atom that characterizes its electron configuration and bonding behavior.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1597,
    "Tool Name": "GetRingSystems",
    "Description": "This function computes the ring systems of a molecule. A ring system is a set of rings that are connected to each other through shared atoms or bonds.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1598,
    "Tool Name": "GetMACCSKeysFingerprint",
    "Description": "This function computes the Molecular ACCess System keys fingerprint for a given molecule. The Molecular ACCess System keys fingerprint is a molecular descriptor that characterizes the presence of specific structural features in a molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1599,
    "Tool Name": "GetMorganFingerprint",
    "Description": "This tool computes the Morgan fingerprint for a given molecule. The Morgan fingerprint is a widely used method to encode molecular structure information. It captures the local chemical environments around each atom up to a specified radius.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1600,
    "Tool Name": "GetTopologicalTorsionFingerprint",
    "Description": "This tool computes the topological torsion fingerprint for a given molecule. The topological torsion fingerprint is a molecular descriptor that characterizes the presence of specific structural features in a molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1601,
    "Tool Name": "GetUSR",
    "Description": "The tool computes the USR (Ultrafast Shape Recognition) descriptor for a given conformer of a molecule and returns it as a list.The USR descriptor is a numerical representation of the shape of a molecule. It captures the 3D shape of a molecule in a compact form, making it particularly useful for comparing molecular shapes efficiently.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1602,
    "Tool Name": "GetUSRCAT",
    "Description": "This function is designed to compute the USRCAT (Ultrafast Shape Recognition with Coordinate Asymmetric Torsions) descriptor for a specified conformer of a molecule. The USRCAT descriptor is a compact representation of the molecular shape, which is useful for various cheminformatics applications such as similarity searching, clustering, and virtual screening.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1603,
    "Tool Name": "AddHydrogens",
    "Description": "This function is used to add hydrogen atoms to the molecular graph of a molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1604,
    "Tool Name": "AddWavyBondsForStereoAny",
    "Description": "This tool adds wavy bonds around double bonds with STEREOANY stereochemistry.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1605,
    "Tool Name": "AssignAtomChiralTagsFromStructure",
    "Description": "This tool sets chiral tags for atoms of the molecular based on the molParity property. This ensures proper definition of the molecule's stereochemistry for further analysis or visualization.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1606,
    "Tool Name": "AssignRadicals",
    "Description": "This tool is used to assign radical counts to atoms within a molecule. It takes a molecule SMILES as input and modifies it, assigning appropriate numbers of radicals to each atom within the molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1607,
    "Tool Name": "AssignStereoChemistry",
    "Description": "This tool is used for assigning Cahn¨CIngold¨CPrelog (CIP) stereochemistry to atoms (R/S) and double bonds (Z/E) within a molecule. Chiral atoms will have a property _CIPCode indicating their chiral code.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1608,
    "Tool Name": "GetAdjacencyMatrix",
    "Description": "This tool is used to obtain the adjacency matrix of a molecule. The adjacency matrix is a mathematical representation of a molecule where rows and columns correspond to atoms, and matrix elements represent whether pairs of atoms are adjacent (connected by a bond) or not.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1609,
    "Tool Name": "GetAllowNontetrahedralChirality",
    "Description": "This tool is used to determine whether recognition of non-tetrahedral chirality from 3D structures is enabled or not.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1610,
    "Tool Name": "GetDistanceMatrix",
    "Description": "The tool computes the topological distance matrix for a given molecule. This matrix provides information about the shortest path between pairs of atoms in the molecular graph, essentially indicating how many bonds need to be traversed to move from one atom to another.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1611,
    "Tool Name": "GetFormalCharge",
    "Description": "This tool is utilized to determine the total formal charge of a given molecule. Formal charge is a concept in chemistry that describes the net charge of an atom or a molecule, considering the redistribution of electrons based on electronegativity differences.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1612,
    "Tool Name": "GetFormalChargeOfAtoms",
    "Description": "This tool is utilized to determine the formal charge of each atom in a given molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1613,
    "Tool Name": "GetMolFrags",
    "Description": "This tool identifies disconnected fragments within a molecule and returns them as atom identifiers or molecules. It allows for flexible representation and manipulation of the fragments in further analysis.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1614,
    "Tool Name": "GetUseLegacyStereoPerception",
    "Description": "This tool is used to determine whether the legacy stereo perception code is being used. The legacy stereo perception code is an older implementation of stereochemistry perception in RDKit, which may be used for compatibility with older versions of the software.\n        smiles: a SMILES string. Please Input SMILES directly without any other characters¡£\n    Returns:\n        str: A markdown string indicating whether legacy stereo perception code is being used, or an error message.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1615,
    "Tool Name": "HapticBondsToDative",
    "Description": "This tool is used to convert a molecule that represents haptic bonds using a dummy atom with a dative bond to a metal atom into a molecule with explicit dative bonds from the atoms of the haptic group to the metal atom.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1616,
    "Tool Name": "HasQueryHs",
    "Description": "This tool is used to check if a molecule contains query H (hydrogen) atoms. Query hydrogens are special types of hydrogen atoms that are used to represent specific chemical environments or constraints in a molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1617,
    "Tool Name": "Kekulize",
    "Description": "This tool is used to perform Kekulization on a molecule. Kekulization is the process of converting aromatic bonds in a molecule to alternating single and double bonds, following the Kekul¨¦ structure representation.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1618,
    "Tool Name": "MergeQueryHs",
    "Description": "This tool is used to merge hydrogen atoms into their neighboring atoms as query atoms. This function is typically used to modify molecules by replacing explicit hydrogen atoms with query atoms, allowing for more flexible substructure searching or atom mapping.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1619,
    "Tool Name": "MurckoDecompose",
    "Description": "This tool is used to perform a Murcko decomposition on a molecule and return the scaffold. The Murcko scaffold represents the core structure of a molecule by removing side chains and retaining the ring system.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1620,
    "Tool Name": "RemoveHydrogens",
    "Description": "This tool is used to remove hydrogen atoms from a molecule's graph. This function is typically used to simplify molecular representations for further analysis or visualization.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1621,
    "Tool Name": "RemoveStereochemistry",
    "Description": "This tool is used to remove all stereochemistry information from a molecule. Stereochemistry information in a molecule refers to the spatial arrangement of atoms or groups around a stereocenter or a double bond.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1622,
    "Tool Name": "SetAromaticity",
    "Description": "This tool is used to perform aromaticity perception on a molecule, which means determining the aromaticity of atoms and bonds in the molecule. Aromaticity is a chemical property that describes the stability and reactivity of certain ring structures in organic molecules.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1623,
    "Tool Name": "SetBondStereoFromDirections",
    "Description": "This tool is used to set the cis/trans stereochemistry on double bonds based on the directions of neighboring bonds.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1624,
    "Tool Name": "IsSubstructof",
    "Description": "This tool is used to check if a molecule(target) is a substructure of another molecule(template). It returns true if the molecule is a substructure of the other molecule and false otherwise. The substructure search is performed by matching the SMARTS pattern of the query molecule to the target molecule.",
    "category": "Computational Tools",
    "Server Name": "SciToolAgent-Chem"
  },
  {
    "IDX": 1625,
    "Tool Name": "generate_presigned_url",
    "Description": "Generate pre-signed URLs for uploading and downloading Alibaba Cloud OSS objects",
    "category": "Computational Tools",
    "Server Name": "SCP-Workflow"
  }
];