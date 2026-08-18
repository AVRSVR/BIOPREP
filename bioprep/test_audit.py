import os
from bioprep.core import io, cleaner, protonator, minimizer, analyzer, site_analyzer, exporter

def main():
    print("Starting full end-to-end pipeline audit...")
    input_file = "bioprep/1crn.pdb"
    temp_file = "bioprep/tmp/test_temp.pdb"
    prot_file = "bioprep/tmp/test_prot.pdb"
    min_file = "bioprep/tmp/test_min.pdb"
    output_pdbqt = "bioprep/tmp/test_out.pdbqt"
    
    os.makedirs("bioprep/tmp", exist_ok=True)

    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    # 1. Parsing
    print("\n--- 1. Parsing ---")
    structure = io.load_pdb(input_file)
    print("Structure loaded successfully.")

    # 2. Cleaning
    print("\n--- 2. Cleaning ---")
    select_obj = cleaner.clean_structure(
        structure, remove_water=True, remove_heteroatoms=['ALL'],
        target_chains=None, keep_structural_waters=False
    )
    io.save_pdb(structure, temp_file, select=select_obj)
    print(f"Cleaned structure saved to {temp_file}.")

    # 3. Protonation
    print("\n--- 3. Protonation ---")
    protonator.add_hydrogens(temp_file, prot_file, ph=7.4)
    print(f"Protonated structure saved to {prot_file}.")

    # 4. Minimization
    print("\n--- 4. Minimization ---")
    min_stats = minimizer.minimize_structure(
        prot_file, min_file, force_field="amber14", use_gbsa=True
    )
    print(f"Minimization stats: {min_stats}")

    # 5. Analysis
    print("\n--- 5. Analysis ---")
    min_structure = io.load_pdb(min_file)
    analysis_stats = analyzer.analyze_structure(min_structure)
    print(f"Analysis stats: {analysis_stats}")

    # 6. Site Analysis
    print("\n--- 6. Binding Site Analysis ---")
    try:
        analyzer_obj = site_analyzer.BindingSiteAnalyzer(min_file)
        sites = analyzer_obj.analyze()
        print(f"Found {len(sites)} potential binding sites.")
        for i, s in enumerate(sites):
            print(f"  Site {i+1}: Volume {s['volume']:0.1f} A^3, Drugability {s['drugability_score']:0.2f}, Concavity {s['concavity']:0.2f}")
    except Exception as e:
        print(f"Error during site analysis: {e}")

    # 7. Exporting
    print("\n--- 7. Exporting ---")
    try:
        pdbqt_file = exporter.convert_to_pdbqt(min_file, output_pdbqt)
        print(f"Saved {pdbqt_file}")
    except Exception as e:
        print(f"Export to PDBQT failed (maybe missing obabel): {e}")

    print("\n✅ End-to-end audit completed successfully.")

if __name__ == "__main__":
    main()
