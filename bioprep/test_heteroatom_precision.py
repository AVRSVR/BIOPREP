
import os
import sys
from pathlib import Path
import logging

# Add project root to path
sys.path.append(str(Path(__file__).parent.parent))

from bioprep.core.analyzer import analyze_structure
from bioprep.core.cleaner import clean_structure
from bioprep.core.io import load_pdb, save_pdb
from bioprep.core.protonator import add_hydrogens

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def create_mock_pdb(filepath, include_residues):
    """Creates a mock PDB file with specific residues for testing."""
    with open(filepath, 'w') as f:
        # Standard protein (ALA)
        f.write("ATOM      1  N   ALA A   1      29.132  18.044  16.123  1.00  0.00           N  \n")
        f.write("ATOM      2  CA  ALA A   1      29.132  19.044  16.123  1.00  0.00           C  \n")
        
        # Heteroatom 1: SO4 (Commonly removed)
        if "SO4" in include_residues:
            f.write("HETATM    3  S   SO4 A   2      35.123  20.123  15.000  1.00  0.00           S  \n")
            f.write("HETATM    4  O1  SO4 A   2      35.123  21.123  15.000  1.00  0.00           O  \n")
            f.write("HETATM    5  O2  SO4 A   2      36.123  20.123  15.000  1.00  0.00           O  \n")

        # Heteroatom 2: OLC (Long chain lipid, should be kept)
        if "OLC" in include_residues:
            f.write("HETATM    6  C1  OLC A   3      40.000  10.000  10.000  1.00  0.00           C  \n")
            f.write("HETATM    7  C2  OLC A   3      41.000  10.000  10.000  1.00  0.00           C  \n")

        f.write("END\n")

def test_surgical_removal():
    test_dir = Path("tmp_test_hets")
    test_dir.mkdir(exist_ok=True)
    
    input_pdb = test_dir / "input_test.pdb"
    output_pdb = test_dir / "output_test_cleaned.pdb"
    
    # 1. Create mock data
    logger.info("Step 1: Creating mock PDB with ALA, SO4, and OLC...")
    create_mock_pdb(input_pdb, ["SO4", "OLC"])
    
    # 2. Analyze
    logger.info("Step 2: Testing Analyzer...")
    struct = load_pdb(str(input_pdb))
    summary = analyze_structure(struct)
    detected_hets = summary['heteroatoms']
    logger.info(f"Detected Heteroatoms: {detected_hets}")
    
    assert "SO4" in detected_hets, "FAIL: SO4 not detected"
    assert "OLC" in detected_hets, "FAIL: OLC not detected"
    
    # 3. Clean (Surgical Removal: Remove SO4, Keep OLC)
    logger.info("Step 3: Testing Surgical Cleaner (Remove SO4, Keep OLC)...")
    select_obj = clean_structure(
        struct,
        target_chains=['A'],
        remove_water=True,
        remove_heteroatoms=['SO4'] # Only remove SO4
    )
    save_pdb(struct, str(output_pdb), select=select_obj)
    
    # 4. Verify results
    logger.info("Step 4: Verifying results in output PDB...")
    with open(output_pdb, 'r') as f:
        content = f.read()
        
    has_so4 = "SO4" in content
    has_olc = "OLC" in content
    
    if not has_so4 and has_olc:
        logger.info("✅ SUCCESS: SO4 was removed, OLC was preserved surgically.")
    else:
        logger.error(f"❌ FAILURE: SO4 present: {has_so4}, OLC present: {has_olc}")
        sys.exit(1)

    # 5. Pipeline Check (Protonation)
    logger.info("Step 5: Testing Protonation Pipeline integrity...")
    try:
        final_pdb = test_dir / "final_test_prot.pdb"
        add_hydrogens(str(output_pdb), str(final_pdb), ph=7.0)
        with open(final_pdb, 'r') as f:
            final_content = f.read()
            has_olc_final = "OLC" in final_content
            if has_olc_final:
                logger.info("✅ SUCCESS: OLC survived the protonation pipeline.")
            else:
                logger.error("❌ FAILURE: OLC was stripped during protonation!")
                sys.exit(1)
    except Exception as e:
        logger.error(f"❌ Pipeline Crash: {e}")
        sys.exit(1)

    logger.info("\n✨ ALL HETEROATOM PRECISION TESTS PASSED! ✨")

if __name__ == "__main__":
    test_surgical_removal()
