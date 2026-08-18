import argparse
import os
import tempfile
from bioprep.core.io import load_pdb, save_pdb
from bioprep.core.cleaner import clean_structure
from bioprep.core.protonator import add_hydrogens

def main():
    parser = argparse.ArgumentParser(description="BioPrep: Prepare protein structures for molecular docking.")
    parser.add_argument("--input", required=True, help="Path to input PDB file")
    parser.add_argument("--output", required=True, help="Path to save cleaned PDB file")
    parser.add_argument("--chain", help="Optional: specific chain to keep (e.g., A)", default=None)
    parser.add_argument("--ph", type=float, default=7.4, help="pH for adding hydrogens (default: 7.4)")

    args = parser.parse_args()

    print(f"[*] Loading {args.input}...")
    try:
        structure = load_pdb(args.input)
    except Exception as e:
        print(f"[!] Error loading PDB file: {e}")
        return

    print("[*] Cleaning structure (removing waters, heteroatoms)...")
    if args.chain:
        print(f"[*] Filtering to keep only Chain {args.chain}")
    select_obj = clean_structure(structure, target_chain=args.chain)

    # Save to a temporary file
    fd, temp_pdb_path = tempfile.mkstemp(suffix=".pdb")
    os.close(fd)
    
    try:
        save_pdb(structure, temp_pdb_path, select=select_obj)
        
        print(f"[*] Adding hydrogens at pH {args.ph}...")
        add_hydrogens(temp_pdb_path, args.output, ph=args.ph)
        
        print(f"[*] Structure saved to {args.output}")
        print("[*] BioPrep finished successfully!")
    except Exception as e:
        print(f"[!] Error during processing: {e}")
    finally:
        # Clean up temporary file
        if os.path.exists(temp_pdb_path):
            os.remove(temp_pdb_path)

if __name__ == "__main__":
    main()
