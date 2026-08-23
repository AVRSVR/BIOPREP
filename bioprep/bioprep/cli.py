import argparse
import os
import shutil
import sys
import tempfile

from bioprep.core.io import load_pdb, save_pdb, ensure_pdb
from bioprep.core.cleaner import clean_structure
from bioprep.core.protonator import add_hydrogens
from bioprep.core.minimizer import minimize_structure, STATUS_FULL


def main():
    parser = argparse.ArgumentParser(
        description="BioPrep: Prepare protein structures for molecular docking."
    )
    parser.add_argument("--input", required=True, help="Path to input PDB file")
    parser.add_argument("--output", required=True, help="Path to save the prepared PDB file")
    parser.add_argument("--chain", action="append", metavar="ID",
                        help="Chain to keep; repeat the flag for several (e.g. --chain A --chain B)")
    parser.add_argument("--ph", type=float, default=7.4,
                        help="pH for adding hydrogens (default: 7.4)")
    parser.add_argument("--keep-water", action="store_true",
                        help="Keep all water molecules")
    parser.add_argument("--keep-structural-water", action="store_true",
                        help="Keep only waters within 4A of the protein")
    parser.add_argument("--remove-het", action="append", default=[], metavar="RESNAME",
                        help="Heteroatom to remove; use ALL to remove every heteroatom")
    parser.add_argument("--protect", action="append", default=[], metavar="RESNAME",
                        help="Ligand to keep even when --remove-het ALL is given")
    parser.add_argument("--minimize", action="store_true",
                        help="Run energy minimization after protonation")
    parser.add_argument("--force-field", default="amber14", choices=["amber14", "charmm36"])
    parser.add_argument("--no-gbsa", action="store_true",
                        help="Disable GBSA implicit solvent during minimization")

    args = parser.parse_args()

    # One working directory for every intermediate, removed on the way out.
    # ensure_pdb would otherwise make its own temporary directory for an mmCIF
    # conversion that nothing ever cleaned up.
    workdir = tempfile.mkdtemp(prefix="bioprep_cli_")

    try:
        return _run(args, workdir)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def _run(args, workdir):
    print(f"[*] Loading {args.input}...")
    try:
        args.input, converted_from = ensure_pdb(args.input, workdir)
        if converted_from:
            print(f"    Converted {converted_from.upper()} input to PDB")
        structure = load_pdb(args.input)
    except Exception as e:
        print(f"[!] Error loading PDB file: {e}")
        return 1

    if args.chain:
        print(f"[*] Keeping only chain(s): {', '.join(args.chain)}")
    print("[*] Cleaning structure...")
    select_obj = clean_structure(
        structure,
        target_chains=args.chain,
        remove_water=not args.keep_water,
        remove_heteroatoms=args.remove_het,
        keep_structural_waters=args.keep_structural_water,
        protect_ligands=args.protect,
    )

    cleaned_path = os.path.join(workdir, "cleaned.pdb")
    protonated_path = os.path.join(workdir, "protonated.pdb")

    try:
        save_pdb(structure, cleaned_path, select=select_obj,
                 source_pdb=args.input)

        print(f"[*] Adding hydrogens at pH {args.ph}...")
        prot = add_hydrogens(cleaned_path, protonated_path, ph=args.ph)
        if prot['ligands_preserved']:
            print(f"    Ligands held out of PDBFixer: {', '.join(prot['ligands_preserved'])}")
        for warning in prot['warnings']:
            print(f"[!] {warning}")

        if args.minimize:
            print(f"[*] Minimizing with {args.force_field}...")
            stats = minimize_structure(
                protonated_path, args.output,
                force_field=args.force_field, use_gbsa=not args.no_gbsa,
            )
            if stats['error']:
                print(f"[!] Minimization failed: {stats['error']}")
                print("[!] The structure was written WITHOUT minimization.")
            else:
                print(f"    {stats['energy_before_kJ_mol']} -> "
                      f"{stats['energy_after_kJ_mol']} kJ/mol "
                      f"(status: {stats['status']})")
                for warning in stats['warnings']:
                    print(f"[!] {warning}")
        else:
            shutil.copy2(protonated_path, args.output)

        print(f"[*] Structure saved to {args.output}")
        print("[*] BioPrep finished successfully!")
        return 0

    except Exception as e:
        print(f"[!] Error during processing: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
