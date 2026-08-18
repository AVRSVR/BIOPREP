import urllib.request
import sys
sys.path.insert(0, 'bioprep')
from bioprep.core import io, cleaner, protonator, minimizer, site_analyzer

print('Downloading 1STP...')
try:
    urllib.request.urlretrieve('https://files.rcsb.org/download/1STP.pdb', 'bioprep/tmp/1stp.pdb')
except Exception as e:
    print('Failed to download 1STP:', e)
    sys.exit(1)

print('\n========== TEST 1: KEEP LIGAND (BTN) ==========')
structure1 = io.load_pdb('bioprep/tmp/1stp.pdb')
select_obj1 = cleaner.clean_structure(structure1, remove_water=True, remove_heteroatoms=[], keep_structural_waters=False)
io.save_pdb(structure1, 'bioprep/tmp/1stp_keep.pdb', select=select_obj1)
protonator.add_hydrogens('bioprep/tmp/1stp_keep.pdb', 'bioprep/tmp/1stp_keep_prot.pdb', add_missing_atoms=True)

try:
    stats1 = minimizer.minimize_structure('bioprep/tmp/1stp_keep_prot.pdb', 'bioprep/tmp/1stp_keep_min.pdb')
    print('Minimizer Result (Keep BTN):', stats1.get('error', 'Success!'))
except Exception as e:
    print('Minimizer Error:', e)

print('\n========== TEST 2: REMOVE LIGAND (BTN) ==========')
structure2 = io.load_pdb('bioprep/tmp/1stp.pdb')
select_obj2 = cleaner.clean_structure(structure2, remove_water=True, remove_heteroatoms=['ALL'], keep_structural_waters=False)
io.save_pdb(structure2, 'bioprep/tmp/1stp_remove.pdb', select=select_obj2)
protonator.add_hydrogens('bioprep/tmp/1stp_remove.pdb', 'bioprep/tmp/1stp_remove_prot.pdb', add_missing_atoms=True)

try:
    stats2 = minimizer.minimize_structure('bioprep/tmp/1stp_remove_prot.pdb', 'bioprep/tmp/1stp_remove_min.pdb')
    print('Minimizer Result (Remove BTN):')
    print('  Before:', stats2['energy_before_kJ_mol'])
    print('  After:', stats2['energy_after_kJ_mol'])
    
    print('\nSite Analysis on Empty Pocket:')
    analyzer = site_analyzer.BindingSiteAnalyzer('bioprep/tmp/1stp_remove_min.pdb')
    sites = analyzer.analyze()
    print(f'Found {len(sites)} binding pockets!')
    if sites:
        print(f'Top pocket volume: {sites[0]["volume"]} A^3 (This is the empty Biotin pocket!)')

except Exception as e:
    print('Minimizer Error:', e)
