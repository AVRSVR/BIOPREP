"""
Grid-based binding-site detection and characterisation.

Pockets are found by marking grid points that sit in the shell just outside
the protein surface, keeping those enclosed by protein in at least three of
six axis directions (a LIGSITE-style test), and clustering what survives.
"""

import numpy as np
from Bio.PDB import PDBParser, NeighborSearch
from typing import List, Dict, Any

from .residues import is_water

# Residue property classes used for pocket characterisation.
RESIDUE_PROPS = {
    'ALA': 'HYDROPHOBIC', 'VAL': 'HYDROPHOBIC', 'LEU': 'HYDROPHOBIC', 'ILE': 'HYDROPHOBIC',
    'PRO': 'HYDROPHOBIC', 'PHE': 'HYDROPHOBIC', 'TRP': 'HYDROPHOBIC', 'MET': 'HYDROPHOBIC',
    'SER': 'POLAR', 'THR': 'POLAR', 'CYS': 'POLAR', 'TYR': 'POLAR', 'ASN': 'POLAR', 'GLN': 'POLAR',
    'ASP': 'CHARGED_NEG', 'GLU': 'CHARGED_NEG',
    'LYS': 'CHARGED_POS', 'ARG': 'CHARGED_POS', 'HIS': 'CHARGED_POS',
    'GLY': 'NEUTRAL',
}

SIDECHAIN_ACCEPTORS = {'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH'}
SIDECHAIN_DONORS = {'NZ', 'NH1', 'NH2', 'ND1', 'NE2', 'ND2', 'NE1'}
BACKBONE_ACCEPTOR = 'O'
BACKBONE_DONOR = 'N'
HYDROPHOBIC_RES = {'ALA', 'VAL', 'LEU', 'ILE', 'MET'}
HYDROPHOBIC_ATOMS = {'CB', 'CG', 'CD', 'CE'}
AROMATIC_RES = {'PHE', 'TRP', 'TYR', 'HIS'}
AROMATIC_RING_ATOMS = {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CH2', 'NE1'}

MAX_PHARMACOPHORES = 40
MIN_POCKET_VOLUME = 50.0
MAX_POCKETS = 5


def _enclosure_directions():
    """
    The directions the enclosure test scans, as unit vectors.

    LIGSITE scans seven axes - the three Cartesian axes plus the four cubic
    body diagonals - in both senses, giving fourteen rays. Scanning only the
    six axis-aligned rays makes the result depend on how the molecule happens
    to be oriented in the file: the same cavity, rotated, is sampled by rays
    that meet its walls at different angles. Measured on streptavidin, a 45
    degree rotation changed total detected cavity volume by nearly half.
    """
    directions = []
    for axis in np.eye(3):
        directions.extend([axis, -axis])
    for signs in ((1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)):
        diagonal = np.array(signs, dtype=float)
        diagonal /= np.linalg.norm(diagonal)
        directions.extend([diagonal, -diagonal])
    return np.array(directions)


ENCLOSURE_DIRECTIONS = _enclosure_directions()      # 14 rays over 7 axes

#: A point counts as enclosed when this fraction of the rays is blocked.
#: Half matches the previous 3-of-6 behaviour on the axis-aligned subset.
ENCLOSURE_FRACTION = 0.5
MIN_BLOCKED_DIRECTIONS = int(round(len(ENCLOSURE_DIRECTIONS) * ENCLOSURE_FRACTION))

# Feature ranking used when trimming to MAX_PHARMACOPHORES. Sidechain and
# aromatic features characterise a pocket; backbone N/O occur in every residue
# and would otherwise crowd everything else out of the list.
PRIORITY = {'AROMATIC': 0, 'SIDECHAIN': 1, 'HYDROPHOBIC': 2, 'BACKBONE': 3}


class BindingSiteAnalyzer:
    def __init__(self, pdb_path: str, ignore_hydrogens: bool = True,
                 ignore_waters: bool = True):
        """
        Load a structure for pocket analysis.

        Hydrogens and waters are excluded by default. Pocket geometry is
        defined by heavy atoms; including hydrogens shrinks every measured
        cavity, and retained waters fill the very pockets being looked for —
        which matters because this runs on protonated, prepared structures.
        """
        self.pdb_path = pdb_path
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure("protein", pdb_path)

        model = next(iter(self.structure), None)
        source = model if model is not None else self.structure

        atoms = []
        for atom in source.get_atoms():
            if ignore_hydrogens and atom.element == 'H':
                continue
            if ignore_hydrogens and atom.get_name().strip().startswith('H'):
                continue
            if ignore_waters and is_water(atom.get_parent().get_resname()):
                continue
            atoms.append(atom)

        self.atoms = atoms
        self.coords = (np.array([a.get_coord() for a in atoms])
                       if atoms else np.empty((0, 3)))
        self.RESIDUE_PROPS = RESIDUE_PROPS

    # ── public API ───────────────────────────────────────────────────────────

    def analyze(self) -> List[Dict[str, Any]]:
        """Detect pockets and characterise each one."""
        if len(self.coords) < 10:
            return []

        pockets, grid_res = self._detect_pockets()

        results = []
        for pocket in pockets:
            info = self._analyze_specific_site(np.array(pocket['points']), grid_res)
            results.append(info)

        results = [r for r in results if r['volume'] >= MIN_POCKET_VOLUME]
        results.sort(key=lambda r: r['drugability_score'], reverse=True)
        results = results[:MAX_POCKETS]

        for index, result in enumerate(results):
            result['id'] = index + 1
        return results

    def get_summary(self, results: List[Dict[str, Any]]) -> Dict[str, Any]:
        """High-level summary across all discovered sites."""
        if not results:
            return {
                "text": "No significant binding pockets were detected.",
                "site_count": 0,
                "primary_volume": 0.0,
                "total_volume": 0.0,
                "total_features": 0,
            }

        primary = results[0]
        total_volume = sum(r['volume'] for r in results)
        total_features = sum(len(r['pharmacophore_points']) for r in results)

        text = (
            f"Analysis complete. Found **{len(results)} distinct pockets** "
            f"with a **Total Cavity Space of {total_volume:.1f} Å³**. "
            f"The **Primary Binding Site** has a volume of "
            f"**{primary['volume']:.1f} Å³** and a drugability score of "
            f"**{primary['drugability_score']}**."
        )

        return {
            "text": text,
            "site_count": len(results),
            "primary_volume": primary['volume'],
            "total_volume": total_volume,
            "total_features": total_features,
        }

    # ── pocket detection ─────────────────────────────────────────────────────

    def _detect_pockets(self, grid_res: float = 1.5, probe_radius: float = 2.8):
        """
        Find candidate pocket clusters.

        Returns ``(pockets, grid_res)``. The resolution is returned because it
        may have been coarsened for a large structure, and the volume of a
        pocket depends on it — assuming 1.5 Å regardless understated every
        volume for exactly the structures that needed scaling.
        """
        from scipy.spatial import KDTree
        from sklearn.cluster import DBSCAN

        min_coords = np.min(self.coords, axis=0) - 5
        max_coords = np.max(self.coords, axis=0) + 5

        # Coarsen the grid rather than exhaust memory on very large complexes.
        box_volume = float(np.prod(max_coords - min_coords))
        target_max_points = 250000
        if (box_volume / (grid_res ** 3)) > target_max_points:
            grid_res = max(1.5, round((box_volume / target_max_points) ** (1 / 3.0), 1))

        axes = [np.arange(min_coords[i], max_coords[i], grid_res) for i in range(3)]
        grid_points = np.array(np.meshgrid(*axes)).T.reshape(-1, 3)

        tree = KDTree(self.coords)

        # Candidates sit in the shell outside the surface: past the probe
        # radius, but not out in bulk solvent.
        min_dists, _ = tree.query(grid_points, k=1)
        candidates = grid_points[(min_dists > probe_radius) & (min_dists < 7.5)]
        if len(candidates) == 0:
            return [], grid_res

        neighbour_lists = tree.query_ball_point(candidates, 13.0)

        enclosed = []
        max_dist = 12.0
        directions = ENCLOSURE_DIRECTIONS
        needed = MIN_BLOCKED_DIRECTIONS
        remaining_after = len(directions) - np.arange(len(directions)) - 1

        for i, point in enumerate(candidates):
            neighbours = neighbour_lists[i]
            if len(neighbours) < 10:
                continue

            vectors = self.coords[neighbours] - point
            dist_sq = np.sum(vectors ** 2, axis=1)

            hits = 0
            for index, direction in enumerate(directions):
                # General projection onto a unit vector, so a diagonal ray is
                # measured the same way an axis-aligned one is.
                projections = vectors @ direction
                # Perpendicular distance from the ray, via Pythagoras.
                perpendicular_sq = dist_sq - projections ** 2
                blocked = ((projections > 1.5) & (projections < max_dist)
                           & (perpendicular_sq < 6.25))
                if np.any(blocked):
                    hits += 1
                    if hits >= needed:
                        enclosed.append(point)
                        break
                elif hits + remaining_after[index] < needed:
                    # Cannot reach the threshold with the rays that are left.
                    break

        if not enclosed:
            return [], grid_res

        points = np.array(enclosed)
        labels = DBSCAN(eps=grid_res * 1.5, min_samples=10).fit(points).labels_

        pockets = []
        for label in set(labels):
            if label == -1:
                continue
            cluster = points[labels == label]
            pockets.append({
                'centroid': np.mean(cluster, axis=0).tolist(),
                'points': cluster.tolist(),
            })

        pockets.sort(key=lambda p: len(p['points']), reverse=True)
        return pockets, grid_res

    # ── pocket characterisation ──────────────────────────────────────────────

    def _analyze_specific_site(self, points: np.ndarray, grid_res: float) -> Dict[str, Any]:
        """Characterise one pocket cluster."""
        centroid = np.mean(points, axis=0)
        volume = len(points) * (grid_res ** 3)

        neighbour_search = NeighborSearch(self.atoms)
        nearby_atoms = neighbour_search.search(centroid, 10.0)

        residues = {atom.get_parent() for atom in nearby_atoms}
        props_counts = {'HYDROPHOBIC': 0, 'POLAR': 0, 'CHARGED': 0}
        labelled = []

        for residue in residues:
            name = residue.get_resname().strip()
            prop = RESIDUE_PROPS.get(name)
            if prop is None:
                continue  # ligands, ions and anything non-standard
            labelled.append((
                float(np.linalg.norm(residue.child_list[0].get_coord() - centroid)),
                f"{name}{residue.get_id()[1]}",
            ))
            if prop == 'HYDROPHOBIC':
                props_counts['HYDROPHOBIC'] += 1
            elif prop == 'POLAR':
                props_counts['POLAR'] += 1
            elif prop.startswith('CHARGED'):
                props_counts['CHARGED'] += 1

        # Report the closest lining residues rather than an arbitrary slice.
        labelled.sort()
        residue_list = [name for _, name in labelled[:15]]

        pharmacophores = self._predict_pharmacophores(nearby_atoms, centroid)

        # ── Drugability, five weighted factors ───────────────────────────────
        if volume < 100:
            volume_score = 0.1
        elif volume < 300:
            volume_score = 0.4 + (volume - 100) / 500
        elif volume <= 1000:
            volume_score = 1.0
        else:
            volume_score = max(0.3, 1.0 - (volume - 1000) / 3000)

        prop_diversity = len([v for v in props_counts.values() if v > 0]) / 3.0

        total_props = sum(props_counts.values()) + 1
        hydrophobic_ratio = props_counts['HYDROPHOBIC'] / total_props
        balance_score = max(0.0, 1.0 - abs(hydrophobic_ratio - 0.4) * 2)

        distances = np.linalg.norm(points - centroid, axis=1)
        mean_dist = float(np.mean(distances))
        concavity_score = (
            min(1.0, max(0.0, 1.0 - (float(np.std(distances)) / mean_dist)))
            if mean_dist > 0 else 0.0
        )

        pharm_score = min(1.0, len(pharmacophores) / 10.0)

        drugability = (
            volume_score * 0.30
            + prop_diversity * 0.20
            + balance_score * 0.20
            + concavity_score * 0.20
            + pharm_score * 0.10
        )
        drugability = min(0.99, max(0.05, drugability))

        return {
            'centroid': centroid.tolist(),
            'volume': round(volume, 1),
            'grid_resolution': grid_res,
            'residues': residue_list,
            'properties': props_counts,
            'pharmacophore_points': pharmacophores,
            'drugability_score': round(drugability, 2),
            'concavity': round(concavity_score, 2),
            # The five weighted factors, exposed so a ranking can be argued
            # with rather than taken on faith. This is an unvalidated
            # heuristic: on HIV-1 protease (1HSG) the real inhibitor site is
            # found but ranks last of five, because it is larger than the
            # 300-1000 A^3 band the volume term rewards.
            'drugability_factors': {
                'volume': {'score': round(volume_score, 3), 'weight': 0.30},
                'property_diversity': {'score': round(prop_diversity, 3),
                                       'weight': 0.20},
                'hydrophobic_balance': {'score': round(balance_score, 3),
                                        'weight': 0.20},
                'concavity': {'score': round(concavity_score, 3), 'weight': 0.20},
                'pharmacophore_density': {'score': round(pharm_score, 3),
                                          'weight': 0.10},
            },
        }

    def _predict_pharmacophores(self, nearby_atoms, centroid) -> List[Dict[str, Any]]:
        """
        Predict pharmacophore features from the atoms lining a pocket.

        Features are ranked before trimming to ``MAX_PHARMACOPHORES``. Backbone
        N and O appear in every residue, so trimming in discovery order used to
        fill the list with backbone atoms and drop the aromatic ring centroids
        entirely — they were appended last.
        """
        features = []
        aromatic_rings = {}

        for atom in nearby_atoms:
            distance = float(np.linalg.norm(atom.get_coord() - centroid))
            if distance > 8.5:
                continue

            name = atom.get_name().strip()
            residue = atom.get_parent()
            res_name = residue.get_resname().strip()
            label = f"{res_name}{residue.id[1]}"

            if name in SIDECHAIN_ACCEPTORS:
                features.append((PRIORITY['SIDECHAIN'], distance, {
                    'type': 'ACCEPTOR', 'coords': atom.get_coord().tolist(),
                    'label': f"Acc-{label}"}))
            elif name in SIDECHAIN_DONORS:
                features.append((PRIORITY['SIDECHAIN'], distance, {
                    'type': 'DONOR', 'coords': atom.get_coord().tolist(),
                    'label': f"Don-{label}"}))
            elif name == BACKBONE_ACCEPTOR:
                features.append((PRIORITY['BACKBONE'], distance, {
                    'type': 'ACCEPTOR', 'coords': atom.get_coord().tolist(),
                    'label': f"Acc-{label}(bb)"}))
            elif name == BACKBONE_DONOR:
                features.append((PRIORITY['BACKBONE'], distance, {
                    'type': 'DONOR', 'coords': atom.get_coord().tolist(),
                    'label': f"Don-{label}(bb)"}))
            elif res_name in HYDROPHOBIC_RES and name in HYDROPHOBIC_ATOMS:
                features.append((PRIORITY['HYDROPHOBIC'], distance, {
                    'type': 'HYDROPHOBIC', 'coords': atom.get_coord().tolist(),
                    'label': f"Hphob-{label}"}))

            if res_name in AROMATIC_RES and name in AROMATIC_RING_ATOMS:
                ring = aromatic_rings.setdefault(
                    residue.get_full_id(), {'coords': [], 'label': f"Arom-{label}"})
                ring['coords'].append(atom.get_coord())

        # One feature per aromatic ring, at the true geometric centroid.
        for ring in aromatic_rings.values():
            if len(ring['coords']) >= 3:
                ring_centroid = np.mean(ring['coords'], axis=0)
                features.append((
                    PRIORITY['AROMATIC'],
                    float(np.linalg.norm(ring_centroid - centroid)),
                    {'type': 'AROMATIC', 'coords': ring_centroid.tolist(),
                     'label': ring['label']},
                ))

        features.sort(key=lambda item: (item[0], item[1]))
        return [feature for _, _, feature in features[:MAX_PHARMACOPHORES]]
