import numpy as np
import os
from Bio.PDB import PDBParser, NeighborSearch, Selection
from typing import List, Dict, Any


class BindingSiteAnalyzer:
    def __init__(self, pdb_path: str):
        self.pdb_path = pdb_path
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure("protein", pdb_path)
        self.atoms = list(self.structure.get_atoms())
        self.coords = np.array([atom.get_coord() for atom in self.atoms])

        # Residue properties mapping
        self.RESIDUE_PROPS = {
            'ALA': 'HYDROPHOBIC', 'VAL': 'HYDROPHOBIC', 'LEU': 'HYDROPHOBIC', 'ILE': 'HYDROPHOBIC',
            'PRO': 'HYDROPHOBIC', 'PHE': 'HYDROPHOBIC', 'TRP': 'HYDROPHOBIC', 'MET': 'HYDROPHOBIC',
            'SER': 'POLAR', 'THR': 'POLAR', 'CYS': 'POLAR', 'TYR': 'POLAR', 'ASN': 'POLAR', 'GLN': 'POLAR',
            'ASP': 'CHARGED_NEG', 'GLU': 'CHARGED_NEG',
            'LYS': 'CHARGED_POS', 'ARG': 'CHARGED_POS', 'HIS': 'CHARGED_POS',
            'GLY': 'NEUTRAL'
        }

    def analyze(self) -> List[Dict[str, Any]]:
        """Main analysis pipeline."""
        pockets = self._detect_pockets()

        results = []
        for i, pocket in enumerate(pockets):
            pocket_points = np.array(pocket['points'])
            site_info = self._analyze_specific_site(pocket_points)
            site_info['id'] = i + 1
            results.append(site_info)

        # Sort by drugability score descending
        results.sort(key=lambda x: x['drugability_score'], reverse=True)
        
        # Filter out tiny pockets and cap at top 5
        results = [r for r in results if r['volume'] >= 50]
        results = results[:5]
        
        # Normalize IDs to be sequential
        for idx, res in enumerate(results):
            res['id'] = idx + 1
            
        return results

    def get_summary(self, results: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Generates a high-level summary of all discovered sites."""
        if not results:
            return {
                "text": "No significant binding pockets were detected.",
                "site_count": 0,
                "primary_volume": 0.0,
                "total_volume": 0.0,
                "total_features": 0
            }

        # Primary site is the first one in sorted results
        primary_site = results[0]
        primary_vol = primary_site['volume']
        primary_score = primary_site['drugability_score']
        
        total_cavity_vol = sum(r['volume'] for r in results)
        total_features = sum(len(r['pharmacophore_points']) for r in results)

        summary_text = (
            f"Analysis complete. Found **{len(results)} distinct pockets** "
            f"with a **Total Cavity Space of {total_cavity_vol:.1f} Å³**. "
            f"The **Primary Binding Site** has a volume of **{primary_vol:.1f} Å³** "
            f"and a drugability score of **{primary_score}**, suggesting a "
            "strong host-guest fit for typical small molecules."
        )

        return {
            "text": summary_text,
            "site_count": len(results),
            "primary_volume": primary_vol,
            "total_volume": total_cavity_vol,
            "total_features": total_features
        }


    def _detect_pockets(self, grid_res: float = 1.5, probe_radius: float = 2.8) -> List[np.ndarray]:
        """
        Grid-based pocket detection using geometric enclosure filtering.
        """
        from scipy.spatial import KDTree
        from sklearn.cluster import DBSCAN
        import numpy as np

        min_coords = np.min(self.coords, axis=0) - 5
        max_coords = np.max(self.coords, axis=0) + 5

        # DYNAMIC GRID SCALING: Protect against MemoryError on huge complexes (e.g. virus capsids)
        box_vol = np.prod(max_coords - min_coords)
        target_max_points = 250000
        if (box_vol / (grid_res**3)) > target_max_points:
            grid_res = (box_vol / target_max_points) ** (1/3.0)
            grid_res = max(1.5, round(grid_res, 1))

        x = np.arange(min_coords[0], max_coords[0], grid_res)
        y = np.arange(min_coords[1], max_coords[1], grid_res)
        z = np.arange(min_coords[2], max_coords[2], grid_res)
        grid_points = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)

        tree = KDTree(self.coords)
        
        # Pre-filter: grid points must be near at least one protein atom (between probe and 7.5A)
        # But for enclosure check, we need all neighbors within a larger radius.
        min_dists, _ = tree.query(grid_points, k=1)
        potential_mask = (min_dists > probe_radius) & (min_dists < 7.5)
        potential_points = grid_points[potential_mask]
        
        if len(potential_points) == 0:
            return []

        # Find neighbors for all potential points in one go if possible, or per-point
        # Dist-based search for enclosure
        all_neighbor_idxs = tree.query_ball_point(potential_points, 13.0)

        # ── ENCLOSURE CHECK (LIGSITE style - Optimized) ────────────────────────
        enclosed_points = []
        max_dist = 12.0
        
        dimensions = [
            (0, 1.0), (0, -1.0), # X+, X-
            (1, 1.0), (1, -1.0), # Y+, Y-
            (2, 1.0), (2, -1.0)  # Z+, Z-
        ]

        for i, point in enumerate(potential_points):
            nb_idx = all_neighbor_idxs[i]
            if len(nb_idx) < 10: continue
            
            nearby_coords = self.coords[nb_idx]
            # Use distance squared to avoid unnecessary sqrt operations
            vectors = nearby_coords - point
            dist_sq = np.sum(vectors**2, axis=1)
            
            hit_dirs = 0
            for axis_idx, sign in dimensions:
                # Fast projection for axis-aligned rays
                projs = vectors[:, axis_idx] * sign
                
                # Pythagorean: perp_dist^2 = dist^2 - proj^2
                # We need perp_dist < 2.5 => perp_dist^2 < 6.25
                perp_sq = dist_sq - projs**2
                
                mask = (projs > 1.5) & (projs < max_dist) & (perp_sq < 6.25)
                if np.any(mask):
                    hit_dirs += 1
                    # Performance: 3 directions is the threshold for a pocket point.
                    # Early exit once we cross the threshold for this grid point.
                    if hit_dirs >= 3:
                        enclosed_points.append(point)
                        break

        if not enclosed_points:
            return []

        # ── CLUSTERING ───────────────────────────────────────────────────────
        points_np = np.array(enclosed_points)
        clustering = DBSCAN(eps=grid_res * 1.5, min_samples=10).fit(points_np)
        labels = clustering.labels_
        
        pockets = []
        unique_labels = set(labels)
        for label in unique_labels:
            if label == -1: continue
            
            cluster_points = points_np[labels == label]
            volume = len(cluster_points) * (grid_res ** 3)
            
            # Standard drug-like pocket is 400-1200 A^3
            drugability = 0.0
            if volume > 300:
                drugability = min(0.95, 0.4 + (volume / 2000))
            else:
                drugability = volume / 750
            
            centroid = np.mean(cluster_points, axis=0)
            pockets.append({
                'id': int(label) + 1,
                'centroid': centroid.tolist(),
                'volume': round(volume, 1),
                'points': cluster_points.tolist(),
                'drugability_score': round(drugability, 2)
            })

        return sorted(pockets, key=lambda x: x['volume'], reverse=True)

    def _analyze_specific_site(self, points: np.ndarray) -> Dict[str, Any]:
        """Detailed analysis of a single pocket cluster, with improved drugability scoring."""
        centroid = np.mean(points, axis=0)
        
        # FIX: More accurate volume: count grid cells with the actual grid res (1.5Å)
        # Each grid voxel = 1.5³ = 3.375 Å³
        voxel_volume = 1.5 ** 3
        volume = len(points) * voxel_volume

        # Find nearby residues (10Å from centroid)
        ns = NeighborSearch(self.atoms)
        nearby_atoms = ns.search(centroid, 10.0)
        nearby_residues = list(set([a.get_parent() for a in nearby_atoms]))

        res_list = []
        props_counts = {'HYDROPHOBIC': 0, 'POLAR': 0, 'CHARGED': 0}

        for res in nearby_residues:
            res_name = res.get_resname()
            prop = self.RESIDUE_PROPS.get(res_name, 'OTHER')
            res_list.append(f"{res_name}{res.get_id()[1]}")

            if 'HYDROPHOBIC' in prop:
                props_counts['HYDROPHOBIC'] += 1
            elif 'POLAR' in prop:
                props_counts['POLAR'] += 1
            elif 'CHARGED' in prop:  # Catches CHARGED_POS and CHARGED_NEG both
                props_counts['CHARGED'] += 1

        # Pharmacophore features
        pharmacophores = self._predict_pharmacophores(nearby_atoms, centroid)

        # ── IMPROVED Drugability Score ────────────────────────────────────────
        # Factor 1: Volume (target: 300–1000 Å³ for drug-like molecules)
        volume_score = 0.0
        if volume < 100:
            volume_score = 0.1  # Too small
        elif volume < 300:
            volume_score = 0.4 + (volume - 100) / 500
        elif volume <= 1000:
            volume_score = 1.0
        else:
            volume_score = max(0.3, 1.0 - (volume - 1000) / 3000)  # Too large = less specific

        # Factor 2: Property diversity (having all 3 property types is best)
        prop_diversity = len([v for v in props_counts.values() if v > 0]) / 3.0

        # Factor 3: Hydrophobic balance (~40% hydrophobic is ideal)
        total_props = sum(props_counts.values()) + 1
        hphob_ratio = props_counts['HYDROPHOBIC'] / total_props
        balance_score = max(0.0, 1.0 - abs(hphob_ratio - 0.4) * 2)

        # Factor 4: Concavity (enclosure) — FIX: NEW — pocket should be geometrically enclosed
        # Estimate using standard deviation of pocket point distances from centroid.
        # A tightly clustered enclosed pocket has lower std dev relative to centroid distance.
        dists_from_centroid = np.linalg.norm(points - centroid, axis=1)
        mean_dist = np.mean(dists_from_centroid)
        std_dist = np.std(dists_from_centroid)
        # Concavity proxy: enclosed = mean_dist high, std low (tight sphere-like pocket)
        if mean_dist > 0:
            concavity_score = min(1.0, max(0.0, 1.0 - (std_dist / mean_dist)))
        else:
            concavity_score = 0.0

        # Factor 5: Pharmacophore density (more pharmacophore features in pocket = better)
        pharm_score = min(1.0, len(pharmacophores) / 10.0)

        # Weighted final score
        drugability = (
            volume_score    * 0.30 +
            prop_diversity  * 0.20 +
            balance_score   * 0.20 +
            concavity_score * 0.20 +
            pharm_score     * 0.10
        )
        drugability = min(0.99, max(0.05, drugability))

        return {
            'centroid': centroid.tolist(),
            'volume': round(volume, 1),
            'residues': res_list[:15],
            'properties': props_counts,
            'pharmacophore_points': pharmacophores,
            'drugability_score': round(drugability, 2),
            'concavity': round(concavity_score, 2)
        }

    def _predict_pharmacophores(self, nearby_atoms, centroid) -> List[Dict[str, Any]]:
        """
        Predict pharmacophore features from atoms lining the pocket.
        
        Improvements:
        - Aromatic centroid now computed as average of ring carbons (not just CG)
        - Backbone acceptors/donors also included
        """
        points = []
        AROMATIC_RES = {'PHE', 'TRP', 'TYR', 'HIS'}
        AROMATIC_RING_ATOMS = {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'CH2', 'NE1'}

        # Group aromatic residues first to compute ring centroids
        aromatic_ring_coords = {}  # {residue_id: [coords, ...]}

        for atom in nearby_atoms:
            dist = np.linalg.norm(atom.get_coord() - centroid)
            if dist > 8.5:
                continue

            at_name = atom.get_name().strip()
            res = atom.get_parent()
            res_name = res.get_resname().strip()
            res_label = f"{res_name}{res.id[1]}"

            # H-bond Acceptors: Sidechain + backbone oxygens
            if at_name in ['OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'O']:
                points.append({
                    'type': 'ACCEPTOR',
                    'coords': atom.get_coord().tolist(),
                    'label': f"Acc-{res_label}"
                })

            # H-bond Donors: Sidechain nitrogens + backbone NH
            elif at_name in ['NZ', 'NH1', 'NH2', 'ND1', 'NE2', 'ND2', 'NE1', 'N']:
                points.append({
                    'type': 'DONOR',
                    'coords': atom.get_coord().tolist(),
                    'label': f"Don-{res_label}"
                })

            # Hydrophobic: aliphatic residues
            elif res_name in ['ALA', 'VAL', 'LEU', 'ILE', 'MET'] and at_name in ['CB', 'CG', 'CD', 'CE']:
                points.append({
                    'type': 'HYDROPHOBIC',
                    'coords': atom.get_coord().tolist(),
                    'label': f"Hphob-{res_label}"
                })

            # FIX: Aromatic — collect ring atom coords for centroid calculation
            elif res_name in AROMATIC_RES and at_name in AROMATIC_RING_ATOMS:
                rid = res.get_full_id()
                if rid not in aromatic_ring_coords:
                    aromatic_ring_coords[rid] = {'coords': [], 'label': f"Arom-{res_label}"}
                aromatic_ring_coords[rid]['coords'].append(atom.get_coord())

        # FIX: Now add one HYDROPHOBIC point per aromatic residue at true ring centroid
        for rid, data in aromatic_ring_coords.items():
            if len(data['coords']) >= 3:  # Need at least 3 ring atoms for a meaningful centroid
                ring_centroid = np.mean(data['coords'], axis=0)
                points.append({
                    'type': 'HYDROPHOBIC',
                    'coords': ring_centroid.tolist(),
                    'label': data['label']
                })

        return points[:40]
