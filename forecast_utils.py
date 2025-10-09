import os
import re
from typing import Optional
from evcouplings.compare import DistanceMap
import numpy as np
import pandas as pd

class PrepTools:
    @staticmethod
    def load_protein_sequences(path: str) -> pd.DataFrame:
        """Load protein sequences catalog CSV."""
        return pd.read_csv(path)

    @staticmethod
    def load_who_catalog(path: str, sheet="Mutation_catalogue") -> pd.DataFrame:
        """Load WHO catalog Excel and flatten multi-index columns."""
        df = pd.read_excel(path, sheet_name=sheet, header=[0, 1])
        df.columns = ["_".join(filter(None, col)).strip() for col in df.columns.values]
        return df

    @staticmethod
    def preprocess_catalog(df: pd.DataFrame, genes_of_interest: list) -> pd.DataFrame:
        """Filter, clean, and convert catalog variants to one-letter format."""
        catalog_df = df[[
            'drug_Unnamed: 0_level_1',
            'variant (common_name)_Unnamed: 2_level_1',
            'DATASET_Present_R',
            'DATASET_Present_S',
            'FINAL CONFIDENCE GRADING_Unnamed: 51_level_1'
        ]].copy()

        catalog_df.columns = ['drug', 'variant', 'Present_R', 'Present_S', 'confidence']
        catalog_df = catalog_df.dropna(subset=['drug', 'variant'])

        # Subset by target genes
        filtered = catalog_df[catalog_df['variant'].str.startswith(tuple(genes_of_interest))]

        # Remove indels
        filtered = filtered[~filtered['variant'].str.contains('ins|del')]

        # Extract one-letter mutation
        filtered['one_letter_mutation'] = filtered['variant'].str.extract(r'_(\w+\d+\w+)')
        filtered['gene'] = filtered['variant'].str.split('_').str[0]

        # Keep valid mutations only
        valid_pattern = r'^\w+_[A-Z]\d+[A-Z]$'
        filtered = filtered[filtered['variant'].str.match(valid_pattern)]
        filtered = filtered.dropna(subset=['one_letter_mutation']).drop_duplicates()

        # Reorder
        first_cols = ['drug', 'gene', 'variant']
        remaining_cols = [c for c in filtered.columns if c not in first_cols]
        filtered = filtered[first_cols + remaining_cols]

        return filtered

    @staticmethod
    def save_preprocessed(filtered: pd.DataFrame, out_path: str):
        """Save processed catalog to CSV."""
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        filtered.to_csv(out_path, index=False)
        return out_path


class MutationTools:
    @staticmethod
    def apply_mutation(wildtype_seq: str, mutation: str, log: Optional[list] = None):
        """
        Apply a point mutation (e.g. 'A36V').
        Returns tuple (orig, new, pos, mutated_seq) or None if invalid.
        """
        if len(mutation) < 3 or not mutation[1:-1].isdigit():
            if log is not None:
                log.append({"type": "bad_format", "mutation": mutation})
            return None

        orig, pos, new = mutation[0], int(mutation[1:-1]), mutation[-1]

        if not (1 <= pos <= len(wildtype_seq)):
            if log is not None:
                log.append({"type": "out_of_range", "mutation": mutation, "pos": pos, "seq_len": len(wildtype_seq)})
            return None

        if wildtype_seq[pos-1] != orig:
            if log is not None:
                log.append({"type": "ref_mismatch", "mutation": mutation,
                            "expected": orig, "found": wildtype_seq[pos-1], "pos": pos})
            return None

        mutated = wildtype_seq[:pos-1] + new + wildtype_seq[pos:]
        return orig, new, pos, mutated

    @staticmethod
    def generate_mutated_fastas(protein_df: pd.DataFrame, catalog_df: pd.DataFrame, output_dir: str) -> list:
        """
        Write FASTA files with wildtype + mutated sequences for each gene.
        Returns a list of successful mutations.
        """
        os.makedirs(output_dir, exist_ok=True)
        results, log = [], []
        for _, row in protein_df.iterrows():
            gene, rv_id, wt_seq = row['gene'], row['RV'], row['protein_sequence']
            out_file = os.path.join(output_dir, f"{rv_id}_{gene}.fasta")
            with open(out_file, "w") as f:
                f.write(f">{rv_id}|{gene}|Wildtype\n{wt_seq}\n")
                for _, mut_row in catalog_df[catalog_df['gene'] == gene].iterrows():
                    mutation = mut_row['one_letter_mutation']
                    result = MutationTools.apply_mutation(wt_seq, mutation, log=log)
                    if result:
                        orig, new, pos, mut_seq = result
                        label = f"{gene}_p.{orig}{pos}{new}"
                        f.write(f">{rv_id}|{gene}|{label}\n{mut_seq}\n")
                        results.append(result)
        return results, log

class CatalogNormalizer:
    three_to_one = {
        'Ala':'A','Cys':'C','Asp':'D','Glu':'E','Phe':'F','Gly':'G','His':'H',
        'Ile':'I','Lys':'K','Leu':'L','Met':'M','Asn':'N','Pro':'P','Gln':'Q',
        'Arg':'R','Ser':'S','Thr':'T','Val':'V','Trp':'W','Tyr':'Y'
    }

    # @staticmethod
    # def convert_three_to_one(mutation: str) -> Optional[str]:
    #     """Convert WHO-2023 style (three-letter) to one-letter (e.g., p.Ser450Leu -> p.S450L)."""
    #     if not isinstance(mutation, str):
    #         return None
    #     m = re.match(r"[a-zA-Z0-9]+_p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", mutation)
    #     if not m:
    #         return None
    #     orig3, pos, new3 = m.groups()
    #     try:
    #         return f"p.{CatalogNormalizer.three_to_one[orig3]}{pos}{CatalogNormalizer.three_to_one[new3]}"
    #     except KeyError:
    #         return None

    @staticmethod
    def convert_three_to_one(mutation: str) -> Optional[str]:
        """Convert WHO-2023 style (three-letter) to one-letter, e.g., rpoB_p.Ser450Leu → S450L."""
        if not isinstance(mutation, str):
            return None
        m = re.match(r"[a-zA-Z0-9]+_p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", mutation)
        if not m:
            return None
        orig3, pos, new3 = m.groups()
        try:
            return f"{CatalogNormalizer.three_to_one[orig3]}{pos}{CatalogNormalizer.three_to_one[new3]}"
        except KeyError:
            return None


    @staticmethod
    def load_and_standardize(path: str, year: int, genes_of_interest: list) -> pd.DataFrame:
        """Load WHO catalog (2021 or 2023), standardize to common format with one-letter mutations."""
        if year == 2021:
            # load and flatten
            df = pd.read_excel(path, sheet_name="Mutation_catalogue", header=[0, 1])
            df.columns = ["_".join(filter(None, col)).strip() for col in df.columns.values]
            df = df[[
                'drug_Unnamed: 0_level_1',
                'variant (common_name)_Unnamed: 2_level_1',
                'DATASET_Present_R',
                'DATASET_Present_S',
                'FINAL CONFIDENCE GRADING_Unnamed: 51_level_1'
            ]].copy()
            df.columns = ['drug', 'variant', 'Present_R', 'Present_S', 'confidence']
            # filter genes
            filtered = df.dropna(subset=['drug', 'variant'])
            filtered = filtered[filtered['variant'].str.startswith(tuple(genes_of_interest))]
            filtered = filtered[~filtered['variant'].str.contains('ins|del')]
            # extract one-letter mutation from pattern "gene_A123B"
            filtered['one_letter_mutation'] = filtered['variant'].str.extract(r'_(\w+\d+\w+)')
            filtered['gene'] = filtered['variant'].str.split('_').str[0]

        elif year == 2023:
            df = pd.read_excel(path, sheet_name="Catalogue_master_file", header=2)
            df = df[['drug','gene','mutation','variant','effect','Present_R','Present_S','FINAL CONFIDENCE GRADING']].copy()
            df = df.rename(columns={'FINAL CONFIDENCE GRADING':'confidence'})
            filtered = df.dropna(subset=['drug','variant'])
            filtered = filtered[filtered['variant'].str.startswith(tuple(genes_of_interest))]
            filtered = filtered[~filtered['variant'].str.contains('ins|del')]
            # normalize mutation string
            filtered['one_letter_mutation'] = filtered['variant'].apply(CatalogNormalizer.convert_three_to_one)

        else:
            raise ValueError("Year must be 2021 or 2023")

        # drop bad
        filtered = filtered.dropna(subset=['one_letter_mutation']).drop_duplicates()
        # reorder
        first_cols = ['drug', 'gene', 'variant']
        remaining_cols = [c for c in filtered.columns if c not in first_cols]
        filtered = filtered[first_cols + remaining_cols]
        return filtered


class FeatureTools:
    @staticmethod
    def add_frequency(df: pd.DataFrame) -> pd.DataFrame:
        """Add frequency = R / (R+S), drop raw counts."""
        df = df.copy()
        if "Present_R" in df and "Present_S" in df:
            df["frequency"] = df["Present_R"] / (df["Present_R"] + df["Present_S"])
            df = df.drop(columns=["Present_R","Present_S"])
        return df

    @staticmethod
    def compute_delta_z(fasta_dir: str, model_name="esm2_t6_8M_UR50D",
                        device=None, layer=6, out_path: Optional[str]=None) -> pd.DataFrame:
        """
        Compute delta-z for mutated FASTA files: Euclidean distance between
        wildtype and mutant embeddings.
        """
        import esm, torch
        from Bio import SeqIO

        device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
        model = model.eval().to(device)
        batch_converter = alphabet.get_batch_converter()

        # def get_embedding(seq: str):
        #     seq = seq.upper().replace("*","")
        #     _, _, toks = batch_converter([("seq", seq)])
        #     toks = toks.to(device)
        #     with torch.no_grad():
        #         res = model(toks, repr_layers=[layer])
        #     emb = res["representations"][layer][0].cpu().numpy()
        #     return emb.mean(axis=0)
        def get_embedding(seq: str):
            seq = seq.upper().replace("*", "")
            # keep only valid AA chars
            seq = "".join([aa for aa in seq if aa in alphabet.standard_toks])
            _, _, toks = batch_converter([("seq", seq)])
            toks = toks.to(device)
            with torch.no_grad():
                res = model(toks, repr_layers=[layer])
            emb = res["representations"][layer][0].cpu().numpy()
            return emb.mean(axis=0)


        results = []
        for fname in os.listdir(fasta_dir):
            if not fname.endswith(".fasta"):
                continue
            file_path = os.path.join(fasta_dir, fname)
            records = list(SeqIO.parse(file_path, "fasta"))
            rv_id, gene_name = fname.replace(".fasta","").split("_")
            wild = next((r for r in records if "Wildtype" in r.description), None)
            if wild is None: continue
            wt_emb = get_embedding(str(wild.seq))
            for r in records:
                if r.description == wild.description: continue
                mut_emb = get_embedding(str(r.seq))
                delta_z = np.linalg.norm(mut_emb - wt_emb)
                results.append({
                    "filename": fname,
                    "mutation": r.description,
                    "delta_z": delta_z,
                    "rv": rv_id,
                    "gene": gene_name
                })
        df = pd.DataFrame(results)
        if out_path:
            df.to_csv(out_path, index=False)
        return df

    @staticmethod
    def merge_delta_z(catalog: pd.DataFrame, delta_z_df: pd.DataFrame) -> pd.DataFrame:
        """Join delta-z with catalog on gene+mutation (normalized)."""
        df = delta_z_df.copy()
        # strip gene_ and p. from mutation field
        df["mutation"] = df["mutation"].apply(lambda x: x.split("_")[1] if "_" in x else x)
        df["mutation"] = df["mutation"].str.replace("p.","", regex=False)
        merged = pd.merge(df, catalog, left_on=["gene","mutation"],
                          right_on=["gene","one_letter_mutation"], how="left")
        merged = merged.drop(columns=["variant","rv","filename"], errors="ignore").drop_duplicates()
        return merged

    @staticmethod
    def add_llr_contributions(catalog: pd.DataFrame, protein_sequences: pd.DataFrame,
                              model_path: str, top_k: int = 10, device=None) -> pd.DataFrame:
        """
        Compute LLR scores + top-K dimension contributions for each variant.
        """
        import esm, torch
        device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        model, alphabet = esm.pretrained.load_model_and_alphabet_local(model_path)
        model = model.to(device).eval()
        bc = alphabet.get_batch_converter()
        mask_idx = alphabet.mask_idx
        get_id = alphabet.get_idx

        llr_contrib_data = []
        catalog = catalog.copy()
        if "Wildtype_AA" not in catalog.columns or "Mutated_AA" not in catalog.columns:
            # Parse mutations like A123B
            catalog["Wildtype_AA"] = catalog["one_letter_mutation"].str[0]
            catalog["position"] = catalog["one_letter_mutation"].str.extract(r"(\d+)").astype(int)
            catalog["Mutated_AA"] = catalog["one_letter_mutation"].str[-1]

        for gene, subset in catalog.groupby("gene"):
            prot_row = protein_sequences[protein_sequences["gene"]==gene]
            if prot_row.empty: continue
            seq = prot_row["protein_sequence"].iloc[0]
            _, _, toks = bc([("seq", seq)])
            toks = toks.to(device)
            for idx, row in subset.iterrows():
                pos = row["position"]
                wt, mt = row["Wildtype_AA"], row["Mutated_AA"]
                try:
                    wt_id, mt_id = get_id(wt), get_id(mt)
                except KeyError: continue
                masked = toks.clone()
                masked[0,pos+1] = mask_idx
                with torch.no_grad():
                    out = model(masked, repr_layers=[model.num_layers])
                logits = out["logits"][0,pos+1]
                probs = logits.softmax(dim=0).log()
                llr = probs[mt_id] - probs[wt_id]
                h = out["representations"][model.num_layers][0,pos+1]
                contrib = (model.lm_head.weight[mt_id]-model.lm_head.weight[wt_id])*h
                llr_contrib_data.append({
                    "row_index": idx,
                    "llr_score": llr.item(),
                    "llr_contributions": contrib.tolist()
                })

        contribs_df = pd.DataFrame(llr_contrib_data).set_index("row_index")
        out = catalog.join(contribs_df, how="left")

        # expand top-K dims
        contrib_matrix = np.array(out["llr_contributions"].dropna().tolist())
        mean_abs = np.mean(np.abs(contrib_matrix), axis=0)
        top_dims = np.argsort(mean_abs)[::-1][:top_k]
        for dim in top_dims:
            out[f"expanded_llr_dim_{dim}"] = out["llr_contributions"].apply(
                lambda v: v[dim] if isinstance(v,list) else np.nan
            )
        return out

class ProximityTools:
    @staticmethod
    def map_confidence(confidence: str) -> str:
        if confidence in ['1) Assoc w R', '2) Assoc w R - Interim']:
            return "Resistant"
        elif confidence in ['4) Not assoc w R - Interim', '5) Not assoc w R']:
            return "Susceptible"
        else:
            return "Unknown"

    @staticmethod
    def adjust_number(index: int) -> int:
        """Adjust index offset for rpoB (WHO vs PDB numbering)."""
        return index + 6

    @staticmethod
    def proximity_1D(gene_df: pd.DataFrame, r_positions: list, gene: str) -> pd.DataFrame:
        results, nearest = [], []
        for _, row in gene_df.iterrows():
            pos = row['WHO_Adjusted_Position'] if gene == "rpoB" else row['position']
            dists = [(abs(pos - r), r) for r in r_positions if r != pos]
            if dists:
                dist, nearest_pos = min(dists, key=lambda x: x[0])
            else:
                dist, nearest_pos = np.nan, np.nan
            results.append(dist)
            nearest.append(nearest_pos)
        gene_df['Proximity_1D'] = results
        gene_df['Nearest_1D_Index'] = nearest
        return gene_df

    @staticmethod
    def proximity_3D(gene_df: pd.DataFrame, r_positions: list, dist_map, gene: str) -> pd.DataFrame:
        results, nearest = [], []
        for _, row in gene_df.iterrows():
            pos = row['WHO_Adjusted_Position'] if gene == "rpoB" else row['position']
            dists = []
            for r in r_positions:
                if r == pos:
                    continue
                try:
                    d = dist_map.dist(pos, r, raise_na=True)
                    if not np.isnan(d):
                        dists.append((d, r))
                except KeyError:
                    continue
            if dists:
                dist, nearest_pos = min(dists, key=lambda x: x[0])
            else:
                dist, nearest_pos = np.nan, np.nan
            results.append(dist)
            nearest.append(nearest_pos)
        gene_df['Proximity_to_R_Conferring'] = results
        gene_df['Nearest_Mutation_Index'] = nearest
        return gene_df

    @staticmethod
    def add_proximity_features(catalog_df: pd.DataFrame,
                               protein_details_path: str,
                               distmap_dir: str) -> pd.DataFrame:
        """
        Add 1D/3D proximity features only.
        Expects catalog_df to have 'one_letter_mutation' and 'confidence'.
        """
        catalog_df = catalog_df.copy()
        catalog_df['phenotype'] = catalog_df['confidence'].apply(ProximityTools.map_confidence)
        catalog_df['position'] = catalog_df['one_letter_mutation'].str.extract(r'(\d+)').astype(int)

        # Load protein details
        details = pd.read_excel(protein_details_path, sheet_name='Sheet1')
        genes = catalog_df['gene'].unique()
        details = details[details['gene_name'].str.contains('|'.join(genes), case=False, na=False)]

        all_data = []

        for _, row in details.iterrows():
            gene = row['gene_name']
            uniprot = row['Uniprot']
            entry = row['Entry']

            gene_df = catalog_df[catalog_df['gene'] == gene].drop_duplicates()

            # Try loading distmap
            try:
                distmap_path = f"{distmap_dir}/{uniprot}/{entry}"
                dist_map = DistanceMap.from_file(os.path.join(distmap_path))
            except FileNotFoundError:
                print(f" DistMap missing for {gene} ({entry}), filling NaNs.")
                if gene == "rpoB":
                    gene_df["WHO_Adjusted_Position"] = gene_df["position"].apply(ProximityTools.adjust_number)
                else:
                    gene_df["WHO_Adjusted_Position"] = np.nan
                gene_df["Proximity_1D"] = np.nan
                gene_df["Nearest_1D_Index"] = np.nan
                gene_df["Proximity_to_R_Conferring"] = np.nan
                gene_df["Nearest_Mutation_Index"] = np.nan
                all_data.append(gene_df)
                continue

            # Compute R-positions
            if gene == "rpoB":
                gene_df['WHO_Adjusted_Position'] = gene_df['position'].apply(ProximityTools.adjust_number)
                r_positions = gene_df[gene_df['phenotype']=="Resistant"]['WHO_Adjusted_Position'].tolist()
            else:
                r_positions = gene_df[gene_df['phenotype']=="Resistant"]['position'].tolist()

            r_positions = np.unique(r_positions)
            gene_df = ProximityTools.proximity_1D(gene_df, r_positions, gene)
            gene_df = ProximityTools.proximity_3D(gene_df, r_positions, dist_map, gene)

            all_data.append(gene_df)

        final_df = pd.concat(all_data, ignore_index=True)

        # Ensure column always exists
        if "Proximity_to_R_Conferring" not in final_df.columns:
            final_df["Proximity_to_R_Conferring"] = np.nan

        # Zeroing rule
        final_df["Proximity_to_R_Conferring_zeroed"] = final_df["Proximity_to_R_Conferring"]
        for (gene, drug), idxs in final_df.groupby(["gene","drug"]).groups.items():
            sub = final_df.loc[idxs]
            for pos, pidxs in sub.groupby("position").groups.items():
                subdf = sub.loc[pidxs]
                if len(subdf) == 1:
                    continue
                if (subdf['phenotype'] == "Resistant").sum() == 1:
                    final_df.loc[pidxs[1:], "Proximity_to_R_Conferring_zeroed"] = 0
                elif (subdf['phenotype'] == "Resistant").sum() > 1:
                    final_df.loc[pidxs, "Proximity_to_R_Conferring_zeroed"] = 0

        return final_df




class AAIndexTools:
    @staticmethod
    def extract_aa_and_position(mutation: str):
        m = re.match(r"([A-Z])(\d+)([A-Z])", mutation)
        if m:
            return m.groups()
        return None, None, None

    @staticmethod
    def add_rosetta_ddg(catalog_df: pd.DataFrame,
                        rosetta_dir: str,
                        out_dir: Optional[str] = None) -> pd.DataFrame:
        """
        Merge Rosetta ddG scores (thermostability) into catalog_df.
        Handles gyrA Nter/Cter special case.
        """
        all_gene_data = []
        for gene, subset in catalog_df.groupby("gene"):
            gene = str(gene)
            if gene == "gyrA":
                nter = pd.read_csv(os.path.join(rosetta_dir, f"{gene}_Nter_ddG.csv"))
                cter = pd.read_csv(os.path.join(rosetta_dir, f"{gene}_Cter_ddG.csv"))
                rosetta = pd.concat([nter, cter], ignore_index=True)
            else:
                rosetta = pd.read_csv(os.path.join(rosetta_dir, f"{gene}_ddG.csv"))

            rosetta['Wildtype_AA'] = rosetta['variant'].str.extract(r'([A-Za-z])(?=\d)')
            rosetta['position']    = rosetta['variant'].str.extract(r'(\d+)').astype(int)
            rosetta['Mutated_AA']  = rosetta['variant'].str.extract(r'(?<=\d)([A-Za-z])')

            merged = pd.merge(
                rosetta,
                subset,
                left_on=["Wildtype_AA","position","Mutated_AA"],
                right_on=["Wildtype_AA","position","Mutated_AA"],
                how="inner"
            )
            # Drop unnecessary Rosetta cols
            drop_cols = ['fa_intra_rep','fa_intra_sol_xover4','lk_ball_wtd',
                         'pro_close','hbond_sr_bb','hbond_lr_bb',
                         'hbond_bb_sc','hbond_sc','omega','p_aa_pp',
                         'ref','rama_prepro','variant']
            merged = merged.drop(columns=[c for c in drop_cols if c in merged], errors="ignore")
            merged = merged.rename(columns={"score":"thermostability"})
            all_gene_data.append(merged)

        # out = pd.concat(all_gene_data, ignore_index=True).drop_duplicates()
        out = pd.concat(all_gene_data, ignore_index=True)
        
        # Drop problematic list columns before de-duplication
        if "llr_contributions" in out.columns:
            out = out.drop(columns=["llr_contributions"])
        
        out = out.drop_duplicates()

        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
            out.to_csv(os.path.join(out_dir,"catalog_with_thermostability.csv"),index=False)
        return out

    @staticmethod
    def add_aaindex_features(catalog_df: pd.DataFrame,
                             aaindex_path: str,
                             n_cols: int = 8,
                             out_path: Optional[str] = None) -> pd.DataFrame:
        """
        Add AAIndex PCA features: mutant vector and delta (mutant - wildtype).
        """
        aa_df = pd.read_csv(aaindex_path, index_col=0)
        aa_dict = aa_df.set_index(aa_df.index).T.to_dict("list")

        all_data = []
        for gene, subset in catalog_df.groupby("gene"):
            subset = subset.copy()
            muts = subset["one_letter_mutation"]

            mutant_cols = [[] for _ in range(n_cols)]
            diff_cols   = [[] for _ in range(n_cols)]

            for mut in muts:
                wt, pos, mt = AAIndexTools.extract_aa_and_position(mut)
                wt_vec = aa_dict.get(wt, [None]*n_cols)[:n_cols]
                mt_vec = aa_dict.get(mt, [None]*n_cols)[:n_cols]

                for i in range(n_cols):
                    mutant_cols[i].append(mt_vec[i])
                    diff_cols[i].append(
                        mt_vec[i]-wt_vec[i] if (wt_vec[i] is not None and mt_vec[i] is not None) else None
                    )

            for i in range(n_cols):
                subset[f"mut_AAIndex{i+1}"]  = mutant_cols[i]
                subset[f"delta_AAIndex{i+1}"] = diff_cols[i]

            all_data.append(subset)

        out = pd.concat(all_data, ignore_index=True).drop_duplicates()
        if out_path:
            out.to_csv(out_path, index=False)
        return out

class FeaturePipeline:
    @staticmethod
    def run(catalog_df: pd.DataFrame,
            protein_df: pd.DataFrame,
            fasta_dir: str,
            distmap_dir: str,
            protein_details_path: str,
            rosetta_dir: str,
            aaindex_path: str,
            esm_model_path: str,
            out_dir: str,
            top_k_llr: int = 10) -> pd.DataFrame:
        """
        Full feature pipeline:
        - Frequency
        - Delta-Z
        - Proximity (1D, 3D) + phenotype mapping
        - LLR contributions
        - Rosetta ddG (thermostability)
        - AAIndex PCA descriptors
        """
        os.makedirs(out_dir, exist_ok=True)

        # Step 1: Frequency
        freq_df = FeatureTools.add_frequency(catalog_df)
        freq_out = os.path.join(out_dir, "step1_freq.csv")
        freq_df.to_csv(freq_out, index=False)

        # Step 2: Delta-Z
        dz_df = FeatureTools.compute_delta_z(fasta_dir, out_path=os.path.join(out_dir,"step2_delta_z.csv"))
        merged_dz = FeatureTools.merge_delta_z(freq_df, dz_df)
        dz_out = os.path.join(out_dir, "step2_freq_delta_z.csv")
        merged_dz.to_csv(dz_out, index=False)

        # Step 3: Proximity
        prox_df = ProximityTools.add_proximity_features(
            merged_dz,
            protein_details_path=protein_details_path,
            distmap_dir=distmap_dir
        )
        prox_out = os.path.join(out_dir, "step3_proximity.csv")
        prox_df.to_csv(prox_out, index=False)

        # Step 4: LLR contributions
        llr_df = FeatureTools.add_llr_contributions(
            prox_df,
            protein_sequences=protein_df,
            model_path=esm_model_path,
            top_k=top_k_llr
        )
        llr_out = os.path.join(out_dir, "step4_llr.csv")
        llr_df.to_csv(llr_out, index=False)

        # Step 5: Rosetta ddG
        ddg_df = AAIndexTools.add_rosetta_ddg(llr_df, rosetta_dir=rosetta_dir)
        ddg_out = os.path.join(out_dir, "step5_rosetta.csv")
        ddg_df.to_csv(ddg_out, index=False)

        # Step 6: AAIndex PCA
        aaindex_df = AAIndexTools.add_aaindex_features(ddg_df, aaindex_path=aaindex_path, n_cols=8)
        aaindex_out = os.path.join(out_dir, "step6_aaindex.csv")
        aaindex_df.to_csv(aaindex_out, index=False)

        print(f" Feature pipeline finished. Final dataset saved at {aaindex_out}")
        return aaindex_df
