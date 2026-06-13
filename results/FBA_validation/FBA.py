import os
import cobra
import pandas as pd
import logging
from pathlib import Path

# --- 1. Path and Parameter Configuration ---
BASE_DIR = "/data/gan_1/aaaResults/FSEOF_FVSEOF"
MODEL_DIR = "/data/gan_1/add_sbml_models"
INPUT_DATA_DIR = os.path.join(BASE_DIR, "441models")
CLUSTER_FILE = os.path.join(BASE_DIR, "cluster_assignments.xlsx")
OUTPUT_BASE_DIR = os.path.join(BASE_DIR, "FBA")

os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)

# Organic Acid Mapping
ORGANIC_ACID_MAP = {
    'UF02577_E': 'Citrate', 'UF02628_E': 'Fumarate', 'UF02776_E': 'Pyruvate',
    'UF02804_E': 'Succinate', 'UF02912_E': 'Butyrate', 'UF02918_E': 'Hexanoate',
    'UF02967_E': 'Lactate', 'UF02970_E': 'Isocitrate', 'UF03065_E': 'Oxaloacetate',
    'UF03160_E': 'Acetate', 'UF03161_E': 'Acetoin', 'UF03271_E': 'Formate',
    'UF03295_E': 'Glyoxylate', 'UF03420_E': 'Propionate', 'UF03510_E': 'Malate'
}

BIOMASS_ID = "BIOMASS"
GRADIENTS = {'up': [2.0, 5.0, 10.0], 'down': [0.0, 0.3, 0.6]}

logging.getLogger("cobra").setLevel(logging.ERROR)

# ---------- Utility Functions ----------

def is_finished(output_file):
    """Check whether a result file already exists and is non-empty."""
    return os.path.exists(output_file) and os.path.getsize(output_file) > 0

def calculate_metrics(model, product_id, biomass_min=0.01):
    """Calculate growth rate and production rate."""
    model.objective = BIOMASS_ID
    sol_growth = model.optimize()
    growth = sol_growth.objective_value if sol_growth.status == 'optimal' else 0.0

    production = 0.0
    if growth > 1e-6:
        with model:
            model.reactions.get_by_id(BIOMASS_ID).lower_bound = biomass_min
            model.objective = product_id
            sol_prod = model.optimize()
            if sol_prod.status == 'optimal':
                production = sol_prod.objective_value

    return growth, production

# ---------- Main Workflow ----------

def main():
    print(f"Reading cluster assignments from: {CLUSTER_FILE}")
    cluster_df = pd.read_excel(CLUSTER_FILE)

    target_strains = cluster_df[cluster_df['Cluster'].isin([1, 3])]
    strain_list = target_strains['Model'].tolist()

    print(f"Total models identified for Cluster 1 & 3: {len(strain_list)}")

    for i, strain_name in enumerate(strain_list):
        model_path = os.path.join(MODEL_DIR, f"{strain_name}.xml")

        if not os.path.exists(model_path):
            print(f"Skip: {strain_name}.xml not found.")
            continue

        print(f"\n[{i+1}/{len(strain_list)}] Loading model: {strain_name}")
        try:
            model = cobra.io.read_sbml_model(model_path)
        except Exception as e:
            print(f"Error loading {strain_name}: {e}")
            continue

        strain_output_dir = os.path.join(OUTPUT_BASE_DIR, strain_name)
        os.makedirs(strain_output_dir, exist_ok=True)

        for acid_id, acid_name in ORGANIC_ACID_MAP.items():
            if acid_id not in model.reactions:
                continue

            fseof_path = os.path.join(INPUT_DATA_DIR, strain_name, f"{acid_id}.xlsx")
            if not os.path.exists(fseof_path):
                continue

            output_file = os.path.join(strain_output_dir, f"{acid_name}.xlsx")
            if is_finished(output_file):
                print(f"  Skip (already finished): {acid_name}")
                continue

            print(f"  Processing Target: {acid_name}")

            wt_growth, wt_prod = calculate_metrics(model, acid_id)
            validation_data = []

            for op_type in ['up', 'down']:
                try:
                    targets_df = pd.read_excel(fseof_path, sheet_name=op_type)
                except Exception:
                    continue

                for _, row in targets_df.iterrows():
                    rxn_id = str(row['reaction']).strip()
                    if rxn_id not in model.reactions:
                        continue

                    model.objective = BIOMASS_ID
                    base_flux = model.optimize().fluxes.get(rxn_id, 0.0)

                    for f in GRADIENTS[op_type]:
                        with model:
                            target_rxn = model.reactions.get_by_id(rxn_id)

                            if op_type == 'down':
                                limit = abs(base_flux) * f
                                target_rxn.bounds = (
                                    -limit if target_rxn.lower_bound < 0 else 0,
                                    limit
                                )
                            else:
                                boost = max(0.1, abs(base_flux)) * f
                                if base_flux >= 0:
                                    target_rxn.lower_bound = boost
                                else:
                                    target_rxn.upper_bound = -boost

                            new_growth, new_prod = calculate_metrics(model, acid_id)

                            validation_data.append({
                                'Reaction_ID': rxn_id,
                                'Operation_Type': op_type,
                                'Factor': f,
                                'Base_Flux': base_flux,
                                'WT_Growth': wt_growth,
                                'New_Growth': new_growth,
                                'Growth_Change_Rate': (
                                    (new_growth - wt_growth) / wt_growth if wt_growth > 0 else 0
                                ),
                                'WT_Production': wt_prod,
                                'New_Production': new_prod,
                                'Production_Change_Rate': (
                                    (new_prod - wt_prod) / wt_prod if wt_prod > 0 else 0
                                )
                            })

            if validation_data:
                pd.DataFrame(validation_data).to_excel(output_file, index=False)
                print(f"    Saved: {output_file}")

    print(f"\nBatch processing complete. Results stored in: {OUTPUT_BASE_DIR}")

if __name__ == "__main__":
    main()
