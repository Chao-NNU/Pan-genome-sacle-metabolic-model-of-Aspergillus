#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import cobra
import pandas as pd

BASE_DIR = "/data/gan_1/aaaResults/FSEOF_FVSEOF"
MODEL_DIR = "/data/gan_1/add_sbml_models"
IN_DIR = os.path.join(BASE_DIR, "results")
OUT_DIR = os.path.join(BASE_DIR, "results/1")
os.makedirs(OUT_DIR, exist_ok=True)

golden_df = pd.read_excel(os.path.join(IN_DIR, "golden_strains.xlsx"))
rxn_df = pd.read_excel(os.path.join(IN_DIR, "top30_key_reactions.xlsx"))

print(f"Golden strains shape: {golden_df.shape}")
print(f"Reactions df shape: {rxn_df.shape}")
print("\nSample from reactions df:")
print(rxn_df.head())

ACID_TO_EXCHANGE = {
    'Citrate': 'UF02577_E',
    'Fumarate': 'UF02628_E',
    'Pyruvate': 'UF02776_E',
    'Succinate': 'UF02804_E',
    'Butyrate': 'UF02912_E',
    'Hexanoate': 'UF02918_E',
    'Lactate': 'UF02967_E',
    'Isocitrate': 'UF02970_E',
    'Oxaloacetate': 'UF03065_E',
    'Acetate': 'UF03160_E',
    'Acetoin': 'UF03161_E',
    'Formate': 'UF03271_E',
    'Glyoxylate': 'UF03295_E',
    'Propionate': 'UF03420_E',
    'Malate': 'UF03510_E'
}

GRADIENTS = {
    "up": [1.5, 2.0, 3.0],
    "down": [0.7, 0.4, 0.2]
}

records = []

processed_count = 0
skipped_count = 0

for idx, r in rxn_df.iterrows():
    acid_name = r["Organic_Acid"]
    rxn_id = r["Reaction_ID"]
    op = r["Operation_Type"]
    
    print(f"\nProcessing row {idx}: Acid={acid_name}, Reaction={rxn_id}, Operation={op}")

    strains = golden_df[
        golden_df["Organic_Acid"] == acid_name
    ]["Strain"].tolist()
    
    print(f"Found {len(strains)} strains for {acid_name}")
    
    if acid_name in ACID_TO_EXCHANGE:
        exchange_id = ACID_TO_EXCHANGE[acid_name]
        print(f"Exchange reaction ID for {acid_name}: {exchange_id}")
    else:
        print(f"Warning: No exchange reaction mapping found for {acid_name}")
        skipped_count += 1
        continue
    
    for strain in strains:
        model_path = os.path.join(MODEL_DIR, f"{strain}.xml")
        
        if not os.path.exists(model_path):
            print(f"  Skipping {strain}: Model file not found at {model_path}")
            skipped_count += 1
            continue
            
        print(f"  Processing strain: {strain}")
        
        try:
            model = cobra.io.read_sbml_model(model_path)
            
            print(f"    Model reactions count: {len(model.reactions)}")
            print(f"    Checking if reaction {rxn_id} exists in model...")
            
            if rxn_id not in model.reactions:
                print(f"    Skipping: Reaction {rxn_id} not found in model")
                skipped_count += 1
                continue
            
            if exchange_id not in model.reactions:
                print(f"    Skipping: Exchange reaction {exchange_id} not found in model")
                potential_exchanges = [rxn.id for rxn in model.exchanges 
                                      if acid_name.lower() in rxn.name.lower() 
                                      or acid_name.lower() in rxn.id.lower()]
                if potential_exchanges:
                    print(f"    But found potential exchange reactions: {potential_exchanges}")
                skipped_count += 1
                continue
            
            target = model.reactions.get_by_id(rxn_id)
            prod_rxn = model.reactions.get_by_id(exchange_id)
            
            if "BIOMASS" in model.reactions:
                biomass = model.reactions.get_by_id("BIOMASS")
                biomass.lower_bound = 0.01
            else:
                biomass_rxns = [rxn.id for rxn in model.reactions 
                               if "biomass" in rxn.id.lower() or "BIOMASS" in rxn.id]
                if biomass_rxns:
                    biomass = model.reactions.get_by_id(biomass_rxns[0])
                    biomass.lower_bound = 0.01
                else:
                    print(f"    Warning: No biomass reaction found in model")
                    skipped_count += 1
                    continue

            try:
                wt = model.optimize()
                if wt.status != 'optimal':
                    print(f"    Skipping: Wild-type optimization failed with status {wt.status}")
                    skipped_count += 1
                    continue
                    
                wt_prod = wt.fluxes[prod_rxn.id]
                print(f"    Wild-type {acid_name} production: {wt_prod:.4f}")
                
            except Exception as e:
                print(f"    Error during wild-type optimization: {str(e)}")
                skipped_count += 1
                continue

            for g in GRADIENTS[op]:
                with model:
                    original_lb = target.lower_bound
                    original_ub = target.upper_bound
                    
                    print(f"    Applying gradient {g} to {rxn_id} ({op})")
                    print(f"    Original bounds: [{original_lb}, {original_ub}]")
                    
                    if op == "up":
                        target.lower_bound = original_lb * g
                        target.upper_bound = original_ub * g
                    else:  # down
                        new_lb = original_lb * g
                        new_ub = original_ub * g
                        if abs(new_lb) < 1e-6:
                            new_lb = 0
                        if abs(new_ub) < 1e-6:
                            new_ub = 0
                        target.lower_bound = new_lb
                        target.upper_bound = new_ub
                    
                    print(f"    New bounds: [{target.lower_bound}, {target.upper_bound}]")

                    sol = model.optimize()
                    
                    if sol.status != 'optimal':
                        print(f"    Gradient {g}: Optimization failed with status {sol.status}")
                        continue

                    uptake = 0
                    for rxn in model.exchanges:
                        if rxn.flux is not None and rxn.flux < 0:
                            uptake += abs(rxn.flux)
                    
                    prod_flux = sol.fluxes[prod_rxn.id]
                    yield_value = prod_flux / uptake if uptake > 0 else 0
                    
                    biomass_flux = sol.fluxes.get("BIOMASS", 0)
                    
                    print(f"    Gradient {g}: Prod={prod_flux:.4f}, Uptake={uptake:.4f}, Yield={yield_value:.4f}, Biomass={biomass_flux:.4f}")
                    
                    records.append({
                        "Organic_Acid": acid_name,
                        "Strain": strain,
                        "Reaction_ID": rxn_id,
                        "Exchange_Reaction": exchange_id,
                        "Operation": op,
                        "Gradient": g,
                        "Product_Flux": prod_flux,
                        "WT_Product_Flux": wt_prod,
                        "Yield": yield_value,
                        "Uptake": uptake,
                        "Biomass_Flux": biomass_flux,
                        "Target_Lower_Bound": target.lower_bound,
                        "Target_Upper_Bound": target.upper_bound
                    })
                    
                    processed_count += 1
                    
        except Exception as e:
            print(f"    Error processing {strain}: {str(e)}")
            import traceback
            traceback.print_exc()
            skipped_count += 1

print(f"\n=== Summary ===")
print(f"Total processed cases: {processed_count}")
print(f"Total skipped cases: {skipped_count}")

if records:
    df = pd.DataFrame(records)
    output_path = os.path.join(OUT_DIR, "gradient_yield_validation.xlsx")
    df.to_excel(output_path, index=False)
    print(f"\nResults saved to: {output_path}")
    print(f"Output shape: {df.shape}")
    
    print("\n=== Results Summary ===")
    print(f"Unique acids processed: {df['Organic_Acid'].nunique()}")
    print(f"Unique strains processed: {df['Strain'].nunique()}")
    print(f"Unique reactions processed: {df['Reaction_ID'].nunique()}")
    
    print("\nResults by organic acid:")
    acid_stats = df.groupby('Organic_Acid').size().reset_index(name='Count')
    print(acid_stats)
    
    print("\nFirst few rows of results:")
    print(df.head())
else:
    print("\nNo records were generated!")
    print("Possible issues to check:")
    print("1. Model files exist in MODEL_DIR")
    print("2. Reaction IDs in top30_key_reactions.xlsx match those in models")
    print("3. Strain names in golden_strains.xlsx match model filenames")
    print("4. Organic acid names match those in ACID_TO_EXCHANGE mapping")
    
print("\nStep 2 complete: gradient validation finished.")
