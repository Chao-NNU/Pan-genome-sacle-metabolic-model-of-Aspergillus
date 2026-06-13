#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np

BASE_DIR = "/data/gan_1/aaaResults/FSEOF_FVSEOF"
FBA_DIR = os.path.join(BASE_DIR, "FBA")
OUT_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(OUT_DIR, exist_ok=True)

records = []

for strain in os.listdir(FBA_DIR):
    strain_dir = os.path.join(FBA_DIR, strain)
    if not os.path.isdir(strain_dir):
        continue

    for f in os.listdir(strain_dir):
        if not f.endswith(".xlsx"):
            continue

        acid = f.replace(".xlsx", "")
        df = pd.read_excel(os.path.join(strain_dir, f))
        df["Strain"] = strain
        df["Organic_Acid"] = acid
        records.append(df)

all_df = pd.concat(records, ignore_index=True)

all_df = all_df.replace([np.inf, -np.inf], np.nan)
all_df = all_df.dropna(subset=["WT_Growth", "WT_Production"])

all_df["GR"] = all_df["New_Growth"] / (all_df["WT_Growth"] + 1e-9)
all_df["PF"] = all_df["New_Production"] / (all_df["WT_Production"] + 1e-9)

fig4b_pool = all_df[(all_df["GR"] > 0.6) & (all_df["PF"] > 1.0)]
fig4b_pool.to_excel(
    os.path.join(OUT_DIR, "figure4b_feasible_reactions.xlsx"),
    index=False
)

golden_rows = []

for acid, sub in fig4b_pool.groupby("Organic_Acid"):
    strain_stat = (
        sub.groupby("Strain")
        .agg(
            feasible_rxn_count=("Reaction_ID", "nunique"),
            max_PF=("PF", "max")
        )
        .reset_index()
        .sort_values(
            ["feasible_rxn_count", "max_PF"],
            ascending=False
        )
        .head(5)
    )
    strain_stat["Organic_Acid"] = acid
    golden_rows.append(strain_stat)

golden_df = pd.concat(golden_rows, ignore_index=True)
golden_df.to_excel(
    os.path.join(OUT_DIR, "golden_strains.xlsx"),
    index=False
)

top_rxn_rows = []

for acid, sub in fig4b_pool.groupby("Organic_Acid"):
    gold_strains = golden_df[
        golden_df["Organic_Acid"] == acid
    ]["Strain"]

    rxn_stat = (
        sub[sub["Strain"].isin(gold_strains)]
        .groupby(["Reaction_ID", "Operation_Type"])
        .agg(
            mean_PF=("PF", "mean"),
            mean_GR=("GR", "mean"),
            support_strains=("Strain", "nunique")
        )
        .reset_index()
    )

    rxn_stat = rxn_stat[rxn_stat["support_strains"] >= 3]
    rxn_stat = rxn_stat.sort_values(
        ["mean_PF", "mean_GR"],
        ascending=False
    ).head(30)

    rxn_stat["Organic_Acid"] = acid
    top_rxn_rows.append(rxn_stat)

top30_df = pd.concat(top_rxn_rows, ignore_index=True)
top30_df.to_excel(
    os.path.join(OUT_DIR, "top30_key_reactions.xlsx"),
    index=False
)

print("Step 1 complete: golden strains & top30 reactions identified.")
