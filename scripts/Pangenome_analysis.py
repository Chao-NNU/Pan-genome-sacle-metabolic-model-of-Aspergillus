#!/usr/bin/env python3
"""
Combined pangenome analysis workflow without plotting code.

This script merges the original two analysis stages:
1. Gene family parsing and classification from OrthoFinder outputs.
2. COG annotation mapping and pangenome curve analysis.

Plotting code and plotting-library dependencies have been removed. The core
analysis logic and output directories are preserved.
"""

import os
import sys
import logging
import warnings
from pathlib import Path
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

warnings.filterwarnings("ignore")


# Configuration
BASE_DIR = Path("/data/gan_1/aaaResults/characteristic/genome")
ORTHOFINDER_RESULTS_DIR = BASE_DIR
OUTPUT_DIR = BASE_DIR
RESULTS_DIR = BASE_DIR / "results"
PLOTS_DIR = BASE_DIR / "plots"
ANNOTATION_DIR = Path("/data/gan_1/aaaResults/carvefungi/annotation")
FASTA_DIR = Path("/data/gan_1/aaaResults/carvefungi/fasta")
COG_RESULTS_DIR = BASE_DIR / "cog_results"

for dir_path in [OUTPUT_DIR, RESULTS_DIR, PLOTS_DIR, COG_RESULTS_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

CORE_THRESHOLD = 0.99
SOFT_CORE_THRESHOLD = 0.95
SHELL_THRESHOLD = 0.15
CLOUD_THRESHOLD = 0.15

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(BASE_DIR / "combined_analysis.log", mode="a"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


COG_CATEGORIES = {'J': 'Translation, ribosomal structure and biogenesis', 'K': 'Transcription', 'L': 'Replication, recombination and repair', 'D': 'Cell cycle control, cell division, chromosome partitioning', 'O': 'Posttranslational modification, protein turnover, chaperones', 'M': 'Cell wall/membrane/envelope biogenesis', 'N': 'Cell motility', 'U': 'Intracellular trafficking, secretion, and vesicular transport', 'T': 'Signal transduction mechanisms', 'C': 'Energy production and conversion', 'E': 'Amino acid transport and metabolism', 'F': 'Nucleotide transport and metabolism', 'G': 'Carbohydrate transport and metabolism', 'H': 'Coenzyme transport and metabolism', 'I': 'Lipid transport and metabolism', 'P': 'Inorganic ion transport and metabolism', 'Q': 'Secondary metabolites biosynthesis, transport and catabolism', 'R': 'General function prediction only', 'S': 'Function unknown', 'V': 'Defense mechanisms', 'W': 'Extracellular structures', 'Y': 'Nuclear structure', 'Z': 'Cytoskeleton'}


class OrthoFinderParser:

    def __init__(self, orthofinder_dir):
        self.orthofinder_dir = Path(orthofinder_dir)
        self.strain_names = []

    def find_orthogroups_file(self):
        possible_paths = [self.orthofinder_dir / 'Orthogroups.tsv', self.orthofinder_dir / 'Orthogroups' / 'Orthogroups.tsv', self.orthofinder_dir / 'Results_OrthoFinder' / 'Orthogroups' / 'Orthogroups.tsv']
        for path in possible_paths:
            if path.exists():
                logger.info(f'Found Orthogroups file: {path}')
                return path
        for root, dirs, files in os.walk(self.orthofinder_dir):
            for file in files:
                if file == 'Orthogroups.tsv':
                    found_path = Path(root) / file
                    logger.info(f'Found Orthogroups file: {found_path}')
                    return found_path
        logger.error('Orthogroups.tsv file not found')
        return None

    def find_species_tree(self):
        possible_paths = [self.orthofinder_dir / 'SpeciesTree_rooted.txt.contree', self.orthofinder_dir / 'SpeciesTree_rooted.txt', self.orthofinder_dir / 'Species_Tree' / 'SpeciesTree_rooted.txt']
        for path in possible_paths:
            if path.exists():
                logger.info(f'Found species tree file: {path}')
                return path
        logger.warning('Species tree file not found')
        return None

    def parse_orthogroups_file(self, orthogroups_path):
        logger.info('Parsing Orthogroups.tsv file...')
        try:
            df = pd.read_csv(orthogroups_path, sep='\t')
            logger.info(f'Orthogroups data shape: {df.shape}')
            logger.info(f'Column names: {list(df.columns)}')
            self.strain_names = list(df.columns[1:])
            logger.info(f'Detected {len(self.strain_names)} strains from the Orthogroups file')
            logger.info(f'First 10 strains: {self.strain_names[:10]}')
            df.to_csv(RESULTS_DIR / 'orthogroups_raw.tsv', sep='\t', index=False)
            return df
        except Exception as e:
            logger.error(f'Failed to read Orthogroups file: {str(e)}')
            return None

    def create_presence_absence_matrix(self, orthogroups_df):
        logger.info('Creating gene family presence/absence matrix...')
        try:
            df = orthogroups_df.copy()
            orthogroup_col = df.columns[0]
            presence_matrix = []
            orthogroup_ids = df[orthogroup_col].tolist()
            strain_columns = []
            for col in df.columns[1:]:
                strain_name = col.strip()
                strain_columns.append(strain_name)
                presence_row = []
                for genes in df[col]:
                    if pd.isna(genes) or str(genes).strip() == '':
                        presence_row.append(0)
                    else:
                        gene_str = str(genes)
                        if ', ' in gene_str:
                            gene_count = len(gene_str.split(', '))
                        else:
                            gene_count = len(gene_str.split(','))
                        presence_row.append(gene_count)
                presence_matrix.append(presence_row)
            presence_df = pd.DataFrame(presence_matrix, columns=orthogroup_ids, index=strain_columns).T
            logger.info(f'Presence/absence matrix shape: {presence_df.shape}')
            logger.info(f'Total gene families: {presence_df.shape[0]}')
            logger.info(f'Average gene families per strain: {presence_df.sum(axis=0).mean():.0f}')
            presence_df.to_csv(RESULTS_DIR / 'gene_family_presence_absence.tsv', sep='\t')
            binary_df = (presence_df > 0).astype(int)
            binary_df.to_csv(RESULTS_DIR / 'gene_family_binary_matrix.tsv', sep='\t')
            return (presence_df, binary_df)
        except Exception as e:
            logger.error(f'Failed to create presence/absence matrix: {str(e)}')
            import traceback
            traceback.print_exc()
            return (None, None)

    def classify_gene_families(self, binary_df):
        logger.info('Classifying gene families...')
        try:
            n_strains = binary_df.shape[1]
            logger.info(f'Total number of strains: {n_strains}')
            classifications = []
            counts = defaultdict(int)
            for og_id, row in binary_df.iterrows():
                presence_count = row.sum()
                presence_ratio = presence_count / n_strains
                if presence_ratio >= CORE_THRESHOLD:
                    classification = 'Core'
                elif presence_ratio >= SOFT_CORE_THRESHOLD:
                    classification = 'Soft Core'
                elif presence_ratio >= SHELL_THRESHOLD:
                    classification = 'Shell'
                elif presence_ratio > 0:
                    classification = 'Cloud'
                else:
                    classification = 'Absent'
                if presence_count == 1:
                    classification = 'Private'
                classifications.append(classification)
                counts[classification] += 1
            classified_df = binary_df.copy()
            classified_df['Classification'] = classifications
            classified_df['Presence_Count'] = classified_df.iloc[:, :-1].sum(axis=1)
            classified_df['Presence_Ratio'] = classified_df['Presence_Count'] / n_strains
            classified_df.to_csv(RESULTS_DIR / 'gene_families_classified.tsv', sep='\t')
            total_families = len(classified_df)
            stats_data = []
            for class_type in ['Core', 'Soft Core', 'Shell', 'Cloud', 'Private']:
                if class_type in counts:
                    count = counts[class_type]
                    percentage = count / total_families * 100
                    stats_data.append({'Classification': class_type, 'Count': count, 'Percentage': round(percentage, 2)})
            stats_df = pd.DataFrame(stats_data)
            stats_df = stats_df.sort_values('Percentage', ascending=False)
            stats_df.to_csv(RESULTS_DIR / 'gene_classification_statistics.csv', index=False)
            logger.info('Gene family classification statistics:')
            for _, row in stats_df.iterrows():
                logger.info(f"  {row['Classification']}: {row['Count']} ({row['Percentage']}%)")
            self.calculate_strain_gene_distribution(classified_df)
            return (classified_df, stats_df)
        except Exception as e:
            logger.error(f'Failed to classify gene families: {str(e)}')
            import traceback
            traceback.print_exc()
            raise

    def calculate_strain_gene_distribution(self, classified_df):
        logger.info('Calculating gene family distribution for each strain...')
        core_mask = classified_df['Classification'].isin(['Core', 'Soft Core'])
        core_families = classified_df[core_mask].index.tolist()
        accessory_mask = classified_df['Classification'].isin(['Shell', 'Cloud'])
        accessory_families = classified_df[accessory_mask].index.tolist()
        private_mask = classified_df['Classification'] == 'Private'
        private_families = classified_df[private_mask].index.tolist()
        strain_stats = []
        for strain in classified_df.columns[:-3]:
            core_count = classified_df.loc[core_families, strain].sum() if core_families else 0
            accessory_count = classified_df.loc[accessory_families, strain].sum() if accessory_families else 0
            private_count = classified_df.loc[private_families, strain].sum() if private_families else 0
            strain_stats.append({'Strain': strain, 'Core_Genes': core_count, 'Accessory_Genes': accessory_count, 'Private_Genes': private_count, 'Total_Genes': core_count + accessory_count + private_count})
        strain_stats_df = pd.DataFrame(strain_stats)
        strain_stats_df.to_csv(RESULTS_DIR / 'strain_gene_distribution.csv', index=False)
        logger.info(f"Strain gene distribution statistics saved to: {RESULTS_DIR / 'strain_gene_distribution.csv'}")


class COGAnalyzer:

    def __init__(self, annotation_dir):
        self.annotation_dir = Path(annotation_dir)
        self.cog_data = {}
        self.gene_to_cog = {}

    def parse_emapper_annotations(self):
        logger.info('Parsing emapper annotation files...')
        annotation_files = list(self.annotation_dir.glob('*.faa.emapper.annotations'))
        logger.info(f'Found {len(annotation_files)} annotation files')
        if len(annotation_files) == 0:
            logger.error('No annotation files found')
            return False
        total_genes = 0
        cog_counts = Counter()
        for file_path in annotation_files:
            strain_name = file_path.stem.replace('.faa.emapper.annotations', '')
            logger.debug(f'Processing strain: {strain_name}')
            try:
                df = pd.read_csv(file_path, sep='\t', comment='#', header=None, low_memory=False)
                strain_cogs = {}
                genes_with_cog = 0
                for idx, row in df.iterrows():
                    if len(row) < 6:
                        continue
                    gene_id = str(row[0]).strip()
                    cog_category_cell = str(row[5]).strip()
                    cog_category = self._extract_cog_category(cog_category_cell)
                    if cog_category:
                        strain_cogs[gene_id] = cog_category
                        self.gene_to_cog[gene_id] = cog_category
                        cog_counts[cog_category] += 1
                        genes_with_cog += 1
                    elif len(row) > 4:
                        fallback_info = str(row[4])
                        cog_category = self._extract_cog_category_fallback(fallback_info)
                        if cog_category:
                            strain_cogs[gene_id] = cog_category
                            self.gene_to_cog[gene_id] = cog_category
                            cog_counts[cog_category] += 1
                            genes_with_cog += 1
                self.cog_data[strain_name] = strain_cogs
                total_genes += genes_with_cog
                logger.debug(f'  Strain {strain_name}: {genes_with_cog} genes have COG annotations')
            except Exception as e:
                logger.error(f'Failed to parse file {file_path}: {str(e)}')
                logger.error(f"Example file row: {(df.iloc[0].tolist()[:10] if not df.empty else 'empty file')}")
                continue
        logger.info(f'Parsed a total of {total_genes} genes with COG annotations')
        logger.info(f'COG category distribution: {dict(cog_counts.most_common(10))}')
        cog_stats = pd.DataFrame(cog_counts.items(), columns=['COG_Category', 'Count'])
        cog_stats['Percentage'] = cog_stats['Count'] / cog_stats['Count'].sum() * 100
        cog_stats = cog_stats.sort_values('Count', ascending=False)
        cog_stats.to_csv(COG_RESULTS_DIR / 'cog_overall_statistics.csv', index=False)
        return True

    def _extract_cog_category(self, cog_cell):
        if pd.isna(cog_cell) or not cog_cell:
            return ''
        cog_cell = str(cog_cell).strip()
        if len(cog_cell) == 1 and cog_cell.isalpha() and cog_cell.isupper():
            return cog_cell
        for char in cog_cell:
            if char.isalpha() and char.isupper():
                return char
        return ''

    def _extract_cog_category_fallback(self, cog_info):
        if pd.isna(cog_info) or not cog_info:
            return ''
        import re
        match = re.search('([A-Z])OG\\d+', cog_info)
        if match:
            potential_cog = match.group(1)
            if potential_cog in COG_CATEGORIES:
                return potential_cog
        return ''

    def map_cog_to_gene_families(self, gene_families_file):
        logger.info('Mapping COG annotations to gene families...')
        gene_families_df = pd.read_csv(gene_families_file, sep='\t', index_col=0)
        orthogroups_file = RESULTS_DIR / 'orthogroups_raw.tsv'
        if not orthogroups_file.exists():
            logger.error(f'Orthogroups file not found: {orthogroups_file}')
            return None
        orthogroups_df = pd.read_csv(orthogroups_file, sep='\t')
        gene_to_orthogroup = {}
        for idx, row in orthogroups_df.iterrows():
            og_id = row.iloc[0]
            for strain_col in orthogroups_df.columns[1:]:
                genes = str(row[strain_col])
                if pd.isna(genes) or genes == '':
                    continue
                gene_list = [g.strip() for g in genes.split(', ') if g.strip()]
                for gene in gene_list:
                    gene_to_orthogroup[gene] = og_id
        logger.info(f'Created {len(gene_to_orthogroup)} gene-to-orthogroup mappings')
        logger.info(f'Number of genes with available COG annotations: {len(self.gene_to_cog)}')
        og_cog_mapping = {}
        og_cog_counts = {}
        genes_with_mapped_cog = 0
        for og_id in gene_families_df.index:
            cog_list = []
            for gene, gene_og in gene_to_orthogroup.items():
                if gene_og == og_id and gene in self.gene_to_cog:
                    cog_list.append(self.gene_to_cog[gene])
            if cog_list:
                genes_with_mapped_cog += len(cog_list)
                cog_counter = Counter(cog_list)
                most_common_cog = cog_counter.most_common(1)[0][0]
                og_cog_mapping[og_id] = most_common_cog
                og_cog_counts[og_id] = dict(cog_counter)
            else:
                og_cog_mapping[og_id] = 'S'
                og_cog_counts[og_id] = {'S': 1}
        logger.info(f'Assigned COG categories to {len(og_cog_mapping)} gene families')
        logger.info(f'Total mapped genes: {genes_with_mapped_cog}')
        gene_families_with_cog = gene_families_df.copy()
        gene_families_with_cog['COG_Category'] = gene_families_with_cog.index.map(lambda x: og_cog_mapping.get(x, 'S'))
        gene_families_with_cog.to_csv(RESULTS_DIR / 'gene_families_with_cog.tsv', sep='\t')
        og_cog_df = pd.DataFrame.from_dict(og_cog_counts, orient='index')
        og_cog_df.fillna(0, inplace=True)
        og_cog_df.to_csv(COG_RESULTS_DIR / 'orthogroup_cog_details.tsv', sep='\t')
        mapping_stats = {'total_orthogroups': len(gene_families_df), 'orthogroups_with_cog': sum((1 for v in og_cog_mapping.values() if v != 'S')), 'genes_in_orthogroups': len(gene_to_orthogroup), 'genes_with_cog_annotation': len(self.gene_to_cog), 'genes_mapped_to_orthogroups': genes_with_mapped_cog}
        mapping_stats_df = pd.DataFrame([mapping_stats])
        mapping_stats_df.to_csv(COG_RESULTS_DIR / 'cog_mapping_statistics.csv', index=False)
        return gene_families_with_cog

    def analyze_cog_by_gene_class(self, gene_families_with_cog):
        logger.info('Analyzing COG distribution by gene class...')
        df = gene_families_with_cog
        if 'Classification' not in df.columns:
            logger.error('The gene family table does not contain a Classification column')
            return None
        classifications = df['Classification'].unique()
        cog_by_class = {}
        for class_type in classifications:
            class_df = df[df['Classification'] == class_type]
            cog_dist = class_df['COG_Category'].value_counts()
            cog_by_class[class_type] = cog_dist
        all_cogs = sorted(set(df['COG_Category'].unique()))
        percentage_data = []
        for cog in all_cogs:
            if cog not in COG_CATEGORIES:
                continue
            row = {'COG_Category': cog, 'COG_Description': COG_CATEGORIES.get(cog, 'Unknown')}
            for class_type in ['Core', 'Accessory', 'Private']:
                if class_type in cog_by_class and cog in cog_by_class[class_type]:
                    count = cog_by_class[class_type][cog]
                    total = len(df[df['Classification'] == class_type])
                    percentage = count / total * 100
                else:
                    count = 0
                    percentage = 0.0
                row[class_type] = percentage
                row[f'{class_type}_count'] = count
            percentage_data.append(row)
        percentage_df = pd.DataFrame(percentage_data)
        percentage_df = percentage_df.sort_values('Core', ascending=False)
        percentage_df.to_csv(COG_RESULTS_DIR / 'cog_percentage_by_gene_class.csv', index=False)
        count_data = []
        for class_type in classifications:
            if class_type in cog_by_class:
                for cog, count in cog_by_class[class_type].items():
                    count_data.append({'Classification': class_type, 'COG_Category': cog, 'Count': count, 'Percentage': count / len(df[df['Classification'] == class_type]) * 100})
        count_df = pd.DataFrame(count_data)
        count_df.to_csv(COG_RESULTS_DIR / 'cog_counts_by_gene_class.csv', index=False)
        logger.info('COG analysis by class completed')
        return (percentage_df, count_df)


class PangenomeCurveAnalyzer:

    def __init__(self, binary_matrix_file):
        self.binary_matrix_file = Path(binary_matrix_file)
        self.binary_matrix = None

    def load_binary_matrix(self):
        logger.info('Loading gene family binary matrix...')
        if not self.binary_matrix_file.exists():
            logger.error(f'Binary matrix file does not exist: {self.binary_matrix_file}')
            return False
        self.binary_matrix = pd.read_csv(self.binary_matrix_file, sep='\t', index_col=0)
        logger.info(f'Matrix shape: {self.binary_matrix.shape}')
        return True

    def calculate_pangenome_curve(self, n_permutations=100, step_size=10):
        logger.info('Calculating pangenome curves...')
        if self.binary_matrix is None:
            logger.error('Binary matrix has not been loaded')
            return None
        n_strains = self.binary_matrix.shape[1]
        n_orthogroups = self.binary_matrix.shape[0]
        logger.info(f'Number of strains: {n_strains}, gene families: {n_orthogroups}')
        if n_strains <= 50:
            steps = list(range(1, n_strains + 1))
        else:
            steps = list(range(1, min(51, n_strains), 1)) + list(range(55, n_strains + 1, step_size))
            if n_strains not in steps:
                steps.append(n_strains)
        pan_means = []
        pan_stds = []
        core_means = []
        core_stds = []
        strain_counts = []
        for n in steps:
            logger.debug(f'Calculating curve point for {n} strains...')
            pan_sizes = []
            core_sizes = []
            for perm in range(n_permutations):
                selected_strains = np.random.choice(self.binary_matrix.columns, size=min(n, n_strains), replace=False)
                selected_data = self.binary_matrix[selected_strains]
                pan_size = (selected_data.sum(axis=1) > 0).sum()
                core_size = (selected_data.sum(axis=1) == len(selected_strains)).sum()
                pan_sizes.append(pan_size)
                core_sizes.append(core_size)
            pan_means.append(np.mean(pan_sizes))
            pan_stds.append(np.std(pan_sizes))
            core_means.append(np.mean(core_sizes))
            core_stds.append(np.std(core_sizes))
            strain_counts.append(n)
        curve_data = pd.DataFrame({'Strain_Count': strain_counts, 'Pan_Genome_Mean': pan_means, 'Pan_Genome_Std': pan_stds, 'Core_Genome_Mean': core_means, 'Core_Genome_Std': core_stds})
        curve_data.to_csv(RESULTS_DIR / 'pangenome_curve_data.csv', index=False)
        logger.info('Pangenome curve calculation completed')
        return curve_data

    def fit_pangenome_models(self, curve_data):
        logger.info('Fitting pangenome models...')
        x = curve_data['Strain_Count'].values
        y_pan = curve_data['Pan_Genome_Mean'].values
        y_core = curve_data['Core_Genome_Mean'].values

        def heaps_law(N, k, gamma):
            return k * N ** gamma
        try:
            popt_pan, pcov_pan = curve_fit(heaps_law, x, y_pan, p0=[1000, 0.5], bounds=([0, 0], [np.inf, 1]))
            k_pan, gamma_pan = popt_pan
            pan_r2 = self._calculate_r2(y_pan, heaps_law(x, *popt_pan))
            if gamma_pan > 0.1:
                openness = 'Open'
            elif gamma_pan > 0.05:
                openness = 'Likely Open'
            else:
                openness = 'Closed'

            def exp_decay(N, a, b, c):
                return a * np.exp(-b * N) + c
            popt_core, pcov_core = curve_fit(exp_decay, x, y_core, p0=[1000, 0.01, 1000])
            a_core, b_core, c_core = popt_core
            core_r2 = self._calculate_r2(y_core, exp_decay(x, *popt_core))
            fit_results = {'Pan_genome_k': k_pan, 'Pan_genome_gamma': gamma_pan, 'Pan_genome_R2': pan_r2, 'Pan_genome_openness': openness, 'Core_genome_a': a_core, 'Core_genome_b': b_core, 'Core_genome_c': c_core, 'Core_genome_R2': core_r2, 'Final_pan_size': y_pan[-1], 'Final_core_size': y_core[-1], 'Total_strains': x[-1]}
            fit_df = pd.DataFrame([fit_results])
            fit_df.to_csv(RESULTS_DIR / 'pangenome_fit_results.csv', index=False)
            logger.info(f'Pangenome fitting result: gamma = {gamma_pan:.4f} ({openness})')
            logger.info(f'Final pan-genome size: {y_pan[-1]:.0f}')
            logger.info(f'Final core-genome size: {y_core[-1]:.0f}')
            return fit_results
        except Exception as e:
            logger.error(f'Model fitting failed: {str(e)}')
            return None

    def _calculate_r2(self, y_true, y_pred):
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        return 1 - ss_res / ss_tot if ss_tot != 0 else 0



def generate_part1_summary_report(parser, stats_df, presence_df, output_path):
    logger.info("Generating part 1 summary report...")
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("=" * 70 + "\n")
        f.write("Pangenome analysis summary report: gene family clustering and classification\n")
        f.write("=" * 70 + "\n\n")

        f.write("1. Basic information\n")
        f.write("-" * 40 + "\n")
        f.write(f"Analysis time: {pd.Timestamp.now()}\n")
        f.write(f"Number of strains: {len(parser.strain_names)}\n")
        f.write(f"OrthoFinder results directory: {ORTHOFINDER_RESULTS_DIR}\n")
        f.write(f"Output directory: {OUTPUT_DIR}\n\n")

        tree_file = parser.find_species_tree()
        if tree_file:
            f.write(f"Species tree file: {tree_file}\n\n")

        f.write("2. Gene family classification statistics\n")
        f.write("-" * 40 + "\n")
        if stats_df is not None:
            for _, row in stats_df.iterrows():
                f.write(f"{row['Classification']}: {row['Count']} ({row['Percentage']}%)\n")

            total_families = stats_df["Count"].sum()
            f.write(f"\nTotal: {total_families} gene families\n\n")

            core_total = stats_df[stats_df["Classification"].isin(["Core", "Soft Core"])]["Count"].sum()
            accessory_total = stats_df[stats_df["Classification"].isin(["Shell", "Cloud"])]["Count"].sum()
            private_total = stats_df[stats_df["Classification"] == "Private"]["Count"].sum() if "Private" in stats_df["Classification"].values else 0

            f.write("Simplified classification:\n")
            f.write(f"  Core genes (Core + Soft Core): {core_total} ({core_total / total_families * 100:.1f}%)\n")
            f.write(f"  Accessory genes (Shell + Cloud): {accessory_total} ({accessory_total / total_families * 100:.1f}%)\n")
            f.write(f"  Private genes: {private_total} ({private_total / total_families * 100:.1f}%)\n\n")
        else:
            f.write("No classification statistics were generated.\n\n")

        f.write("3. Gene family size statistics per strain\n")
        f.write("-" * 40 + "\n")
        if presence_df is not None:
            strain_counts = presence_df.sum(axis=0)
            f.write(f"Minimum: {strain_counts.min():.0f} gene families\n")
            f.write(f"Maximum: {strain_counts.max():.0f} gene families\n")
            f.write(f"Mean: {strain_counts.mean():.0f} gene families\n")
            f.write(f"Median: {strain_counts.median():.0f} gene families\n")
            f.write(f"Standard deviation: {strain_counts.std():.0f}\n")
            f.write(f"Q1: {strain_counts.quantile(0.25):.0f}\n")
            f.write(f"Q3: {strain_counts.quantile(0.75):.0f}\n\n")
        else:
            f.write("No gene family size statistics were generated.\n\n")

        f.write("4. Generated result files\n")
        f.write("-" * 40 + "\n")
        f.write(f"Results directory: {RESULTS_DIR}\n")
        for result_file in sorted(RESULTS_DIR.glob("*")):
            size_kb = os.path.getsize(result_file) / 1024 if os.path.exists(result_file) else 0
            f.write(f"  - {result_file.name} ({size_kb:.1f} KB)\n")

        f.write("\n5. Analysis parameters\n")
        f.write("-" * 40 + "\n")
        f.write(f"Core threshold: {CORE_THRESHOLD * 100}% of strains\n")
        f.write(f"Soft-core threshold: {SOFT_CORE_THRESHOLD * 100}% of strains\n")
        f.write(f"Shell threshold: {SHELL_THRESHOLD * 100}% of strains\n")
        f.write(f"Cloud threshold: <{CLOUD_THRESHOLD * 100}% of strains\n\n")

        f.write("=" * 70 + "\n")
        f.write("End of report\n")
        f.write("=" * 70 + "\n")
    logger.info(f"Part 1 summary report saved: {output_path}")


def generate_part2_summary_report(fit_results, percentage_df, output_path):
    logger.info("Generating part 2 summary report...")
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("=" * 70 + "\n")
        f.write("Pangenome analysis summary report: COG annotation and pangenome curve analysis\n")
        f.write("=" * 70 + "\n\n")

        f.write("1. Basic information\n")
        f.write("-" * 40 + "\n")
        f.write(f"Analysis time: {pd.Timestamp.now()}\n")
        f.write(f"Annotation directory: {ANNOTATION_DIR}\n")
        f.write(f"COG results directory: {COG_RESULTS_DIR}\n\n")

        f.write("2. COG annotation statistics\n")
        f.write("-" * 40 + "\n")
        cog_stats_file = COG_RESULTS_DIR / "cog_overall_statistics.csv"
        if cog_stats_file.exists():
            cog_stats = pd.read_csv(cog_stats_file)
            f.write(f"Total annotated genes: {cog_stats['Count'].sum():,}\n")
            f.write(f"Number of COG categories: {len(cog_stats)}\n\n")
            f.write("Top 10 COG categories:\n")
            for _, row in cog_stats.head(10).iterrows():
                desc = COG_CATEGORIES.get(row["COG_Category"], "Unknown")
                f.write(f"  {row['COG_Category']}: {row['Count']:,} ({row['Percentage']:.1f}%) - {desc}\n")
            f.write("\n")

        f.write("3. COG distribution by gene class\n")
        f.write("-" * 40 + "\n")
        if percentage_df is not None:
            core_top5 = percentage_df.nlargest(5, "Core")
            f.write("Top COG categories in core genes:\n")
            for _, row in core_top5.iterrows():
                desc = COG_CATEGORIES.get(row["COG_Category"], "Unknown")
                f.write(f"  {row['COG_Category']}: {row['Core']:.1f}% - {desc}\n")
            f.write("\n")

            accessory_top5 = percentage_df.nlargest(5, "Accessory")
            f.write("Top COG categories in accessory genes:\n")
            for _, row in accessory_top5.iterrows():
                desc = COG_CATEGORIES.get(row["COG_Category"], "Unknown")
                f.write(f"  {row['COG_Category']}: {row['Accessory']:.1f}% - {desc}\n")
            f.write("\n")

        f.write("4. Pangenome curve analysis\n")
        f.write("-" * 40 + "\n")
        if fit_results:
            f.write(f"Pangenome type: {fit_results['Pan_genome_openness']}\n")
            f.write(f"Heap's law gamma: {fit_results['Pan_genome_gamma']:.4f}\n")
            f.write(f"Heap's law R2: {fit_results['Pan_genome_R2']:.4f}\n")
            f.write(f"Final pan-genome size: {fit_results['Final_pan_size']:.0f} gene families\n")
            f.write(f"Final core-genome size: {fit_results['Final_core_size']:.0f} gene families\n")
            f.write(f"Core/pan-genome ratio: {(fit_results['Final_core_size'] / fit_results['Final_pan_size'] * 100):.1f}%\n\n")

        f.write("5. Generated result files\n")
        f.write("-" * 40 + "\n")
        f.write(f"COG results directory: {COG_RESULTS_DIR}\n")
        f.write("  - cog_overall_statistics.csv\n")
        f.write("  - cog_percentage_by_gene_class.csv\n")
        f.write("  - cog_counts_by_gene_class.csv\n")
        f.write("  - gene_families_with_cog.tsv\n\n")
        f.write(f"Pangenome curve files directory: {RESULTS_DIR}\n")
        f.write("  - pangenome_curve_data.csv\n")
        f.write("  - pangenome_fit_results.csv\n\n")

        f.write("=" * 70 + "\n")
        f.write("End of report\n")
        f.write("=" * 70 + "\n")
    logger.info(f"Part 2 summary report saved: {output_path}")


def run_part1():
    logger.info("=" * 60)
    logger.info("Starting part 1: gene family clustering and classification")
    logger.info("=" * 60)

    parser = OrthoFinderParser(ORTHOFINDER_RESULTS_DIR)

    orthogroups_path = parser.find_orthogroups_file()
    if not orthogroups_path:
        logger.error("Orthogroups.tsv file not found. Check the input path.")
        logger.error(f"Checked directory: {ORTHOFINDER_RESULTS_DIR}")
        return False, None, None, None

    orthogroups_df = parser.parse_orthogroups_file(orthogroups_path)
    if orthogroups_df is None:
        logger.error("Failed to parse Orthogroups file")
        return False, None, None, None

    presence_df, binary_df = parser.create_presence_absence_matrix(orthogroups_df)
    if presence_df is None or binary_df is None:
        logger.error("Failed to create presence/absence matrices")
        return False, parser, None, None

    classified_df, stats_df = parser.classify_gene_families(binary_df)
    if classified_df is None or stats_df is None:
        logger.error("Failed to classify gene families")
        return False, parser, presence_df, None

    report_path = OUTPUT_DIR / "analysis_summary_part1.txt"
    generate_part1_summary_report(parser, stats_df, presence_df, report_path)

    logger.info("=" * 60)
    logger.info("Part 1 completed")
    logger.info("=" * 60)

    print("\n" + "=" * 60)
    print("Part 1 completed")
    print("=" * 60)
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Total gene families: {len(classified_df)}")
    print(f"Number of strains: {len(parser.strain_names)}")

    print("\nGene family classification statistics:")
    for _, row in stats_df.iterrows():
        print(f"  {row['Classification']}: {row['Count']} ({row['Percentage']}%)")

    strain_counts = presence_df.sum(axis=0)
    print("\nGene family counts per strain:")
    print(f"  Minimum: {strain_counts.min():.0f}")
    print(f"  Maximum: {strain_counts.max():.0f}")
    print(f"  Mean: {strain_counts.mean():.0f}")
    print(f"  Median: {strain_counts.median():.0f}")

    return True, parser, presence_df, stats_df


def run_part2():
    logger.info("=" * 60)
    logger.info("Starting part 2: COG annotation and pangenome curve analysis")
    logger.info("=" * 60)

    required_files = [
        RESULTS_DIR / "gene_families_classified.tsv",
        RESULTS_DIR / "gene_family_binary_matrix.tsv"
    ]

    missing_files = [path for path in required_files if not path.exists()]
    if missing_files:
        logger.error(f"Missing part 1 result files: {missing_files}")
        logger.error("Run part 1 before part 2")
        return False

    logger.info("Step 1: COG annotation analysis")
    cog_analyzer = COGAnalyzer(ANNOTATION_DIR)

    if not cog_analyzer.parse_emapper_annotations():
        logger.error("Failed to parse COG annotations")
        return False

    gene_families_file = RESULTS_DIR / "gene_families_classified.tsv"
    gene_families_with_cog = cog_analyzer.map_cog_to_gene_families(gene_families_file)

    if gene_families_with_cog is None:
        logger.error("Failed to map COG annotations")
        return False

    percentage_df, count_df = cog_analyzer.analyze_cog_by_gene_class(gene_families_with_cog)

    if percentage_df is None or count_df is None:
        logger.error("Failed to analyze COG categories by gene class")
        return False

    logger.info("Step 2: Pangenome curve analysis")
    curve_analyzer = PangenomeCurveAnalyzer(RESULTS_DIR / "gene_family_binary_matrix.tsv")

    if not curve_analyzer.load_binary_matrix():
        logger.error("Failed to load binary matrix")
        return False

    curve_data = curve_analyzer.calculate_pangenome_curve(n_permutations=50, step_size=20)

    if curve_data is None:
        logger.error("Failed to calculate pangenome curves")
        return False

    fit_results = curve_analyzer.fit_pangenome_models(curve_data)

    report_path = BASE_DIR / "analysis_summary_part2.txt"
    generate_part2_summary_report(fit_results, percentage_df, report_path)

    logger.info("=" * 60)
    logger.info("Part 2 completed")
    logger.info("=" * 60)

    print("\n" + "=" * 60)
    print("Part 2 completed")
    print("=" * 60)

    cog_stats_file = COG_RESULTS_DIR / "cog_overall_statistics.csv"
    if cog_stats_file.exists():
        cog_stats = pd.read_csv(cog_stats_file)
        print("\nCOG annotation statistics:")
        print(f"  Total annotated genes: {cog_stats['Count'].sum():,}")
        print(f"  Number of COG categories: {len(cog_stats)}")
        print("  Most common COG categories:")
        for _, row in cog_stats.head(5).iterrows():
            desc = COG_CATEGORIES.get(row["COG_Category"], "Unknown")[:30]
            print(f"    {row['COG_Category']}: {row['Count']:,} ({row['Percentage']:.1f}%) - {desc}")

    if fit_results:
        print("\nPangenome analysis:")
        print(f"  Type: {fit_results['Pan_genome_openness']}")
        print(f"  Gamma: {fit_results['Pan_genome_gamma']:.4f}")
        print(f"  Final pan-genome size: {fit_results['Final_pan_size']:.0f}")
        print(f"  Final core-genome size: {fit_results['Final_core_size']:.0f}")
        print(f"  Core/pan-genome ratio: {(fit_results['Final_core_size'] / fit_results['Final_pan_size'] * 100):.1f}%")

    print(f"\nOutput directory: {BASE_DIR}")
    print(f"COG results directory: {COG_RESULTS_DIR}")
    print(f"Summary report: {report_path}")

    return True


def main():
    part1_success, _, _, _ = run_part1()
    if not part1_success:
        return False

    part2_success = run_part2()
    if not part2_success:
        return False

    return True


if __name__ == "__main__":
    success = main()
    if success:
        print("\nCombined analysis completed successfully.")
    else:
        print("\nCombined analysis failed. Check the log messages above.")
        sys.exit(1)
