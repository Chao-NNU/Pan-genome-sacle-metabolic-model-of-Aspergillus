import os
import pandas as pd
from pathlib import Path
import glob

def process_strain_data():
    fseof_base_path = "/data/gan_1/FSEOF/model"
    fvseof_base_path = "/data/gan_1/FVSEOF/model"
    output_base_path = "/data/gan_1/addModels/FSEOF_FVSEOF"

    os.makedirs(output_base_path, exist_ok=True)

    strain_folders = [f for f in os.listdir(fseof_base_path) 
                     if os.path.isdir(os.path.join(fseof_base_path, f))]

    print(f"Found {len(strain_folders)} strain folders")

    for strain in strain_folders:
        print(f"Processing strain: {strain}")

        strain_output_path = os.path.join(output_base_path, strain)
        os.makedirs(strain_output_path, exist_ok=True)

        fseof_strain_path = os.path.join(fseof_base_path, strain)
        fvseof_strain_path = os.path.join(fvseof_base_path, strain)

        if not os.path.exists(fseof_strain_path):
            print(f"Warning: FSEOF path does not exist: {fseof_strain_path}")
            continue

        if not os.path.exists(fvseof_strain_path):
            print(f"Warning: FVSEOF path does not exist: {fvseof_strain_path}")
            continue

        fseof_files = glob.glob(os.path.join(fseof_strain_path, "AmplificationTargets_*.xlsx"))

        for fseof_file in fseof_files:
            filename = os.path.basename(fseof_file)
            acid_reaction = filename.replace("AmplificationTargets_", "").replace(".xlsx", "")

            print(f"  Processing organic acid: {acid_reaction}")

            fvseof_up_file = os.path.join(fvseof_strain_path, f"fvseof_up_targets_{acid_reaction}.csv")
            fvseof_down_file = os.path.join(fvseof_strain_path, f"fvseof_down_targets_{acid_reaction}.csv")

            if not os.path.exists(fvseof_up_file) or not os.path.exists(fvseof_down_file):
                print(f"    Warning: FVSEOF files do not exist for {acid_reaction}")
                continue

            try:
                fseof_df = pd.read_excel(fseof_file, header=0)
                if fseof_df.shape[1] >= 2:
                    fseof_df.columns = ['reaction', 'q_slope'] + list(fseof_df.columns[2:])
                else:
                    print(f"    FSEOF file has insufficient columns: {fseof_file}")
                    continue

                fseof_up_targets = set(fseof_df[fseof_df['q_slope'] > 0]['reaction'].tolist())
                fseof_down_targets = set(fseof_df[fseof_df['q_slope'] < 0]['reaction'].tolist())

                fvseof_up_df = pd.read_csv(fvseof_up_file, header=None)
                fvseof_up_targets = set(fvseof_up_df[0].tolist())

                fvseof_down_df = pd.read_csv(fvseof_down_file, header=None)
                fvseof_down_targets = set(fvseof_down_df[0].tolist())

                common_up_targets = sorted(list(fseof_up_targets & fvseof_up_targets))
                common_down_targets = sorted(list(fseof_down_targets & fvseof_down_targets))

                output_file = os.path.join(strain_output_path, f"{acid_reaction}.xlsx")

                with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                    up_df = pd.DataFrame(common_up_targets, columns=['reaction'])
                    up_df.to_excel(writer, sheet_name='up', index=False)

                    down_df = pd.DataFrame(common_down_targets, columns=['reaction'])
                    down_df.to_excel(writer, sheet_name='down', index=False)

                print(f"    Completed: {acid_reaction} - up intersection: {len(common_up_targets)}, down intersection: {len(common_down_targets)}")

            except Exception as e:
                print(f"    Error processing {acid_reaction}: {str(e)}")
                continue

    print("All strains processed")

def main():
    try:
        process_strain_data()
    except KeyboardInterrupt:
        print("\nProgram interrupted by user")
    except Exception as e:
        print(f"Program execution error: {str(e)}")

if __name__ == "__main__":
    main()
