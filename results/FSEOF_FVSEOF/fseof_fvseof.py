import os
import pandas as pd
from pathlib import Path
import glob

def process_strain_data():
    # Define file paths
    fvseof_base_path = "/data/gan_1/FVSEOF/model"
    fseof_base_path = "/data/gan_1/FSEOF/model"
    output_base_path = "/data/gan_1/aaaResults/FSEOF_FVSEOF/441models"
    
    # Create output directory
    os.makedirs(output_base_path, exist_ok=True)
    
    # Get all strain folders
    strain_folders = [f for f in os.listdir(fseof_base_path) 
                     if os.path.isdir(os.path.join(fseof_base_path, f))]
    
    print(f"Found {len(strain_folders)} strain folders")
    
    # Process each strain
    for strain in strain_folders:
        print(f"Processing strain: {strain}")
        
        # Create output directory for current strain
        strain_output_path = os.path.join(output_base_path, strain)
        os.makedirs(strain_output_path, exist_ok=True)
        
        # Get all organic acid files under FSEOF directory for current strain
        fseof_strain_path = os.path.join(fseof_base_path, strain)
        acid_files = glob.glob(os.path.join(fseof_strain_path, "AmplificationTargets_*.xlsx"))
        
        # Process each organic acid file
        for acid_file in acid_files:
            # Extract reaction name of organic acid
            acid_name = os.path.basename(acid_file).replace("AmplificationTargets_", "").replace(".xlsx", "")
            print(f"  Processing organic acid: {acid_name}")
            
            try:
                # Read FSEOF data - robust column handling
                fseof_df = pd.read_excel(acid_file, header=None)
                
                # Check column number and keep only the first two columns
                if fseof_df.shape[1] >= 2:
                    fseof_df = fseof_df.iloc[:, :2]
                    fseof_df.columns = ['reaction', 'q_slope']
                else:
                    print(f"    Warning: {acid_file} has less than 2 columns, skipped")
                    continue
                
                # Ensure q_slope is numeric type
                fseof_df['q_slope'] = pd.to_numeric(fseof_df['q_slope'], errors='coerce')
                
                # Separate up-regulated and down-regulated targets
                fseof_up = set(fseof_df[fseof_df['q_slope'] > 0]['reaction'].dropna().tolist())
                fseof_down = set(fseof_df[fseof_df['q_slope'] < 0]['reaction'].dropna().tolist())
                
                # Build file paths for FVSEOF results
                fvseof_up_file = os.path.join(fvseof_base_path, strain, f"fvseof_up_targets_{acid_name}.csv")
                fvseof_down_file = os.path.join(fvseof_base_path, strain, f"fvseof_down_targets_{acid_name}.csv")
                
                # Read FVSEOF up-regulated targets
                fvseof_up = set()
                if os.path.exists(fvseof_up_file):
                    try:
                        fvseof_up_df = pd.read_csv(fvseof_up_file, header=None)
                        if fvseof_up_df.shape[1] >= 1:
                            fvseof_up = set(fvseof_up_df.iloc[:, 0].dropna().tolist())
                    except Exception as e:
                        print(f"    Failed to read FVSEOF up-regulated file: {str(e)}")
                
                # Read FVSEOF down-regulated targets
                fvseof_down = set()
                if os.path.exists(fvseof_down_file):
                    try:
                        fvseof_down_df = pd.read_csv(fvseof_down_file, header=None)
                        if fvseof_down_df.shape[1] >= 1:
                            fvseof_down = set(fvseof_down_df.iloc[:, 0].dropna().tolist())
                    except Exception as e:
                        print(f"    Failed to read FVSEOF down-regulated file: {str(e)}")
                
                # Calculate intersection of target sets
                common_up = fseof_up.intersection(fvseof_up)
                common_down = fseof_down.intersection(fvseof_down)
                
                # Define output file path
                output_file = os.path.join(strain_output_path, f"{acid_name}.xlsx")
                
                with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                    # Write up-regulated targets sheet
                    if common_up:
                        up_data = []
                        for reaction in common_up:
                            # Query corresponding q_slope value
                            match = fseof_df[fseof_df['reaction'] == reaction]
                            if not match.empty:
                                q_slope = match['q_slope'].values[0]
                                up_data.append({'reaction': reaction, 'q_slope': q_slope})
                        if up_data:
                            up_df = pd.DataFrame(up_data)
                            up_df.to_excel(writer, sheet_name='up', index=False)
                        else:
                            pd.DataFrame(columns=['reaction', 'q_slope']).to_excel(writer, sheet_name='up', index=False)
                    else:
                        # Create empty sheet if no common up-regulated targets
                        pd.DataFrame(columns=['reaction', 'q_slope']).to_excel(writer, sheet_name='up', index=False)
                    
                    # Write down-regulated targets sheet
                    if common_down:
                        down_data = []
                        for reaction in common_down:
                            # Query corresponding q_slope value
                            match = fseof_df[fseof_df['reaction'] == reaction]
                            if not match.empty:
                                q_slope = match['q_slope'].values[0]
                                down_data.append({'reaction': reaction, 'q_slope': q_slope})
                        if down_data:
                            down_df = pd.DataFrame(down_data)
                            down_df.to_excel(writer, sheet_name='down', index=False)
                        else:
                            pd.DataFrame(columns=['reaction', 'q_slope']).to_excel(writer, sheet_name='down', index=False)
                    else:
                        # Create empty sheet if no common down-regulated targets
                        pd.DataFrame(columns=['reaction', 'q_slope']).to_excel(writer, sheet_name='down', index=False)
                
                print(f"    Generated: {acid_name}.xlsx (Up: {len(common_up)}, Down: {len(common_down)})")
                
            except Exception as e:
                print(f"    Error occurred while processing {acid_name}: {str(e)}")
                # Print debug information
                print(f"    File path: {acid_file}")
                if 'fseof_df' in locals():
                    print(f"    FSEOF data shape: {fseof_df.shape}")
                continue
    
    print("All strains processed!")

def main():
    try:
        process_strain_data()
    except KeyboardInterrupt:
        print("\nProgram interrupted by user")
    except Exception as e:
        print(f"Program execution error: {str(e)}")

if __name__ == "__main__":
    main()
