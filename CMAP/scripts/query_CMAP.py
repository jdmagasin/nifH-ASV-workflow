import os
import argparse
import gzip
import shutil
import subprocess
import pycmap
import pandas as pd

"""
Script Name: CMAP Data Colocalization

Description:
This script colocalizes input data with the CMAP (Simons Collaboration on Ocean Processes and Ecology) database. It takes an input dataset in CSV format, 
colocalizes the data with specified targets in the CMAP database, and saves the 
colocalized dataset as a compressed CSV file.

Usage:
python cmap_colocalization.py input <input_path> output <output_path> --cmap_api_key <your_cmap_api_key>

Arguments:
    input: Path to the input dataset CSV file.
    output: Path to save the colocalized dataset CSV file.
    --cmap_api_key: Your CMAP API key. Get it at https://simonscmap.com/apikeymanagement.
"""

# Get the full path to the currently running script
script_path: str = __file__

# Extract the name of the script for use later
script_name: str = os.path.basename(p=script_path)

print(f"{script_name} is currently running...")


def read_input_file(input_path: str) -> pd.DataFrame:
    """
    Read the input file and return a DataFrame.

    Args:
        input_path (str): Path to the input dataset CSV file.

    Returns:
        pd.DataFrame: DataFrame containing the input data.
    """
    return pd.read_csv(filepath_or_buffer=input_path, sep=",", comment="#")


def setup_colocalization(args: argparse.Namespace) -> dict:
    """
    Setup colocalization parameters.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.

    Returns:
        dict: Dictionary specifying the targets, tables, variables, and tolerances.
    """
    # Set the CMAP API key
    pycmap.API(token=args.cmap_api_key)

    # Define targets
    targets: dict[str, dict[str, list[str] | list[int | float]]] = {
        "tblAltimetry_REP_Signal": {
            "variables": ["sla"],
            "tolerances": [
                1,
                0.25,
                0.25,
                1,
            ],
        },
        "tblAltimetry_NRT_Signal": {
            "variables": ["sla"],
            "tolerances": [
                1,
                0.25,
                0.25,
                1,
            ],
        },
        "tblCHL_REP": {"variables": ["chl"], "tolerances": [4, 0.25, 0.25, 1]},
        "tblModis_AOD_REP": {
            "variables": ["AOD"],
            "tolerances": [
                15,
                0.5,
                0.5,
                1,
            ],  # Sources of these particles include: volcanic ash, wildfire smoke, windblown sand and dust
        },
        "tblModis_PAR": {"variables": ["PAR"], "tolerances": [15, 0.5, 0.5, 1]},
        "tblSSS_NRT": {"variables": ["sss"], "tolerances": [1, 0.25, 0.25, 1]},
        "tblSST_AVHRR_OI_NRT": {
            "variables": ["sst"],
            "tolerances": [1, 0.25, 0.25, 1],
        },
        # "tblWind_NRT_hourly": {
        #     "variables": [
        #         # wind_speed", # no longer available
        #         # "wind_curl",  # FIXME: Exists in CMAP docs of Feb 2021 but causes KeyError crash
        #         "stress_curl",
        #     ],
        #     "tolerances": [1, 0.25, 0.25, 1],
        # },
        "tblDarwin_Nutrient": {
            "variables": ["DIN", "PO4", "FeT", "O2", "SiO2"],
            "tolerances": [2, 0.5, 0.5, 5],
        },
        "tblDarwin_Ecosystem": {
            "variables": [
                "phytoplankton_diversity_shannon_index",
                "phytoplankton",
                "zooplankton",
                "CHL",
                "primary_production",
            ],
            "tolerances": [2, 0.5, 0.5, 5],
        },
        "tblDarwin_Phytoplankton": {
            "variables": [
                "diatom",
                "coccolithophore",
                "mixotrophic_dinoflagellate",
                "picoeukaryote",
                "picoprokaryote",
            ],
            "tolerances": [2, 0.5, 0.5, 5],
        },
        "tblDarwin_Nutrient_Climatology": {
            "variables": [
                "DIC_darwin_clim",
                "NH4_darwin_clim",
                "NO2_darwin_clim",
                "NO3_darwin_clim",
                "PO4_darwin_clim",
                "SiO2_darwin_clim",
                "FeT_darwin_clim",
                "DOC_darwin_clim",
                "DON_darwin_clim",
                "DOP_darwin_clim",
                "DOFe_darwin_clim",
                "POC_darwin_clim",
                "PON_darwin_clim",
                "POSi_darwin_clim",
                "POFe_darwin_clim",
                "PIC_darwin_clim",
                "ALK_darwin_clim",
                "O2_darwin_clim",
                "CDOM_darwin_clim",
            ],
            "tolerances": [2, 0.5, 0.5, 5],
        },
        "tblPisces_NRT": {
            "variables": ["NO3", "PO4", "Fe", "O2", "Si", "PP", "PHYC", "CHL"],
            "tolerances": [4, 0.5, 0.5, 5],
        },
        # 	"tblArgoMerge_REP": {
        #                                  "variables": ["argo_merge_salinity_adj", "argo_merge_temperature_adj", "argo_merge_O2_adj", "argo_merge_turbidity_adj", "argo_merge_chl_adj", "argo_merge_cdom_adj", "argo_merge_NO3", "argo_merge_NO3_adj",
        # 		"argo_merge_ph_adj"],
        #                                  "tolerances": [1, 1, 1, 5]
        #                                  },
        "tblArgo_MLD_Climatology": {
            "variables": [
                "mls_da_argo_clim",
                "mls_dt_argo_clim",
                "mlt_da_argo_clim",
                "mlt_dt_argo_clim",
                "mlpd_da_argo_clim",
                "mlpd_dt_argo_clim",
                "mld_da_mean_argo_clim",
                "mld_dt_mean_argo_clim",
                "mld_da_median_argo_clim",
                "mld_dt_median_argo_clim",
            ],
            "tolerances": [1, 1, 1, 5],
        },
        "tblGlobalDrifterProgram": {
            "variables": ["sst"],
            "tolerances": [1, 0.25, 0.25, 1],
        },
        # 	"tblMercator_MLD_NRT": {
        #                                  "variables": ["mld_nrt"],
        #                                  "tolerances": [1, 1, 1, 5]
        #                                  },
        "tblWOA_2018_MLD_qrtdeg_Climatology": {
            "variables": ["M_an_clim", "M_mn_clim"],
            "tolerances": [1, 0.5, 0.5, 5],
        },
        "tblWOA_2018_MLD_1deg_Climatology": {
            "variables": ["M_an_clim", "M_mn_clim"],
            "tolerances": [1, 0.5, 0.5, 5],
        },
        "tblWOA_2018_qrtdeg_Climatology": {
            "variables": [
                "C_an_clim",
                "C_mn_clim",
                "s_an_clim",
                "s_mn_clim",
                "t_an_clim",
                "t_mn_clim",
                "I_an_clim",
                "I_mn_clim",
            ],
            "tolerances": [1, 0.5, 0.5, 5],
        },
        "tblWOA_2018_1deg_Climatology": {
            "variables": [
                "C_an_clim",
                "C_mn_clim",
                "s_an_clim",
                "s_mn_clim",
                "t_an_clim",
                "t_mn_clim",
                "A_an_clim",
                "A_mn_clim",
                "O_an_clim",
                "O_mn_clim",
                "I_an_clim",
                "I_mn_clim",
                "n_an_clim",
                "n_mn_clim",
                "p_an_clim",
                "p_mn_clim",
                "si_an_clim",
                "si_mn_clim",
            ],
            "tolerances": [1, 1, 1, 5],
        },
        "tblWOA_Climatology": {
            "variables": [
                "sea_water_temp_WOA_clim",
                "density_WOA_clim",
                "salinity_WOA_clim",
                "nitrate_WOA_clim",
                "phosphate_WOA_clim",
                "silicate_WOA_clim",
                "oxygen_WOA_clim",
                "AOU_WOA_clim",
                "o2sat_WOA_clim",
                "conductivity_WOA_clim",
            ],
            "tolerances": [1, 1, 1, 5],
        },
    }

    return targets


def colocalize_data(input_data: pd.DataFrame, targets: dict) -> pd.DataFrame:
    """
    Colocalize input data with CMAP.

    Args:
        input_data (pd.DataFrame): Input dataset.
        targets (dict): Dictionary specifying the targets, tables, variables, and tolerances.

    Returns:
        pd.DataFrame: Colocalized dataset.
    """
    # Print message indicating colocalization process initiation
    print(
        f"""\nThis should take a while depending on the size of your dataframe and number of cores.\nColocalizing data with CMAP.\n"""
    )

    colocalized_dfs: list[Any] = []
    for target_name, target_info in targets.items():
        colocalized_df: Any = pycmap.Sample(
            source=input_data,
            targets={target_name: target_info},
            replaceWithMonthlyClimatolog=True,
        )
        colocalized_dfs.append(colocalized_df)

    # Concatenate colocalized DataFrames
    colocalized_data: pd.DataFrame = pd.concat(objs=colocalized_dfs, axis=1)

    return colocalized_data


def compress_and_save(df: pd.DataFrame, output_path: str) -> None:
    """
    Compresses the DataFrame to a gzip file and saves it.

    Args:
        df (pd.DataFrame): DataFrame to be compressed and saved.
        output_path (str): Path to save the compressed file.
    """
    # Save the colocalized dataset as a temporary CSV file
    temp_csv_path = "temp_colocalized_df.csv"
    df.to_csv(path_or_buf=temp_csv_path, index=False)

    # Compress the temporary CSV file with gzip
    with open(file=temp_csv_path, mode="rb") as f_in, gzip.open(
        filename=output_path, mode="wb"
    ) as f_out:
        shutil.copyfileobj(fsrc=f_in, fdst=f_out)

    print(f"\nColocalized dataset saved and compressed to {output_path}")

    # Remove the temporary CSV file
    subprocess.run(args=["rm", temp_csv_path])

    print("Temporary CSV file removed.\n")


def main() -> None:
    """
    Main function to execute the script. This function handles the following steps:
    1. Parsing command-line arguments using argparse.
    2. Setting up the colocalization parameters.
    3. Colocalizing the input data with the CMAP database.
    4. Compressing and saving the colocalized dataset.

    Returns:
        None
    """
    try:
        # Parse command-line arguments

        # Create an ArgumentParser object to handle command-line arguments
        parser = argparse.ArgumentParser(description="Process input and output paths.")
        # Add arguments for input, output, and CMAP API key
        parser.add_argument(
            "input", type=str, help="Path to the input dataset CSV file."
        )
        parser.add_argument(
            "output", type=str, help="Path to save the colocalized dataset CSV file."
        )
        parser.add_argument(
            "--cmap_api_key",
            type=str,
            help="""
            Your specific CMAP API key. This can be obtained at https://simonscmap.com/apikeymanagement
            """,
        )
        # Parse the provided arguments
        args, unknown_args = parser.parse_known_args()

        # Check for any unknown arguments and raise an error if found
        if unknown_args:
            raise argparse.ArgumentError(
                argument=None,
                message=f"UnreTESTcognized arguments: {' '.join(unknown_args)}",
            )

        # Check if all required arguments are provided
        if args.input is None and args.output is None and args.cmap_api_key is None:
            parser.print_help()
            exit()
        # return

    except argparse.ArgumentError as e:
        print(f"Error: {e}")
        print("Failed to parse command-line arguments.")
    except Exception as e:
        print(f"Error: {e}")
        print(f"Failed to execute script: '{script_name}'.")

    # Pass success flag
    success = False

    try:

        # Check if help message is requested
        # Exit gracefully without setting success flag

        # Setup colocalization parameters
        targets = setup_colocalization(args=args)

        # Print Pycmap version
        print(f"pycmap version: {pycmap.__version__}")

        # Load data from the specified input path
        input_data: DataFrame = read_input_file(input_path=args.input)

        # Colocalize data
        colocalized_data: DataFrame = colocalize_data(
            input_data=input_data, targets=targets
        )

        # Compress and save the colocalized dataset
        compress_and_save(df=colocalized_data, output_path=args.output)

        # Change success flag to true with successful completion
        success = True

    # Exceptions and Finally
    except UnboundLocalError as e:
        print(f"Error: {e}")
        print("Failed to parse command-line arguments.")
    except FileNotFoundError:
        print(f"Input file '{args.input}' not found")
        # raise FileNotFoundError(f"Input file '{args.input}' not found")
    except pd.errors.EmptyDataError as e:
        print(f"Error: {e}")
        print("Failed to colocalize data: Input data file is empty.")
    except KeyError as e:
        print(
            f"""Error: Column '{e.args[0]}' not found in the input DataFrame.
If error states "'None' is missing", an argument was not pass. Please see usage under --help and check your values.
"""
        )
    except Exception as e:
        print(f"Error: {e}")
        print("Failed to colocalize data and save the dataset.")
    finally:
        # success block with print statements
        if success:
            print(f"Script '{script_name}' executed successfully!!")
        else:
            print(
                f"""Script '{script_name}' exited with an error. No consensus ID
was made and no output was not written! See logs for details...
"""
            )


# Conditional block ensures that the `main()` function is executed only when
# called from the command line
if __name__ == "__main__":
    main()
