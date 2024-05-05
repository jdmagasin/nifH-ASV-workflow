#!/usr/bin/env python3
"""
query_cmap.py 

This module queries the CMAP data repository to colocalize environmental metadata with user provided sample location metadata.

Overview:
This script colocalizes sample collection metadata with the Simon's Collaborative Marine Atlas Project (CMAP) repository to retrieve environmental data of interest to enhance downstream analysis. Input collection data, as a CSV, is supplied along with targets and parameters which are used to query the CMAP catalog, returning a colocalized dataset as a compressed CSV file. 

Collection data **MUST** be supplied as a CSV that is formatted specifically for the CMAP query to execute properly. The required columns headers and their expected data types for the sample collection data:
      "header":   data type
      - "lat":    float
      - "lon":    float
      - "time":   str
      - "depth":  float
To avoid data typing issues, the script will first attempt to convert values to the proper data type prior to attempting to colocalize the data. If this conversation fails, the entire row (sample) is removed. Please make sure your data if formatted properly prior to running this script. The script ignores all other column headers so these restrictions do not apply to them. These column can be used for other relevant metadata, e.g., sample IDs.

The targets and their parameters have been set and reflect those used to generate the nifH ASV database CMAP data file. However, the tables, targets, and parameters can all be adjusted based at the users discretion by editing the targets dictionary in the setup_colocalization() function below. See the CMAP webpage and catalog for more information. https://simonscmap.com

A ***valid API key is required** to query the CMAP repository. The script will fail if a valid API key is not passed. Your API key can be obtained here (https://simonscmap.com/apikeymanagement). If the accompanied Snakefile is used to manage the this script within the CMAP stage of the nifH workflow, your API key needs to be added to the associated config file. If using this script directly, see Usage below. 

#_# Specific to the nifH-ASV-workflow #_#
The CMAP stage is part of the nifH-ASV-workflow, a customizable bioinformatic workflow that compiles multiple datasets and produces high quality ASVs  through a series of comprehensive stages. The CMAP stage is managed by a Snakefile that creates, activates, and deactivates a separate, self-contained conda environment upon the completion of the stage, having no effect on the conda environments required for the DADA2 and post-pipeline workflows. The main output from the GatherMetadata stage, metadata.cmap.csv, a table of the collection coordinates, dates at local noon, and depths from all the samples, is passed to this script (query_cmap.py) to query the CMAP data portal for co-localized environmental data.  
#_# Specific to the nifH-ASV-workflow #_#

Usage:
python3 query_cmap.py <input_path> <output_path> <your_cmap_api_key> -h,--help

Arguments:
    input: Path to the input dataset CSV file.
    output: Path to save the compressed colocalized dataset CSV file.
    cmap_api_key: **Your CMAP API key. Get it at https://simonscmap.com/apikeymanagement.**
    -h, --help: Show this help message and exit
"""

import os
import sys
import argparse
import gzip
import shutil
import subprocess
import pycmap
import pandas as pd
from argparse import RawTextHelpFormatter


__author__ = "Michael (Mo) Morando"
__copyright__ = "Copyright 2023"
__maintainer__ = "Michael (Mo) Morando"
__email__ = "mikemo@gmail.com"
__status__ = "Stable"


# Get the full path to the currently running script
script_path: str = __file__

# Extract the name of the script for use later
script_name: str = os.path.basename(p=script_path)

# General print statement to show script was called
print(f"{script_name} is currently running...")


# _# Script starts #_#
def read_input_file(input_path: str) -> pd.DataFrame | None:
    """
    Read the input dataset from a CSV file and return it as a DataFrame.

    Args:
        input_path (str): Path to the input dataset CSV file.

    Returns:
        pd.DataFrame: DataFrame containing the input data.
    """
    try:
        print("Attempting to read input file...")
        df: pd.DataFrame = pd.read_csv(
            filepath_or_buffer=input_path, sep=",", comment="#"
        )
        # -# check if pd.DataFrame is 'empty'
        #!# however, it appears that if file is completely empty, pd.read_csv()
        #!# does not read it in so it is not caught by this if statement but
        #!# instead dealt with separately
        if not df.empty:
            required_headers: list[str] = ["lat", "lon", "time", "depth"]
            missing_headers: list[str] = [
                head for head in required_headers if head not in df.columns
            ]
            # Check if there is a problem with the headers
            if missing_headers:
                missing_headers_str: str = ", ".join(missing_headers)
                raise ValueError(
                    f"Missing required columns: {missing_headers_str}.\nExpected columns are: {', '.join(required_headers)}\nPlease check your column headers on file: '{os.path.basename(p=input_path)}'"
                )
            # -# Check and remove entire rows with incorrect data types
            # supply required column headers and their data types
            # ONLY these will be checked, all other headers will be ignored
            #!# Major issues can occur downstream when executing query if the
            #!# expected values are not correct, e.g,. supplying 'surface' as a
            #!# depth. To handle this, the entire row is removed if the data
            #!# type of a single value is not what is expected
            expected_types: dict[str, type[float] | type[str]] = {
                "lat": float,
                "lon": float,
                "time": str,
                "depth": float,
            }
            # Print out the expected data types and the start of the validation process
            print(f"Expecting {expected_types}")
            print("Validating the data type of the values of these columns...")
            # Initialize lists to store:
            # valid rows
            # count of removed rows
            valid_rows = []
            rows_removed_count = 0
            # Loop through each row in the DataFrame
            for index, row in df.iterrows():
                # Set row to valid
                valid_row = True
                # Loop through each column and its expected data type
                for column, data_type in expected_types.items():
                    try:
                        # Attempt to convert the value in the row to the
                        # expected data type
                        row[column] = data_type(row[column])
                    # If conversion fails, print ValueError message then..
                    except ValueError:
                        print(
                            f"Row {index}: Value '{row[column]}' in column '{column}' is not convertible to {data_type.__name__}, removing row..."
                        )
                        # Mark row as invalid
                        valid_row = False
                        # Break out of inner loop because we don't want this
                        # row anymore
                        break
                # If the row is valid after all columns are evaluated, add it to the list of valid rows
                if valid_row:
                    valid_rows.append(row)
                # If the row is invalid, increment the count of removed rows
                # variable
                else:
                    rows_removed_count += 1
            # Print to indicate how many rows were removed
            print(
                f"Total number of rows REMOVED that could not be converted: {rows_removed_count}"
            )
            # Create a new DataFrame containing only the valid rows
            # the index needs to be reset or downstream concatenation fail
            # If the conversion process removes all the rows leaving the df empty, we can check for that here:
            if not len(valid_rows) <= 1:
                df = pd.DataFrame(data=valid_rows).reset_index(drop=True)
            else:
                print(
                    f"\nError: The input data file '{input_path}' values are not the expected types. The data type conversion checking process removed all rows from the DataFrame. Please see in above in log to determine which rows and values need to be adjusted!"
                )
                print("Contents of converted input file:")
                print(valid_rows)
                sys.exit(10)
            # IF everything checks out, return the usable df
            return df
        # If empty...
        else:
            print(f"Input data file {input_path} headers are there but rows are empty!")
            print("Contents of input file:")
            print(df)
            sys.exit(10)

    except FileNotFoundError:
        print(f"Input file '{input_path}' not found")
        sys.exit(1)
    except pd.errors.EmptyDataError as e:
        print(f"Error: {e}")
        print(f"Input data file {input_path} is empty!")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def contains_letters_and_numbers(s: str) -> bool:
    """
    Check if a string contains both letters and numbers.

    Args:
        s (str): Input string to check.

    Returns:
        bool: True if the string contains both letters and numbers, False otherwise.
    """
    try:
        if not s:
            raise ValueError("Input string cannot be empty")

        return any(char.isalpha() for char in s) and any(char.isdigit() for char in s)

    except TypeError as e:
        print(f"Error: {e}")
        return False
    except ValueError as e:
        print(f"Error: {e}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False


def setup_colocalization(cmap_api_key: str) -> dict:
    """
    Setup colocalization parameters.

    Args:
        cmap_api_key (str): CMAP API key.

    Returns:
        dict: Dictionary specifying the targets, tables, variables, and tolerances.
    """

    # Check if cmap_api_key contains both alphanumeric characters to help
    # validate it is an actual API key
    # Passing a single number or letter does not trip pycmap's error handling
    # for invalid API keys so this helps with that
    if not contains_letters_and_numbers(cmap_api_key):
        raise ValueError(
            "cmap_api_key must contain both alphabetic and numeric characters"
        )

    # Set the CMAP API key
    pycmap.API(token=cmap_api_key)

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
        #         # wind_speed", #! no longer available in catalog
        #         # "wind_curl",  #! FIXME: Exists in CMAP docs of Feb 2021 but causes KeyError crash
        #         "stress_curl", #! Exists in CMAP docs of Feb 2021 but causes KeyError crash
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
    Colocalize input data with CMAP repository.
    Required Columns: [lat, lon, time, depth]

    Args:
        input_data (pd.DataFrame): Input dataset.
        targets (dict): Dictionary specifying the targets, tables, variables, and tolerances.

    Returns:
        pd.DataFrame: Colocalized dataset.
    """
    # Print message indicating colocalization process initiation
    print(
        f"""\nColocalizing data with CMAP.\nThis should take a while depending on the size of your dataframe and number of cores.\n"""
    )

    # Colocalize df with CMAP repository
    colocalized_data: pd.DataFrame = pycmap.Sample(
        source=input_data, targets=targets, replaceWithMonthlyClimatolog=True
    )

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
    Main function to execute the script when called from the command line. This function handles the following steps:
    1. Parsing command-line arguments using argparse.
    2. Setting up the colocalization parameters.
    3. Colocalizing the input data with the CMAP database.
    4. Compressing and saving the colocalized dataset.

    Returns:
        None
    """
    # Parse command-line arguments

    # Create an ArgumentParser object to handle command-line arguments
    parser = argparse.ArgumentParser(
        description=(
            f"""
{script_name} colocalizes sample collection data with the Simon's Collaborative Marine Atlas Project (CMAP) repository to retrieve environmental data of interest to be used in analysis. Input collection data, as a CSV, is supplied along with targets and parameters which are used to query the CMAP catalog, returning a colocalized dataset as a compressed CSV file.     
      
Collection data **MUST** be supplied as a CSV that is formatted specifically for the CMAP query to execute properly. The required columns headers and their expected data types for the sample collection data:
      "header":   data type
      - "lat":    float
      - "lon":    float
      - "time":   str
      - "depth":  float
The script will attempt to convert values to the proper data type prior to querying CMAP. If this conversation fails, the entire row (sample) is removed. Please make sure your data if formatted properly prior to running this script. All other column headers can be any value and in any format as the script ignores them. These column can be used for other relevant metadata, e.g.,  sample IDs.
      
The targets and their parameters have been set and reflect those used to generate the nifH ASV database CMAP data file. However, the tables, targets, and parameters can all be adjusted based at the users discretion by editing the targets dictionary in the setup_colocalization() function below. See the CMAP webpage and catalog for more information. https://simonscmap.com
      
A ***valid API key is required** to query the CMAP repository. The script will fail if a valid API key is not passed. Your API key can be obtained here (https://simonscmap.com/apikeymanagement). If using this script directly, see Usage below. If the accompanied Snakefile is used to manage the this script within the CMAP stage of the nifH workflow, your API key needs to be added to the associated config file. 

"""
        ),
        formatter_class=RawTextHelpFormatter,
        usage=f"{script_name} [-h,--help] <input_path> <output_path> <your_cmap_api_key>",
    )
    # Add arguments for input, output, and CMAP API key
    parser.add_argument(
        "input",
        type=str,
        help="Path to the input dataset CSV file.",
    )
    parser.add_argument(
        "output",
        type=str,
        help="Path to save the colocalized dataset CSV file.",
    )
    parser.add_argument(
        "cmap_api_key",
        type=str,
        help="""
        Your specific CMAP API key. This can be obtained at:
        https://simonscmap.com/apikeymanagement
        """,
    )
    # Parse the provided arguments
    args, unknown_args = parser.parse_known_args()

    try:
        # if len(sys.argv) != 4:
        #   raise argparse.ArgumentError(
        #       argument=None,
        #       message=f"Total number of args is not correct. Expected 3 arguments, got {len(sys.argv) - 1}!\nArguments passed were:\n{'\n'.join((sys.argv)[1:])}",
        #   )
        # Check for any unknown arguments and raise an error if found
        if unknown_args:
            raise argparse.ArgumentError(
                argument=None,
                message=f"Unrecognized arguments: {' '.join(unknown_args)}",
            )

        # Check if all required arguments are provided
        if args.input is None and args.output is None and args.cmap_api_key is None:
            parser.print_help()
            exit()
        # return

    except argparse.ArgumentError as e:
        print(f"Error: {e}")
        print("Failed to parse command-line arguments.")
        sys.exit(5)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        print(f"Failed to execute script: '{script_name}'.")

    # Pass success flag
    success = False

    try:

        # Check if help message is requested
        # Exit gracefully without setting success flag

        # Setup colocalization parameters
        targets = setup_colocalization(cmap_api_key=args.cmap_api_key)

        # Print Pycmap version
        print(f"pycmap version: {pycmap.__version__}\n")

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
        print("Failed to parse * command-line arguments.")
    except FileNotFoundError:
        print(f"Input file '{args.input}' not found")
        # raise FileNotFoundError(f"Input file '{args.input}' not found")
    # except pd.errors.EmptyDataError as e:
    #     print(f"Error: {e}")
    #     print("Failed to colocalize data: Input data file is empty.")
    except KeyError as e:
        print(
            f"""Error: Column '{e.args[0]}' not found in the input DataFrame.
If error states "'None' is missing", an argument was not pass. Please see usage under --help and check your values.
"""
        )
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        print("Failed to colocalize data and save the dataset.")
    finally:
        # success block with print statements
        if success:
            print(f"\nScript '{script_name}' executed successfully!!")
        else:
            print(
                f"""\nScript '{script_name}' exited with an error.
The CMAP repository was not queried and no output was not written!
"""
            )


# Conditional block ensures that the `main()` function is executed only when
# called from the command line
if __name__ == "__main__":
    main()
