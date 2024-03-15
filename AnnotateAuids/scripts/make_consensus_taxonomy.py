"""
make_consensus_taxonomy.py

This script is part of the AnnotateAuids stage of the nifH ASV workflow. It is managed by a Makefile and depends on the merged annotation table (auids.annot.tsv) created by this stage. 

Overview:
    - reads an input annotation table tsv file
    - cleans up columns in the DataFrame
    - creates a consensus taxonomy ID -> consensus_id
    - writes the updated annotation table with new column consensus_id added
A new column "consensus_id" is added to the input table (auids.annot.tsv), representing the most informative taxonomic ID aggregated across the merged annotation table supplied and new thresholds passed to this script. This ID should be used in any discussion of the AUID. Each AUID is assigned a new ID hierarchically. Oligotypes for UCYN-A are always used when available. NCD/cyano ID are consider next followed by Genome879 tax ID that pass new threshold supplied to the script (default is 97.0%). If none of these specific ID are available, the more general nifH cluster is used. 

A main function executes and manages the entire process of the script. Below is a breakdown of its functionality:

    1. Parse Command-Line Arguments:
        - Uses argparse to parse the command-line arguments provided by the user.
        - Requires the input path to the annotation table, output path for the updated table, and optionally, the minimum threshold percentage identity for Genome879 ID.
    
    2. Read Input Data:
        - Calls the 'read_tsv()' function to read the input annotation table into a pandas DataFrame.
        - Handles errors related to file not found or empty DataFrame.

    3. Clean Columns:
        - Invokes the 'clean_columns()' function to perform necessary cleanup operations on the DataFrame.
        - Ensures that the DataFrame is properly formatted and ready for consensus ID generation.

    4. Generate Consensus IDs:
        - Utilizes the 'make_consensus_id()' function to create consensus taxonomy IDs based on specified criteria.
        - Consensus IDs are assigned hierarchically, considering various taxonomic IDs and thresholds.

    5. Write Output Data:
        - Calls the 'write_tsv()' function to write the updated DataFrame with consensus IDs to a TSV file.
        - Provides feedback to the user upon successful completion or error during file writing.

Usage:
    python make_consensus_taxonomy.py <path_to_input_table> <path_to_output_table> [--min_pid_genome879 <threshold_pid>] [-h,--help]

Arguments:
    annotation_table: Path to the input annotation table.
    output_table: Path to the output annotation table with consensus ID added.
    --min_pid_genome879: Minimum threshold percentage identity to consider Genome879.id in consensus ID. If not provided, default value is 97.0.

Dependencies:
- pandas

Returns:
    None

"""

import os
import argparse
import sys
import pandas as pd
from argparse import RawTextHelpFormatter


# Get the full path to the currently running script
script_path: str = __file__

# Extract the name of the script for use later
script_name: str = os.path.basename(p=script_path)

# General print statement to show script was called
print(f"{script_name} is currently running...")


# _# Script starts #_#
def setup_argparse() -> argparse.ArgumentParser:
    """
    Set up argparse for command-line argument parsing.

    Returns:
        argparse.ArgumentParser: An ArgumentParser object configured with script-specific command-line argument options.
    """
    parser = argparse.ArgumentParser(
        description=f"""
Overview of {script_name}:
    - reads an input annotation table tsv file
    - cleans up columns in the DataFrame
    - creates a consensus taxonomy ID -> consensus_id
    - writes the updated annotation table with new column consensus_id added
A new column "consensus_id" is added to the input table (auids.annot.tsv), representing the most informative taxonomic ID aggregated across the merged annotation table supplied and new thresholds passed to this script. This ID should be used in any discussion of the AUID. Each AUID is assigned a new ID hierarchically. Oligotypes for UCYN-A are always used when available. NCD/cyano ID are consider next followed by Genome879 tax ID that pass new threshold supplied to the script (default is 97.0%). If none of these specific ID are available, the more general nifH cluster is used. 

	A main function executes and manages the entire process of the script. Below is a breakdown of its functionality:

	1. Parse Command-Line Arguments:
	- Uses argparse to parse the command-line arguments provided by the user.
	- Requires the input path to the annotation table, output path for the updated table, and optionally, the minimum threshold percentage identity for Genome879 ID.

	2. Read Input Data:
	- Calls the 'read_tsv()' function to read the input annotation table into a pandas DataFrame.
	- Handles errors related to file not found or empty DataFrame.

	3. Clean Columns:
	- Invokes the 'clean_columns()' function to perform necessary cleanup operations on the DataFrame.
	- Ensures that the DataFrame is properly formatted and ready for consensus ID generation.

	4. Generate Consensus IDs:
	- Utilizes the 'make_consensus_id()' function to create consensus taxonomy IDs based on specified criteria.
	- Consensus IDs are assigned hierarchically, considering various taxonomic IDs and thresholds.

	5. Write Output Data:
	- Calls the 'write_tsv()' function to write the updated DataFrame with consensus IDs to a TSV file.
	- Provides feedback to the user upon successful completion or error during file writing.

""",
        formatter_class=RawTextHelpFormatter,
        usage=f"{script_name} <path_to_input_table> <path_to_output_table> [--min_pid_genome879 <threshold_pid>] [-h,--help]",
    )
    parser.add_argument(
        "annotation_table",
        help="Path to the input annotation table",
    )
    parser.add_argument(
        "output_table",
        help="Path to output annotation table with consensus ID added",
    )
    parser.add_argument(
        "--min_pid_genome879",
        type=float,
        default=97.0,
        help="Minimum threshold pid to consider Genome879.id in consensus ID",
    )
   
    return parser

def read_tsv(annotation_table: str) -> pd.DataFrame:
    """
    Read the input file (merged annotation table) as a tsv and return a pandas DataFrame.

    Args:
        annotation_table (str): Path to the input annotation table.

    Returns:
        pd.DataFrame: The DataFrame read from the input file
    """
    try:
        print(
            f"Attempting to read input annotation table '{os.path.basename(p=annotation_table)}'"
        )

        df: pd.DataFrame = pd.read_csv(
            filepath_or_buffer=annotation_table,
            sep="\t",
        )

        return df

    except FileNotFoundError:
        print(f"Error: File not found at '{annotation_table}'!")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: DataFrame: '{annotation_table}' is empty")
        sys.exit(1)
    except pd.errors.ParserError:
        print(f"Error: There was an error parsing the csv file: '{annotation_table}'.")
        sys.exit(1)


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean up columns in the DataFrame to work with consensus ID function.

    Args:
        df (pd.DataFrame): DataFrame containing the input data.

    Returns:
        pd.DataFrame: The cleaned DataFrame.
    """
    try:
        print("Cleaning columns...")
        if not df.empty or None:
            required_headers: list[str] = [
                "subcluster",
                "cluster",
                "MarineDiazo.description",
                "MarineDiazo.subject",
                "Genome879.pctId",
                "UCYNAoligos.description",
                "Genome879.tax",
            ]
            missing_headers: list[str] = [
                col for col in required_headers if col not in df.columns
            ]
            # Check if there are any missing headers
            if missing_headers:
                missing_headers_str: str = ",".join(missing_headers)
                raise ValueError(
                    f"There is a problem with the expected columns:\n{required_headers}\nMissing columns: {missing_headers_str}.\nPlease check your column headers on the input file"
                )

            # Check if any required column is all NA values
            columns_with_na_values: list[str] = [
                col for col in required_headers if df[col].isna().all()
            ]
            if columns_with_na_values:
                raise ValueError(
                    f"The following required columns contain only NA values: {', '.join(columns_with_na_values)}"
                )

                # Clean up columns
            df["subcluster"] = df["subcluster"].fillna("")
            df["cluster"] = df["cluster"].apply(
                lambda x: str(object=int(x)) if not pd.isna(x) else x
            ).fillna("")
            # Make new ID columns
            df["MarineDiazo.id"] = (
                df["MarineDiazo.description"] + ";" + df["MarineDiazo.subject"]
            )
            df["UCYNAoligos.id"] = (
                "UCYN-" + df["UCYNAoligos.description"].str.split("_").str[1]
            )

            return df

        # If empty
        else:
            raise pd.errors.EmptyDataError(
                f"DataFrame: '{sys.argv[1]}' headers are there but rows are empty!"
            )

    except pd.errors.EmptyDataError as e:
        print(f"Error: {e}")
        print("Contents of input file:")
        print(df)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


def make_consensus_id(
    df: pd.DataFrame, min_pid_genome879: float = 97.0,) -> pd.DataFrame:
    """
    Create the consensus ID in the DataFrame.

    Args:
        df (pd.DataFrame): DataFrame containing the input data.
        min_pid_genome879 (float): Minimum threshold percentage identity to consider Genome879.id in consensus ID.

    Returns:
        pd.DataFrame: The DataFrame with the consensus ID added.
    """
    try:
        print("Creating consensus ID...")
        if not df.empty or None:
            # Create pid flag
            df.loc[df["Genome879.pctId"] >= min_pid_genome879, "Genome879.pid_flag"] = 1
            condition_pid_flag: int = df["Genome879.pid_flag"] == 1

            # Initial new column for consensus ID
            df["consensus_id"] = "unknown"

            # Assign consensus ID
            df["consensus_id"] = (
                df["UCYNAoligos.id"]
                .fillna(df["MarineDiazo.id"])
                .fillna(df["Genome879.tax"].where(condition_pid_flag))
                .fillna(
                    df.apply(
                        lambda row: (
                        "unknown" + row["subcluster"] if (not row["subcluster"] == "")
                        else "unknown" + row["cluster"] if (not row["cluster"] == "")
                        else "unknown"
                        ),
                        axis=1,
                    )
                )
            )

            # Remove pid flag
            df.drop(columns=["Genome879.pid_flag"], inplace=True)

            return df
        else:
            sys.exit(1)
    except Exception as e:
        print(f"Error creating consensus ID: {e}")
        sys.exit(1)


def write_tsv(df: pd.DataFrame, output_path: str) -> bool:
    """
    Write the DataFrame with newly added consensus ID column to a TSV file.

    Args:
        df (pd.DataFrame): DataFrame to be written to the file.
        output_path (str): Path to the output file.
    Returns:
    Boolean. This value that is used in main() function to demonstrate the code has executed properly. If file is written, True is returned. If not file write, False is returned.
    """
    if not df.empty or None:
        try:
            df.to_csv(output_path, sep="\t", index=False)
            print(f"Output file '{output_path}' written as tsv.")
            return True
        except Exception as e:
            print(f"Error writing file: {e}")
            sys.exit(1)
    return False


def main():
    """
    Main function to execute the script.

    This function serves as the entry point for the script execution. It orchestrates the entire process of reading input data, processing it, generating consensus taxonomy IDs, and writing the updated output data. Below is a breakdown of its functionality:

    1. setup_argparse(): Set up argparse for command-line argument parsing.
    2. read_tsv(input_table: str) -> pd.DataFrame: Read the input file (merged annotation table) as a tsv and return a pandas DataFrame.
    3. clean_columns(df: pd.DataFrame) -> pd.DataFrame: Clean up columns in the DataFrame to work with consensus ID function.
    4. make_consensus_id(df: pd.DataFrame, min_pid_genome879: float) -> pd.DataFrame: Create the consensus ID in the DataFrame.
    5. write_tsv(df: pd.DataFrame, output_path: str) -> None: Write the DataFrame to a TSV file.

    Returns:
        None
    """
    # Set success flag to False
    success = False

    try:
        parser: argparse.ArgumentParser = setup_argparse()
        # args: argparse.Namespace = parser.parse_args()
        # Parse the provided arguments
        args, unknown_args = parser.parse_known_args()
        # Check for any unknown arguments and raise an error if found
        if unknown_args:
            raise argparse.ArgumentError(
                argument=None,
                message=f"Unrecognized arguments: {' '.join(unknown_args)}",
            )

        # Check if all required arguments are provided
        if args.annotation_table is None or args.output_table is None:
            parser.print_help()
            sys.exit()
        # return

    except argparse.ArgumentError as e:
        print(f"Error: {e}")
        print("Failed to parse command-line arguments.")
        sys.exit(5)
    except Exception as e:
        print(f"Error: {e}")
        print(f"Failed to execute script: '{script_name}'.")
        sys.exit(11)

    try:
        df: pd.DataFrame = read_tsv(annotation_table=args.annotation_table)
        if df is not df.empty:
            df = clean_columns(df=df)
            df = make_consensus_id(
                df=df,
                min_pid_genome879=args.min_pid_genome879,
            )

            success: bool = write_tsv(df=df, output_path=args.output_table)

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    finally:
			  # success block with print statements
        if success:
            print(f"\nScript '{script_name}' executed successfully!")
        else:
            print(
                f"""\nScript '{script_name}' exited with an error.
No consensus ID was made and no output was written! :(
"""
            )

# Conditional block ensures that the `main()` function is executed only when
# called from the command line
if __name__ == "__main__":
    main()
