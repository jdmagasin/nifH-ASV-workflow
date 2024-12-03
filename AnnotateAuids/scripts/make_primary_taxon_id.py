#!/usr/bin/env python3
"""
make_primary_taxon_id.py

This module adds a primary_id column to an annotation table and is part of the AnnotateAuids stage of the nifH-ASV-workflow. 

Overview:
    - Reads an input annotation table tsv file
    - Cleans up columns in the DataFrame
    - Creates a primary taxonomy ID -> primary_id
    - Writes the updated annotation table with new column primary_id added

    The script depends on the merged annotation table (auids.annot.tsv) created by the AnnotateAuids stage of the nifH-ASV-workflow.

Dependencies:
    - pandas
    - os
    - sys
    - argparse

Usage:
    python3 make_primary_taxon_id.py <path_to_input_table> <path_to_output_table> [--min_pid_genome879 <threshold_pid>] [-h,--help]

Arguments:
    annotation_table: Path to the input annotation table (as tsv).
    output_table: Path to the output annotation table with primary ID added.
    --min_pid_genome879: Minimum threshold percentage identity to consider. Genome879.id in primary ID. If not provided, default value is 97.0.

Returns:
    None

Exit Codes:
    0: Successful execution
    1: General error (file not found, empty DataFrame, parsing error, etc.)
    5: Failed to parse command-line arguments
    11: Unexpected error during script execution

Input File Format:
    TSV file with required columns: subcluster, cluster, Genome879.pctId, Genome879.tax, MarineDiazo.description, MarineDiazo.subject, UCYNAoligos.description

Output File Format:
    TSV file with all input columns plus an additional 'primary_id' column

Error Handling:
    - Checks for missing or empty required columns
    - Handles file not found, empty DataFrame, and parsing errors
    - Provides informative error messages for various failure scenarios

Author: Michael (Mo) Morando

"""

import os
import argparse
import sys
import pandas as pd
from argparse import RawTextHelpFormatter
from typing import List

__author__ = "Michael (Mo) Morando"
__copyright__ = "Copyright 2023"
__maintainer__ = "Michael (Mo) Morando"
__email__ = "mikemo@gmail.com"
__status__ = "Stable"

# Extract the name of the script for use later
script_name: str = os.path.basename(p=__file__)

def setup_input_args(parser:  argparse.ArgumentParser) -> None:
    """
    Add input arguments to the parser.

    Args:
        parser (argparse.ArgumentParser): The argument parser object to add arguments to.

    Returns:
        None
    """
    parser.add_argument(
        "annotation_table",
        help="Path to the input annotation table as tsv",
    )
    parser.add_argument(
        "output_table",
        help="Path to output annotation table with primary ID added",
    )
    
def setup_optional_args(parser: argparse.ArgumentParser) -> None:
    """
        Add optional arguments to the parser.
    
    Args:
        parser (argparse.ArgumentParser): The argument parser object to add arguments to.

    Returns:
        None
    """
    parser.add_argument(
        "--min_pid_genome879",
        type=float,
        default=97.0,
        help="Minimum threshold pid to consider Genome879.id in primary ID. Default is 97.0 PID.",
    )

def read_tsv(annotation_table: str) -> pd.DataFrame:
    """
    Read the input file (merged annotation table) as a tsv and return a pandas DataFrame.

    Args:
        annotation_table (str): Path to the input annotation table.

    Returns:
        pd.DataFrame: The DataFrame read from the input file.
    
    Raises:
        FileNotFoundError: If the input file is not found
        pd.errors.EmptyDataError: If the DataFrame is empty
        pd.errors.ParserError: If there's an error parsing the TSV file
    """
    print(
        f"Attempting to read input annotation table '{os.path.basename(p=annotation_table)}'"
    )

    try:
        df: pd.DataFrame = pd.read_csv(
            filepath_or_buffer=annotation_table,
            sep="\t",
        )
        
        return df

    except FileNotFoundError:
        print(f"Error: File not found at '{annotation_table}'!")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: Problem reading DataFrame: '{annotation_table}' is empty")
        sys.exit(1)
    except pd.errors.ParserError:
        print(f"Error: There was an error parsing the TSV file: '{annotation_table}'.")
        sys.exit(1)


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean up columns in the DataFrame to work with primary ID function.

    Args:
        df (pd.DataFrame): DataFrame containing the input data.

    Returns:
        pd.DataFrame: The cleaned DataFrame.
        
    Raises:
        ValueError: If required columns are missing or contain only NA values
        pd.errors.EmptyDataError: If the DataFrame is empty.
    """
    print("Cleaning columns...")

    try:
        if not df.empty or None:
            required_headers: List[str] = [
                "subcluster",
                "cluster",
                "Genome879.pctId",
                "Genome879.tax",
            ]
            missing_headers: List[str] = [
                col for col in required_headers if col not in df.columns
            ]
            # Check if there are any missing headers
            if missing_headers:
                missing_headers_str: str = ",".join(missing_headers)
                raise ValueError(
                    f"There is a problem with the expected columns:\n{required_headers}\nMissing columns: {missing_headers_str}.\nPlease check your column headers on the input file"
                )

            # Check if any required column is all NA values
            columns_with_na_values: List[str] = [
                col for col in required_headers if df[col].isna().all()
            ]
            if columns_with_na_values:
                raise ValueError(
                    f"Columns contain only NA values: {', '.join(columns_with_na_values)}"
                )

            # Clean up columns
            df["subcluster"] = df["subcluster"].fillna(value="NA")
            df["cluster"] = (
                df["cluster"]
                .apply(lambda x: str(object=int(x)) if not pd.isna(x) else "NA")
                # .fillna(value="NA")
            )
            # Make new ID columns
            df["MarineDiazo.id"] = df.apply(
                lambda row : f"{row['MarineDiazo.description']};{row['MarineDiazo.subject']}" 
                if pd.notna(row["MarineDiazo.description"]) and pd.notna(row["MarineDiazo.subject"]) 
                else pd.NA,
                axis=1
            )
            df["UCYNAoligos.id"] = (
                df["UCYNAoligos.description"]
                .apply(lambda x: f"UCYN-{x.split('_')[1]}" if pd.notna(x) else pd.NA)
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
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)


def make_primary_id(
    df: pd.DataFrame,
    min_pid_genome879: float = 97.0,
) -> pd.DataFrame:
    """
    Create the primary ID in the DataFrame.

    Args:
        df (pd.DataFrame): DataFrame containing the input data.
        min_pid_genome879 (float): Minimum threshold percentage identity to consider Genome879.id in primary ID. 97.0% is used by default.

    Returns:
        pd.DataFrame: The DataFrame with the primary ID added.
    
    Raises:
        Exception: If an error occurs during primary ID creation.
    """
    try:
        print("Creating primary ID...")
        if not df.empty or None:
            # Create pid flag
            df.loc[df["Genome879.pctId"] >= min_pid_genome879, "Genome879.pid_flag"] = 1
            condition_pid_flag: pd.Series[bool] = df["Genome879.pid_flag"] == 1

            # Initial new column for primary ID
            df["primary_id"] = "unknown"

            # Assign primary ID
            df["primary_id"] = (
                df["UCYNAoligos.id"]
                .fillna(df["MarineDiazo.id"])
                .fillna(df["Genome879.tax"].where(condition_pid_flag))
                .fillna(
                    df.apply(
                        lambda row: (
                            "unknown" + row["subcluster"]
                            if (row["subcluster"] != "NA")
                            else (
                                "unknown" + row["cluster"]
                                if (row["cluster"] != "NA")
                                else "unknown"
                            )
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
        print(f"Error creating primary ID: {e}")
        sys.exit(1)


def write_tsv(df: pd.DataFrame, output_path: str) -> bool:
    """
    Write the DataFrame with newly added primary ID column to a TSV file.

    Args:
        df (pd.DataFrame): DataFrame to be written to the file.
        output_path (str): Path to the output file.
    
    Returns:
        bool: True if file is written successfully, False otherwise
    
    Raises:
        Exception: If an error occurs during file writing.
    """
    if not df.empty or None:
        try:
            df.to_csv(path_or_buf=output_path, sep="\t", index=False, na_rep='NA')
            print(f"Output file '{output_path}' written as tsv.")
            return True
        except Exception as e:
            print(f"Error writing file: {e}")
            sys.exit(1)
    return False


def main():
    """
    Main function to execute the script.

    This function serves as the entry point for the script execution. It orchestrates the entire process of reading input data, processing it, generating primary taxonomy IDs, and writing the updated output data. Below is a breakdown of its functionality:

    1. setup_argparse(): Set up argparse for command-line argument parsing.
    2. read_tsv(input_table: str) -> pd.DataFrame: Read the input file (merged annotation table) as a tsv and return a pandas DataFrame.
    3. clean_columns(df: pd.DataFrame) -> pd.DataFrame: Clean up columns in the DataFrame to work with primary ID function.
    4. make_primary_id(df: pd.DataFrame, min_pid_genome879: float) -> pd.DataFrame: Create the primary ID in the DataFrame.
    5. write_tsv(df: pd.DataFrame, output_path: str) -> None: Write the DataFrame to a TSV file.

    Returns:
        None
    
    Raises:
        argparse.ArgumentError: If there are unrecognized command-line arguments.
        Exception: For any unexpected errors during script execution.
        SystemExit: If an error occurs during execution
    """
    # General print statement to show script was called
    print(f"{script_name} is currently running...")
    
    success = False

    try:
        df: pd.DataFrame = read_tsv(annotation_table=args.annotation_table)
        df = clean_columns(df=df)
        df = make_primary_id(
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
No primary ID was made and no output was written! :(
"""
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter,
        usage=f"{script_name} <path_to_input_table> <path_to_output_table> [--min_pid_genome879 <threshold_pid>] [-h,--help]",
    )
    
    try:
        # parser: argparse.ArgumentParser = setup_argparse()
        setup_input_args(parser)
        setup_optional_args(parser)
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

    except argparse.ArgumentError as e:
        print(f"Error: {e}")
        print("Failed to parse command-line arguments.")
        sys.exit(5)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        print(f"Failed to execute script: '{script_name}'.")
        sys.exit(11)

    main()
