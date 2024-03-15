import os
import argparse
import pandas as pd

"""
make_consensus_taxonomy.py

This script reads an input annotation table, cleans up columns, creates a consensus taxonomy ID, 
and writes the output annotation table with the consensus ID added.

Usage:
    python make_consensus_taxonomy.py --annotation_table <path_to_input_table> --output_table <path_to_output_table> [--min_pid_genome879 <threshold_pid>]

Arguments:
    --annotation_table: Path to the input annotation table.
    --output_table: Path to the output annotation table with consensus ID added.
    --min_pid_genome879: Minimum threshold percentage identity to consider Genome879.id in consensus ID. 
                          If not provided, default value is 97.0.

"""


def setup_argparse():
    """
    Set up argparse for command-line argument parsing.
    """
    parser = argparse.ArgumentParser(
        description="Assign consensus taxonomy ID for annotation file"
    )
    parser.add_argument("--annotation_table", help="Path to the input annotation table")
    parser.add_argument(
        "--output_table", help="Path to output annotation table with consensus ID added"
    )
    parser.add_argument(
        "--min_pid_genome879",
        type=float,
        default=97.0,
        help="Minimum threshold pid to consider Genome879.id in consensus ID",
    )
    return parser


def read_file(file_path):
    """
    Read the input file and return a pandas DataFrame.
    """
    try:
        df: pd.DataFrame = pd.read_csv(
            filepath_or_buffer=file_path,
            sep="\t",
        )
        return df
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}'!")
        return None
    except pd.errors.EmptyDataError:
        print(f"Error: DataFrame: '{file_path}' is empty")
        return None
    except pd.errors.ParserError:
        print(f"Error: There was an error parsing the csv file: '{file_path}'.")
        return None


def clean_columns(df):
    """
    Clean up columns in the DataFrame.
    """
    if df is not None:
        # Clean up columns.  Make sub/cluster strings with "" for unknown.
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


def make_consensus_id(df, min_pid_genome879=97.0):
    """
    Create the consensus ID in the DataFrame.
    """
    if df is not None:
        # Create pid flag
        df.loc[df["Genome879.pctId"] >= min_pid_genome879, "Genome879.pid_flag"] = 1
        condition_pid_flag = df["Genome879.pid_flag"] == 1

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

    return df


def write_file(df, output_path):
    """
    Write the DataFrame to a CSV file.
    """
    if df is not None:
        try:
            df.to_csv(output_path, sep="\t", index=False)
            print(f"Output file '{output_path}' written successfully!")
            return True
        except Exception as e:
            print(f"Error writing file: {e}")
    return False


def main(
    annotation_table,
    output_table,
    min_pid_genome879=97.0,
):
    """
    Main function to execute the script.
    """
    df = read_file(annotation_table)
    if df is not None:
        df = clean_columns(df)
        df = make_consensus_id(
            df,
            min_pid_genome879,
        )
        success = write_file(df, output_table)
        if success:
            print("Script executed successfully!")
        else:
            print("Script exited with an error.")


if __name__ == "__main__":
    parser = setup_argparse()
    args = parser.parse_args()
    main(
        annotation_table=args.annotation_table,
        output_table=args.output_table,
        min_pid_genome879=args.min_pid_genome879,
    )


#  python make_consensus_taxonomy.py --annotation_table ../auids.annot.tsv --output_table TEMP_TEST_ANNOTATION_CON.csv --min_pid_genome879 97.0
