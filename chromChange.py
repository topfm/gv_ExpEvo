import sys
import pandas as pd

def change_column_value(input_file, output_file, new_value):
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(input_file, sep="\t", skiprows= 3)
    
    # Change the values in the "#CHROM" column to the new value
    df['#CHROM'] = new_value
    
    # Write the modified DataFrame back to a CSV file
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py input_file new_value output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    new_value = sys.argv[2]
    output_file = sys.argv[3]
    
    change_column_value(input_file, output_file, new_value)
    print("Values changed successfully.")

