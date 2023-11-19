import netCDF4 as nc
import numpy as np
import netCDF4
import sys

def main():

    args = sys.argv[1:]
    file_name1 = args[0]  # user provided gridcell csv file
    file_name2 = args[1]

    # Open the two netCDF files
    file1 = nc.Dataset(file_name1)
    file2 = nc.Dataset(file_name2)

    # Get the variables from each file
    variables1 = file1.variables
    variables2 = file2.variables

    # Compare the variables
    for var in variables1:
        if var in variables2:
            # If the variable is in both files, compare the data
            print(var)
            if not np.allclose(np.array(variables1[var][:], dtype=float), np.array(variables2[var][:], dtype=float)):
                print(f'Difference in variable: {var}')
                print(f'File1: {variables1[var][:]}')
                print(f'File2: {variables2[var][:]}')
        else:
            # If the variable is not in the second file, print it
            print(f'Variable {var} is not in the second file')

    # Check for variables in the second file that are not in the first
    for var in variables2:
        if var not in variables1:
            print(f'Variable {var} is not in the first file')

    # Close the files
    file1.close()
    file2.close()

if __name__ == '__main__':
    main()
