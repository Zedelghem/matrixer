# matrixer
A tool to quickly generate input files for MrBayes

## FAQ
**Q:** What does Matrixer do?

**A:** It prepares fasta files for MrBayes from the alignments and feature matrices.


**Q:** What are the inputs?

**A:** You need an alignment in FASTA format and a feature matrix in a CSV file.

## Parameters:
```
--alignment, -a         [REQUIRED] Path to the FASTA file with the alignment. Only one alignment is legal.
--features, -f          [REQUIRED] Path(s) (seperated with spaces) to the CSV file(s) with the feature matrices.
--multiplications, -m   An integer telling the script how many times you wish to multiply the included feature matrix 
                        (or multiple matrices). If multiple matrices are provided, the order of values of --multiplications
                        corresponds to the order of paths in --features. Default value is 1.
--binarize, -b          T(rue) or F(alse) for every feature file. The order in --binarize corresponds to the
                        order of values in --features and --multiplications. False by default.
--destination, -d       Path to the destination folder where the output will be saved. The folder must exist. The current folder
                        is the default value. You cannot choose the filename. It will be generated based on the name of the alignment 
                        file and the feature files.
```

**CAREFUL:** If by accident you use both the GNU and UNIX parameter definitions, the script will take the last value.

## Examples:

`./matrixer.py -a two_cores_alignments.fasta -f nc_matrix.csv clans_matrix.csv -m 1 10 -b T F -d fasta`

`./matrixer.py -a two_cores_alignments.fasta --features nc_matrix.csv clans_matrix.csv -m 1 10 --binarize T F -destination fasta`
