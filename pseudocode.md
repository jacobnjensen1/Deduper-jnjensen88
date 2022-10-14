# The problem

PCR in library prep can cause two or more reads to come from the same section of a physical strand of RNA. These duplicates must not be counted, as they would artificially inflate the measured abundance of that gene. Duplicates can be detected by looking for reads that align to the exact same position on a chromosome and have the same UMI. When duplicates have been detected, only one will be retained in the output. This pseudocode will only cover the base assignment, not the challenge goals.

# Examples

- Case 1:
    - input: `test_sorted.sam` sorted copy of test case from template
    - output: `test_out.sam`
- Case 2:
    - covers more cases than Case 1, including soft clipping, not known but no N UMI, and multiple UMIs for one position.
    - input: `test_2_sorted.sam` heavily modified copy of input for case 1.
    - output: `test_2_out.sam`

# Pseudocode

- make a set for known UMIs 
- read through UMI file
    - append each line to the known UMI set
- make a temporary list for storing records
- Read through the sam file line by line, outputting headers
    - If the temporary list is empty, split the current line on tabs and append to temporary list
        - continue reading next line
    - If the current line aligns on the same chromosome and position (accounting for soft clipping) as the stored record(s), split on tabs and append to the temporary list
    - After the chromosome or position changes or the file ends, loop through the stored lines
        - If the list only contains one record, check that it is mapped and in a primary alignment (mapped and primary might not be required, Leslie's test has 0 for all bitflags)
            - If it is, output it and continue reading
            - If not, clear the list and continue
        - Add lines which are mapped and in a primary alignment to a new list (again, might not be required)
        - loop through the new list
            - make a list which only contains the UMIs of the reads
            - make a set from the UMI list
            - get the intersection of the UMI set and the set of known UMIs
            - use list.index() to find the first location of each UMI in the UMI list and write out the corresponding lines
        - Because the chromosome and/or position has changed, clear the temporary list

# Functions

- adjustedPosition()
    - input: list of strings containing one split line
    - output: the listed position - the amount softclipped on the left end