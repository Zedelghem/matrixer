for X in fasta/*; do
    # Just the filename, no folder name
    Y=$(echo $X | sed "s/.*\///")
    # Shortest match from the end (.*)
    Z=${Y%.*}
    # Longest match from the end (.*)
    #Z=${Y%%.*}
    seqmagick convert --output-format nexus --alphabet protein $X nexus/$Z.nexus
done