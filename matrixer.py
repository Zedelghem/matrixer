#!/usr/bin/python3
import sys, copy, os.path
###########################################################################
# PARSING THE ARGUMENTS
###########################################################################

# Loading parameters from the command line
S_ARGUMENTS = " ".join(sys.argv[1:])

# Setting variables accordingly
# Some cleaning
params = [x.strip() for x in S_ARGUMENTS.split("-") if x != ""]

# Check if all the parameters exist on the list of legal parameters.
# If something is wrong, kill the script.
legal_parameters = ["a", "alignment", "t", "typedata", "f", "features", "m", "multiplications", "b", "binarize", "d", "destination", 'e', 'exportdata']
if not all([x.split()[0] in legal_parameters for x in params]):
    print("One of the parameters you provided does not exist. Do consult the README file.")
    sys.exit()

# Check if the required parameters are present.
# If not, kill the script.
param_names = [x.split()[0] for x in params]
if (("a" not in param_names and
    "alignment" not in param_names) or
    ("f" not in param_names and
    "features" not in param_names)):
    print("One of the required arguments is missing. Make sure that you provide the path to alignment and to the features. Do consult the README file.")
    sys.exit()

# Make a dictionary with parameter values
arguments = {}
for param in params:
    # I take the first letter of the parameter name as the dictionary key for that parameter
    # Therefore, if any parameter is duplicated in UNIX and GNU, the last value is taken.
    param_name = param.split()[0][0]
    param_values = param.split()[1:]
    arguments[param_name] = param_values

# If multiplications has not been set, input default values (1)
if "m" not in arguments.keys():
    arguments["m"] = [1]*len(arguments["f"])
# Else: try to make integers out of the input values and check if there is a multiplications factor for every feature matrix.
else:
    try:
        arguments["m"] = [int(x) for x in arguments["m"]]
    except:
        print("Multiplication factors you provided are not integers.")
        sys.exit()

    if len(arguments["m"]) != len(arguments["f"]):
        print("The number of multiplication factors need to be equal to the number of feature matrices you provided.")
        sys.exit()

# If binarize has not been set, set default value (F).
if "b" not in arguments.keys():
    arguments["b"] = [False]*len(arguments["f"])
# Else: make booleans out of the input values and check if the length of binarize is equal to the length of features.
else:
    temp_bn = []
    for bn in arguments["b"]:
        if bn == "F" or bn == "False" or bn == "0":
            temp_bn.append(False)
        elif bn == "T" or bn == "True" or bn =="1":
            temp_bn.append(True)
    arguments["b"] = temp_bn

    if len(arguments["b"]) != len(arguments["f"]):
        print("The number of binarization statements must be equal to the number of feature matrices you provided. Also, make sure you did not make any typos in the input.")
        sys.exit()

# Check the destination folder.
if "d" not in arguments.keys():
    arguments["d"] = "."
else:
    if not os.path.exists(arguments["d"][0]):
        print("The destination directory you provided does not exist. Pick another one.")
        sys.exit()

# Check if exportdata has been flagged
if "e" not in arguments.keys():
    export = False
else:
    # If flagged: require datatypes after -e (restriction, standard, protein, dna)
    export = True
    feature_types = arguments["e"]

# Check if alignment datatype different than protein
if "t" not in arguments.keys():
    alignment_datatype = "protein"
elif len(arguments["t"]) > 1:
    print("Too many alignment datatypes provided. After -t there must be only one argument.")
    sys.exit()
else:
    if arguments["t"][0].lower() not in ["protein", "dna"]:
        print("Wrong alignment datatype.")
        sys.exit()
    else:
        alignment_datatype = arguments["t"][0]

###########################################################################
# MAIN PART
############################################################################

# Load the alignment
try:
    with open(arguments["a"][0], 'r') as plik:
        alpha_plik = [x.split("\n") for x in plik.read().split(">")]
        alignment = []

        for line in alpha_plik:
            # Deleting / and - from the PF name to make sure the nexus file is generated without ' '
            alignment.append([">" + line[0].replace("/", "_").replace("-", "_"), "".join(line[1:])])
except:
    print("Could not parse the alignment file. Make sure the path is right and that the file is properly formatted.")
    sys.exit()

#print(alignment)

# Load the feature matrices
try:
    features = []
    for f_matrix in arguments["f"]:
        with open(f_matrix, 'r') as plik:
            pre_combined = [x.strip().split(",") for x in plik.readlines()]
            combined = [[line[0], "".join(line[1:])] for line in pre_combined if line[0] != ""]
        features.append(combined)
except:
    print("Could not parse the feature matrix files. Make sure the path is right and that the file is properly formatted.")
    sys.exit()


# Splits a string into n-long chunks
def split_len(seq, length):
    return [seq[i:i+length] for i in range(0, len(seq), length)]

# Binarize the indicated feature matrices
features_ready = []
for ind, feature in enumerate(features):
    if arguments["b"][ind]:
        features_ready.append([ [ x[0], "".join([ z if z == "0" else "1" for z in split_len(x[1], 1) ]) ] for x in feature ])
    else:
        features_ready.append(feature)

# Adding Universal matrix to the tail of a core fasta file
# Works also for combined marix
# And for all other matrices
#"""
def addMatrixAtTail(fastaList, uniMatrix, multiplier):
    INTola_uni_core_fasta = []
    for line in fastaList:
        if line[0].startswith(">PF"):
            for mLine in uniMatrix:
                if line[0].startswith(">" + mLine[0]):
                    INTola_uni_core_fasta.append([line[0], line[1] + (multiplier * mLine[1])])
        else:
            INTola_uni_core_fasta.append(line)
    
    return INTola_uni_core_fasta
#"""

def saveFasta(fname, fastaList):
    with open(fname + ".fasta", 'w') as plik:
        for line in fastaList[1:]:
            #if not line[0].startswith(">"):
            #    line[0] = ">" + line[0]
            plik.write("\n".join(line))
            plik.write("\n")


######################
# COMBINE EVERYTHING #
######################

dump = copy.deepcopy(alignment)

for ind, feat in enumerate(features_ready):
    dump = addMatrixAtTail(dump, feat, arguments["m"][ind])

###################
# SAVING TO FILE  #
###################

# Generating descriptive filename
# Combine feature filenames with their respective multiplications and binarizations
nameparts = []
for ind, feat in enumerate(arguments["f"]):
    if arguments["b"][ind]:
        binarized = "-binarized"
    else:
        binarized = ""

    if arguments["m"][ind] > 1:
        multip = str(arguments["m"][ind]) + "_"
    else:
        multip = ""

    nameparts.append(multip + feat.split(".")[0] + binarized)

# Save to file
file_name_blueprint = arguments["a"][0].split(".")[0] + "+" + "+".join(nameparts)
saveFasta(arguments["d"][0] + "/" + file_name_blueprint, dump)

# Convert to nexus
if export:
    os.system("mkdir nexus")
    os.system("seqmagick convert --output-format nexus --alphabet protein " + arguments["d"][0] + "/" + file_name_blueprint + ".fasta" + " nexus/" + file_name_blueprint + ".nexus")

    # Modify nexus header
    with open("nexus/" + file_name_blueprint + ".nexus", "r") as nexus_raw:
        nexus_file = nexus_raw.readlines()

    # Store info about length of a line of every feature matrix
    feature_lengths = [len(x[0][1]) for x in features]

    # Create a list of tuples (dataformat, feature line length * multiplication factor)
    feature_info = []
    for ind, feat in enumerate(features):
        feature_info.append(tuple([feature_types[ind], feature_lengths[ind] * arguments["m"][ind]]))
    
    #print(feature_info)

    format_line = nexus_file[3].split()
    alignment_seq_len = len(alignment[1][1])
    #print(alignment_seq_len)

    # Generate the feature part of the datatype declaration
    # If a feature has binarize set to True, then enforce datatype "restriction"
    feature_ranges = []
    for ind, feat in enumerate(feature_info):
        if ind > 0:
            feature_ranges.append(feat[0] + ": " + str(int(feature_ranges[-1].split("-")[-1]) + 1) + "-" + str(int(feature_ranges[-1].split("-")[-1]) + feat[1]))
        else:
            feature_ranges.append(feat[0] + ": " + str(alignment_seq_len + 1) + "-" + str(alignment_seq_len + feat[1]))
    
    #print(feature_ranges)

    format_line[1] = "datatype=mixed(" + alignment_datatype + ": 1-" + str(alignment_seq_len) + ", " + ", ".join(feature_ranges) + ")"

    # Append the line to the header
    nexus_file[3] = "	" + " ".join(format_line) + "\n"

    # Overwrite the nexus file
    with open("nexus/" + file_name_blueprint + ".nexus", "w") as nexus_out:
        nexus_out.write("".join(nexus_file))