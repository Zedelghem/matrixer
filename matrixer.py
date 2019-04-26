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
legal_parameters = ["a", "alignment", "f", "features", "m", "multiplications", "b", "binarize", "d", "destination"]
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

###########################################################################
# MAIN PART
############################################################################

# Load the alignment
try:
    with open(arguments["a"][0], 'r') as plik:
        alpha_plik = [x.split("\n") for x in plik.read().split(">")]
        alignment = []

        for line in alpha_plik:
            alignment.append([">" + line[0], "".join(line[1:])])
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
            for mLine in uniMatrix[1:]:
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
saveFasta(arguments["d"][0] + "/" + arguments["a"][0].split(".")[0] + "+" + "+".join(nameparts), dump)