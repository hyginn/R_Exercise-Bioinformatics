# .init.R
# Functions to initialize this session
# Boris Steipe
# ====================================================================

# Create a local copy of myScript.R if that hasn't been done yet.
if (! file.exists("myScript.R")) {
    file.copy("tmp.R", "myScript.R")
}

source(".utilities.R")

file.edit("R_Exercise-Bioinformatics.R")

rm(init)  # not needed in Workspace after we've source()'d this file

# [End]
