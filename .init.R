# .init.R
# Functions to initialize this session
# Boris Steipe
# ====================================================================

# Create a local copy of myScript.R if that hasn't been done yet.
if (! file.exists("myScript.R")) {
    file.copy("tmp.R", "myScript.R")
}

file.edit("R_Exercise-Bioinformatics.R")

# [End]
