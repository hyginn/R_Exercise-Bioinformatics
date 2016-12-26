# R_Exercise-Bioinformatics.R
#
# Purpose:
#
# Version: 1.0
#
# Date:    2016  12  26
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.0    First code
#
#
# ==============================================================================
#
#           R E A D    T H I S    B E F O R E    Y O U    B E G I N .
#
#
#
# == STRUCTURE OF THIS PROJECT =================================================

# When this project loads, ".Rprofile" is executed, defines a function init()
# and prompts you to execute the function. The sole purpose of init() is to
# source() the file ".init.R" at the time you are sitting at the console and
# working. ".init.R" in turn loads utility functions into the Workspace by
# source()'ing ".utilities.R", and then loads the file you are reading now,
# which contains usage information and other notes. The actual contents of the
# project is in the Project Files listed below. There may be additional files
# either in the main working directory or in subfolders "assets" and/or "data" -
# as the case may be. All files and folders are listed in the File Pane and you
# can click on them to explore what else is here.
#

# == HOW TO WORK WITH THE FILES IN THIS PROJECT ================================

# This project contains scripts and data for an introduction to bioinformatics
# with R. The project assumes you have worked through the introductory R
# tutorial (http://steipe.biochemistry.utoronto.ca/abc/index.php/R_tutorial) and
# the exercises in the "R_Exercise-BasicSetup" project. If you haven't done this
# yet, complete it first.

# === Study code in the Project Files, write code in "myScript.R"

# In general you study and execute the code in the major Project Files, which
# you do not modify (i.e. don't edit, save, and commit them). This way you can
# update them from GitHub through R's version control interface if necessary.

# However you should still copy, modify and otherwise experiment and play with
# code, so you build _active_ coding skills. I have provided a file "myScript.R"
# for this purpose, but of course you can create any number of scripts on your
# own (prefacing their names with the prefix "my" may be a good idea). Save your
# code in "myScript.R", but don't _commit_ your changes to version control. The
# reason is the same: once you _commit_ your changes, you can no longer simply
# pull updated versions of the project from the GitHub repository but will
# need to do more complicated things.


# ==== In case of error messages when updating files
#
# This probably won't happen, but if you pull from GitHub and get the following
# type of error ...
#     ---------------
#     error: Your local changes to the following files would be
#     overwritten by merge
#     ...
#     Please commit your changes or stash them before you can merge.
#     ---------------
# ... then, you need to bring the offending file into its original state.
# Open the Commit window, select the file, and click on the Revert button.
#
# Of course, you can save a local copy under a different name before you revert,
# in case you want to keep your changes.


# === Working with the project files

# Work through the project files in the sequence listed below. While working
# with the project files, DO NOT SIMPLY  source()  THEM!

# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask me. Don't continue if you don't understand what's going on.
# That's not how it works ...  Every single expression counts.

# == PROJECT FILES =============================================================

#    1: Data.R
#    2: Sequence.R
#    3: Structure.R
#    4: Phylogeny.R
#    5: Function.R




# [END]
