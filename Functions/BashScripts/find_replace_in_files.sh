#!/bin/sh
# *****************************************************************************
# find_replace_in_files.sh
# Performs a recurseive, case-sensitive directory find and replace of files
# For case-insensitive, use the-i switch in the grep call
# Uses a "startdirectory" parameter so that you can run it outside of the
#  specified direcotry --- else this script will modify itself!
# Changes the Internal Field Separator (IFS) to carriage return to accomodate
#  file paths that contain whitespace
# *****************************************************************************

# **********************Change Variables Here *********************************
startdirectory="[REDACTED]"
searchterm="[REDACTED]"
replaceterm="[REDACTED]"
# *****************************************************************************

# IFS (Internal Field Separator) default is whitespace, which is a 
# problem if there are spaces in the file path. Change IFS to carriage
# return instead
IFS_backup=$IFS
IFS=$'\n'

echo "*****************************************"
echo "* Find and Replace in Files Version 0.6 *"
echo "*****************************************"

# get file names and store as loop variable $file
for file in $(grep -l -R $searchterm $startdirectory)
do
	# extract modification date and time of file
	MODTIME=`stat -c %Y "$file"`
	HMODTIME=`date -d @"$MODTIME"`
	# make the replacementwith sed and update file
	sed -i "s/$searchterm/$replaceterm/ig" "$file"
	# change the modification date of the file back to the original date
	touch -d @$MODTIME "$file"
	echo "Modified: " "$file"
	echo "Modification date/time: " $HMODTIME "(sec since epoch: "$MODTIME")"
done

# change IFS back to whitespace
IFS=IFS_backup
echo " *** Yay! All Done! *** "
