#!/usr/bin/env bash
set -euo pipefail
# 96_well_spec_analysis_test.sh
# Copyright Jackson M. Tsuji, 2018
# Neufeld lab, University of Waterloo, Canada
# Created May 31, 2018
# Description: Runs automated test of 96_well_spec_analysis.R

# Hard-coded variables
SCRIPT_VERSION="v0.3-dev" # to match git tag
TARGET_FILES=(example_raw_data.tsv example_unknowns.tsv)

# If no input is provided, exit out and provide help
if [ $# == 0 ]
	then
	printf "$(basename $0): Runs automated test of 96_well_spec_analysis.R - simple check of md5 hashes for key output files.\n"
	printf "Version: ${SCRIPT_VERSION}\n"
	printf "Contact Jackson M. Tsuji (jackson.tsuji@uwaterloo.ca) for bug reports or feature requests.\n\n"
	printf "Usage: $(basename $0) repo_directory\n\n"
	printf "Usage details:\n"
	printf "1. repo_directory: the path to the base directory containing the cloned/downloaded git repo (this script comes inside that folder).\n\n"
	printf "NOTE: This script will check MD5 hashes of:\n"

	for file in ${TARGET_FILES[@]}; do
		printf "* ${file}\n"
	done
	printf "\n"

	exit 1
fi

# Set variable from user input
WORK_DIR=$1

# Hard-coded extensions of user-supplied variable
MD5_DIR=${WORK_DIR}/testing/output_known_md5
TEST_DIR=${WORK_DIR}/testing/output_test

function test_inputs {
	# Description: tests that provided folders and files exist in the proper configuration
	# GLOBAL params: ${MD5_DIR}; ${TEST_DIR}; ${TARGET_FILES}
	# Return: none
	
	if [ ! -d ${WORK_DIR} ]; then
		echo "ERROR: Cannot find user-specified directory at '${WORK_DIR}'. Exiting..."
		exit 1
	fi

	if [ ! -d ${MD5_DIR} ]; then
		echo "ERROR: Cannot find md5sum directory at '${MD5_DIR}'. Exiting..."
		exit 1
	fi
	
	for file in ${TARGET_FILES[@]}; do
		if [ ! -f ${MD5_DIR}/${file}.md5 ]; then
			echo "ERROR: Cannot find md5 file '${file}.md5' in '${MD5_DIR}'. Exiting..."
			exit 1
		fi
	done
		
}

function run_script {
	
	Rscript 96_well_spec_analysis.R -i testing/input/example_raw_plate_data.txt -m testing/input/example_sample_metadata.tsv -o testing/output_test/example > /dev/null

}


function check_md5 {
	# Description: generates md5sum for all test files and compares to test files
	# GLOBAL params: ${TEST_DIR}; ${TARGET_FILES}
	# Return: .md5 file for each entry in ${TARGET_FILES}; variable ${OVERALL_TEST_STATUS} (0 if passed, 1 if at least one file failed)

	# Initialize the ${OVERALL_TEST_STATUS} variable
	OVERALL_TEST_STATUS=0

	for file in ${TARGET_FILES[@]}; do
		# Generate MD5
		md5sum ${TEST_DIR}/${file} > ${TEST_DIR}/${file}.md5
		
		# Check if the MD5 files match. Will be 0 if matching and 1 if not matching
		local test_md5=$(cut -d ' ' -f 1 ${TEST_DIR}/${file}.md5)
		local known_md5=$(cut -d ' ' -f 1 ${MD5_DIR}/${file}.md5)
		
		# Report to user
		if [ ${test_md5} = ${known_md5} ]; then
			echo "* ${file}: MD5 check PASSED."
			
			# Clean up
			rm ${TEST_DIR}/${file}.md5
			
		elif [ ${test_md5} != ${known_md5} ]; then
			echo "* ${file}: MD5 check FAILED. Leaving behind ${TEST_DIR}/${file}.md5 for reference."
			
			# Change ${OVERALL_TEST_STATUS} to 1 if at least one test fails
			OVERALL_TEST_STATUS=1
			
			# For now, do not clean up the md5 if failed.
			
		fi
		
	done

}


function main {
	echo "Running $(basename $0), version ${SCRIPT_VERSION}, on $(date)."
	echo ""
	
	test_inputs
	
	# Run script
	echo "Running R script..."
	mkdir -p ${TEST_DIR}
	run_script
	echo ""

	# Check MD5s
	echo "Testing outputs..."
	check_md5
	echo ""
	
	if [ ${OVERALL_TEST_STATUS} = 0 ]; then
		printf "Overall: PASSED. Intermediate files cleaned up. Good to go!\n\n"
		rm -r ${TEST_DIR} # cleanup
	elif [ ${OVERALL_TEST_STATUS} = 1 ]; then
		printf "Overall: FAILED (see above). Kept intermediate files.\n\n"
	fi

}

main

