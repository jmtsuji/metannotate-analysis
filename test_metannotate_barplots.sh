#!/usr/bin/env bash
set -euo pipefail
# test_metannotate_barplots.sh
# Copyright Jackson M. Tsuji, 2018
# Neufeld lab, University of Waterloo, Canada
# Description: Runs automated test of metannotate_barplots.R

# Hard-coded variables
SCRIPT_VERSION="v0.9" # to match git tag
SOURCE_FILES=(rpoB_0_MetagenomeTest_0_annotations_5z4KAl762541689.tsv dataset_info_template_FILLED.tsv hmm_info_template_FILLED.tsv 171011_barplot_05_custom_plot_template_Family_0.01_FILLED.tsv)
TARGET_FILES=(171011_barplot_01_total_normalized_hits.tsv 171011_barplot_02_collapsed_to_Family.tsv 171011_barplot_03_plotting_table_Family_0.01.tsv)

# If no input is provided, exit out and provide help
if [ $# == 0 ]
	then
	printf "$(basename $0): Runs automated test of metannotate_barplots.R - simple check of md5 hashes for key output files.\n"
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
SOURCE_DIR=${WORK_DIR}/testing/inputs
TEST_DIR=${WORK_DIR}/testing/outputs
MD5_DIR=${WORK_DIR}/testing/output_md5

function test_inputs {
	# Description: tests that provided folders and files exist in the proper configuration
	# GLOBAL params: ${MD5_DIR}; ${TEST_DIR}; ${TARGET_FILES}
	# Return: none
	
	if [ ! -d ${WORK_DIR} ]; then
		echo "ERROR: Cannot find user-specified directory at '${WORK_DIR}'. Exiting..."
		exit 1
	fi

	if [ ! -d ${SOURCE_DIR} ]; then
		echo "ERROR: Cannot find input file directory at '${SOURCE_DIR}'. Exiting..."
		exit 1
	fi
	
	for file in ${SOURCE_FILES[@]}; do
		if [ ! -f ${SOURCE_DIR}/${file}.md5 ]; then
			echo "ERROR: Cannot find input file '${file}' in '${SOURCE_DIR}'. Exiting..."
			exit 1
		fi
	done
	
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
	
	local out_name="${TEST_DIR}/barplot_test" # **LINKED to the prefix of the TARGET_FILES[@]
	Rscript metannotate_barplots.R -i ${SOURCE_DIR}/SOURCE_FILES[0] -d ${SOURCE_DIR}/SOURCE_FILES[1] -m ${SOURCE_DIR}/SOURCE_FILES[2] -t ${SOURCE_DIR}/SOURCE_FILES[3] -o ${out_name} > /dev/null

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

