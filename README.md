
# SmartCut

## Overview

This repository contains the source code for the SmartCut algorithm. The SmartCut algorithm is designed to split input DNA sequences into oligonucleotides using a deep learning-based model. These strands are used for polymerase cycling assembly.

## Requirements

The SmartCut algorithm is based on TensorFlow version 2.6.0. All other necessary packages are listed in `requirements.txt`.

The algorithm has been tested and confirmed to work on the following platforms:

* Windows 10/11
* Ubuntu 20.04
* CentOS 7.5

## Usage

To use the SmartCut algorithm, run the `design_genome.py` script. This script applies the SmartCut algorithm to split the sequences of two chromosomes and two genomes.

## Sequencing Results

The `sequencing_results` folder contains sequencing files for the DNA sequences from `data_for_Table1.xlsx` that were designed and synthesized.

## File Structure

* `design_genome.py`: Script to run the SmartCut algorithm for chromosome and genome sequence splitting.
* `sequencing_results`: Folder containing the sequencing results of designed DNA sequences `data_for_Table1.xlsx` .
* `requirements.txt`: List of required Python packages.
