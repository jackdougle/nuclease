#!/bin/bash
#rust-duk --ref <file> --in <file> --outu <file> --outm <file> --k <int>

intro() {
echo "
Written by Jack Douglass
Last modified August 11, 2025

Purpose:

Input Parameters:

Output Parameters:

Memory Parameters:

Function and usage documentation at https://docs.rs/rust-duk and ./README.md.
Please contact jack.gdouglass@gmail.com if you encounter any problems.
"
}

pushd . > /dev/null

rust-duk "$@"