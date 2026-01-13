#!/bin/sh
# Create a SHA256 hash of all relevant documentation files.
# Useful for checking if the document has changed.

cat ./*.tex $(find Figures -name "*.pdf" | sort) $(find Figures -name "*.png" | sort) doc.bib fassaad.bib  | sha256sum | cut -f1 -d' ' > doc.hash
