#!/bin/bash

# Check if the 'inkscape', 'gs' and 'pdfunite' commands are available
command -v inkscape >/dev/null 2>&1 || { echo >&2 "Inkscape is not installed. Aborting."; exit 1; }
command -v gs >/dev/null 2>&1 || { echo >&2 "Ghostscript is not installed. Aborting."; exit 1; }
command -v pdfunite >/dev/null 2>&1 || { echo >&2 "pdfunite is not installed. Aborting."; exit 1; }

# Iterate over each SVG file in the current directory
for svg_file in *.svg; do
    # Check if there are any SVG files
    if [ -e "$svg_file" ]; then
        # Define the output PDF file name by replacing .svg with .pdf
        pdf_file="${svg_file%.svg}.pdf"
        
        # Convert SVG to PDF using Inkscape
        inkscape --export-filename="$pdf_file" "$svg_file"
        
        # Check if the conversion was successful
        if [ $? -eq 0 ]; then
            echo "Conversion of $svg_file to $pdf_file successful."

            # Optimize the PDF file using Ghostscript to reduce size
            gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/prepress -dNOPAUSE -dQUIET -dBATCH -sOutputFile="${pdf_file%.pdf}_optimized.pdf" "$pdf_file"
            
            # Check if optimization was successful
            if [ $? -eq 0 ]; then
                echo "Optimization of $pdf_file successful."
            else
                echo "Error optimizing $pdf_file."
            fi
        else
            echo "Error converting $svg_file to $pdf_file."
        fi
    else
        echo "No SVG files found in the current directory."
    fi
done

# Merge all the optimized PDF files into a single file
pdfunite *_optimized.pdf merged_output.pdf

# Check if the merging was successful
if [ $? -eq 0 ]; then
    echo "Merging of PDF files successful. Output saved as merged_output.pdf."
else
    echo "Error merging PDF files."
fi

rm *_optimized.pdf
rm *_Figures.pdf
