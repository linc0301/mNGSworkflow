#!/bin/bash
#find and print the sequence with discription of "5S ribosomal"

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_file output_file"
    exit 1
fi

input_file="$1"
output_file="$2"

awk '
BEGIN { found = 0; desc_line = ""; sequence = "" }
/^>/ {
    if (found) {  
        print desc_line > "'$output_file'"
        print sequence > "'$output_file'"
        found = 0
        sequence = ""  
    }
    if ($0 ~ /5S rib/) {  
        found = 1
        desc_line = $0  
    }
    next 
}
{
    if (found) { 
        sequence = (sequence ? sequence ORS : "") $0 
    }
}
END {
    if (found) { 
        print desc_line > "'$output_file'"
        print sequence > "'$output_file'"
    }
}' "$input_file"