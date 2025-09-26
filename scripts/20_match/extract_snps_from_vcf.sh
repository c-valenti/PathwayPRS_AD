#!/usr/bin/env bash

# Script to extract SNPs from dbSNP VCF using clump file rsIDs
# Usage: ./extract_snps_from_vcf.sh clump_file.txt dbsnp_file.vcf.gz [output_file.vcf]

# Check if required arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <clump_file> <dbsnp_vcf.gz> [output_file]"
    echo "Example: $0 clumped_snps.txt dbsnp_151.vcf.gz extracted_snps.vcf"
    exit 1
fi

CLUMP_FILE="$1"
DBSNP_VCF="$2"
OUTPUT_FILE="${3:-extracted_snps.vcf}"

# Check if input files exist
if [ ! -f "$CLUMP_FILE" ]; then
    echo "Error: Clump file '$CLUMP_FILE' not found!"
    exit 1
fi

if [ ! -f "$DBSNP_VCF" ]; then
    echo "Error: dbSNP VCF file '$DBSNP_VCF' not found!"
    exit 1
fi

# Create temporary file for rsIDs
TEMP_RSIDS=$(mktemp)
trap "rm -f $TEMP_RSIDS" EXIT

echo "Extracting rsIDs from clump file..."

# Extract rsIDs from the clump file (3rd column, skip header)
# Handle the case where there might be spaces or tabs
awk 'NR > 1 && $3 != "" {print $3}' "$CLUMP_FILE" | sort -u > "$TEMP_RSIDS"

# Count how many rsIDs we found
RSID_COUNT=$(wc -l < "$TEMP_RSIDS")
echo "Found $RSID_COUNT unique rsIDs in clump file"

# Show first few rsIDs for verification
echo "First few rsIDs:"
head -5 "$TEMP_RSIDS"

echo "Extracting matching SNPs from dbSNP VCF..."

# Method 1: Using bcftools (recommended if available)
if command -v bcftools &> /dev/null; then
    echo "Using bcftools for extraction..."
    
    # First, copy the header
    bcftools view -h "$DBSNP_VCF" > "$OUTPUT_FILE"
    
    # Then extract SNPs matching our rsIDs
    bcftools view -H "$DBSNP_VCF" | \
    awk -v rsids="$TEMP_RSIDS" '
    BEGIN {
        # Read rsIDs into array
        while ((getline rsid < rsids) > 0) {
            wanted[rsid] = 1
        }
        close(rsids)
    }
    {
        if ($3 in wanted) {
            print $0
        }
    }' >> "$OUTPUT_FILE"

# Method 2: Using zcat and grep (fallback)
else
    echo "bcftools not found, using zcat and grep..."
    
    # Copy header
    zcat "$DBSNP_VCF" | grep "^#" > "$OUTPUT_FILE"
    
    # Create grep pattern from rsIDs
    GREP_PATTERN=$(paste -sd'|' "$TEMP_RSIDS")
    
    # Extract matching lines
    zcat "$DBSNP_VCF" | grep -v "^#" | \
    awk -v pattern="$GREP_PATTERN" '
    BEGIN {
        split(pattern, rsids, "|")
        for (i in rsids) {
            wanted[rsids[i]] = 1
        }
    }
    {
        if ($3 in wanted) {
            print $0
        }
    }' >> "$OUTPUT_FILE"
fi

# Count results
EXTRACTED_COUNT=$(grep -v "^#" "$OUTPUT_FILE" | wc -l)
echo "Extracted $EXTRACTED_COUNT SNPs from dbSNP VCF"

# Find missing SNPs
MISSING_FILE="${OUTPUT_FILE%.vcf}_missing.txt"
FOUND_RSIDS=$(mktemp)
trap "rm -f $TEMP_RSIDS $FOUND_RSIDS" EXIT

echo "Identifying missing SNPs..."

# Get rsIDs that were found
grep -v "^#" "$OUTPUT_FILE" | cut -f3 | sort > "$FOUND_RSIDS"

# Find the difference
comm -23 "$TEMP_RSIDS" "$FOUND_RSIDS" > "$MISSING_FILE"

MISSING_COUNT=$(wc -l < "$MISSING_FILE")

# Show summary
echo "Summary:"
echo "- Input rsIDs: $RSID_COUNT"
echo "- Extracted SNPs: $EXTRACTED_COUNT"
echo "- Missing SNPs: $MISSING_COUNT"
echo "- Output file: $OUTPUT_FILE"
echo "- Missing SNPs file: $MISSING_FILE"

# Show first few extracted SNPs
echo ""
echo "First few extracted SNPs:"
grep -v "^#" "$OUTPUT_FILE" | head -3 | cut -f1-5

# Show first few missing SNPs
if [ $MISSING_COUNT -gt 0 ]; then
    echo ""
    echo "First few missing SNPs:"
    head -5 "$MISSING_FILE"
fi

# Optional: Create a summary file with just the key information
SUMMARY_FILE="${OUTPUT_FILE%.vcf}_summary.txt"
echo "Creating summary file: $SUMMARY_FILE"
echo -e "CHROM\tPOS\tID\tREF\tALT" > "$SUMMARY_FILE"
grep -v "^#" "$OUTPUT_FILE" | cut -f1-5 >> "$SUMMARY_FILE"

echo "Done!"
