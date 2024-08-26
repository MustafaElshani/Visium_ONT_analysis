!/bin/bash

# Input BAM file
bam_file="./SR-610_DMSO/tagged.bam"

# Output directory
output_dir="./fastqs"
mkdir -p "$output_dir"

# Excluded reads file
excluded_reads_file="./excluded_reads.csv"
excluded_large_reads_file="./excluded_large_reads.csv"
excluded_mismatched_reads_file="./excluded_mismatched_reads.csv"  # New file for mismatched reads

# Stats file
stats_file="./tagged_stats.txt"

# Output FASTQ files
r1_fastq="$output_dir/SR-610_DMSO_S1_L001_R1_001.fastq"
r2_fastq="$output_dir/SR-610_DMSO_S1_L001_R2_001.fastq"

# Count total number of reads in the BAM file using 35 threads
echo "Counting total number of reads in the BAM file..." | tee "$stats_file"
total_reads=$(samtools view -@ 64 -c "$bam_file")
echo "Total reads in BAM file: $total_reads" | tee -a "$stats_file"

# Initialize counters
processed_reads=0
excluded_reads=0
excluded_large_reads=0
excluded_mismatched_reads=0  # Counter for mismatched reads
barcode_lengths=()
umi_lengths=()

# Initialize progress bar variables
bar_length=50  # Length of the progress bar
# Extract reads from the BAM file and format them into R1 and R2 FASTQ files using 35 threads
echo "Extracting reads from BAM to FASTQ files..." | tee -a "$stats_file"
samtools view -@ 35 "$bam_file" | awk -v r1="$r1_fastq" -v r2="$r2_fastq" -v excluded_reads_file="$excluded_reads_file" -v excluded_large_reads_file="$excluded_large_reads_file" -v excluded_mismatched_reads_file="$excluded_mismatched_re$
BEGIN {
    excluded_reads = 0;
    excluded_large_reads = 0;
    excluded_mismatched_reads = 0;  # Initialize mismatched reads counter
    print "ReadID" > excluded_reads_file;
    print "ReadID" > excluded_large_reads_file;
    print "ReadID" > excluded_mismatched_reads_file;  # Initialize mismatched reads file
}
{
    read_id = $1;
    sequence = $10;
    quality = $11;
    cr = ""; ur = ""; cy = ""; uy = "";
 # Extract the CR, UR, CY, and UY tags (use uncorrected barcode and UMI)
    for(i=12; i<=NF; i++) {
        if ($i ~ /^CR:Z:/) cr = substr($i, 6);  # Uncorrected barcode
        if ($i ~ /^UR:Z:/) ur = substr($i, 6);  # Uncorrected UMI
        if ($i ~ /^CY:Z:/) cy = substr($i, 6);  # Quality scores for uncorrected barcode
        if ($i ~ /^UY:Z:/) uy = substr($i, 6);  # Quality scores for uncorrected UMI
    }

    # Track barcode and UMI lengths
    barcode_lengths[length(cr)]++;
    umi_lengths[length(ur)]++;

    # Exclude reads without 16bp barcode or 12bp UMI
    if (length(cr) != 16 || length(ur) != 12) {
        print read_id >> excluded_reads_file;
        excluded_reads++;
        next;
    }

    # Trim UMI and corresponding quality scores to 10 bases 
    #NB THIS MAY NOT BE Required if -r1-length 26 is used in spaceranger count
    ur = substr(ur, 1, 10);
    uy = substr(uy, 1, 10);
 # Ensure total R1 length is 26bp (16bp barcode + 10bp UMI)
    if (length(cr) + length(ur) != 26 || length(cy) != 16 || length(uy) != 10) {
        print read_id >> excluded_reads_file;
        excluded_reads++;
        next;
    }

    # Exclude sequences longer than 3000 bases in R2
    if (length(sequence) > 3000) {
        print read_id >> excluded_large_reads_file;
        excluded_large_reads++;
        next;
    }

    # Ensure the sequence and quality scores match in length for R2
    if (length(sequence) != length(quality)) {
        print read_id >> excluded_mismatched_reads_file;  # Log mismatched reads
        excluded_mismatched_reads++;
        next;
    }
   # Write R1 entry with uncorrected barcode and trimmed UMI
    printf("@%s 1:N:0:GATAATACCG+TTTACGTGGT\n%s%s\n+\n%s%s\n", read_id, cr, ur, cy, uy) >> r1;

    # Write R2 entry
    printf("@%s 2:N:0:GATAATACCG+TTTACGTGGT\n%s\n+\n%s\n", read_id, sequence, quality) >> r2;

    # Update the progress
    processed_reads++;
    if (processed_reads % 1000 == 0) {
        progress = int((processed_reads / total_reads) * bar_length);
        printf("\r[%-*s] %d%%", bar_length, substr("==================================================", 1, progress), int((processed_reads / total_reads) * 100)) > "/dev/stderr";
    }
}
END {
    printf("\nExcluded reads (incorrect barcode or UMI): %d\n", excluded_reads) > "/dev/stderr";
    printf("Excluded large reads (>3000 bases or quality mismatch): %d\n", excluded_large_reads) > "/dev/stderr";
    printf("Excluded mismatched reads (sequence and quality length mismatch): %d\n", excluded_mismatched_reads) > "/dev/stderr";  # Report mismatched reads
    
    print "\nBarcode" > stats_file;
    for (len in barcode_lengths) {
        printf("Length %d: Uncorrected Barcode = %d\n", len, barcode_lengths[len]) >> stats_file;
    }

    print "\nUMI" >> stats_file;
    for (len in umi_lengths) {
        printf("Length %d: Uncorrected UMI = %d\n", len, umi_lengths[len]) >> stats_file;
    }

    printf("\nExcluded reads (incorrect barcode or UMI): %d\n", excluded_reads) >> stats_file;
    printf("Excluded large reads (>3000 bases or quality mismatch): %d\n", excluded_large_reads) >> stats_file;
    printf("Excluded mismatched reads (sequence and quality length mismatch): %d\n", excluded_mismatched_reads) >> stats_file;  # Log mismatched reads to stats

    printf("\r[%-50s] 100%%\n", "==================================================") > "/dev/stderr";  # Ensure the progress bar completes at 100%
}'

# Report the number of reads that were successfully processed
echo "Counting the number of reads in the R1 FASTQ file..." | tee -a "$stats_file"
processed_reads=$(wc -l < "$r1_fastq")
processed_reads=$((processed_reads / 4))
echo "Number of reads with all required tags in the FASTQ: $processed_reads" | tee -a "$stats_file"
# Gzipping the FASTQ files with multiple cores using pigz
echo "Gzipping the FASTQ files with multiple cores..." | tee -a "$stats_file"
pigz -p 64 "$r1_fastq"
pigz -p 64 "$r2_fastq"
echo "Gzipping completed." | tee -a "$stats_file"

# Final message
echo "FASTQ files have been successfully generated and compressed in the $output_dir directory." | tee -a "$stats_file"
echo "Excluded reads have been saved to $excluded_reads_file, $excluded_large_reads_file, and $excluded_mismatched_reads_file" | tee -a "$stats_file"
