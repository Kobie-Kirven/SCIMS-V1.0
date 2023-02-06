import sys
import pysam
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from numpy import median

def read_scaffolds(scaffolds_file):
    scaffolds = []
    with open(scaffolds_file, "r") as f:
        for line in f:
            scaffolds.append(line.strip())
    return scaffolds

def parse_sam(sam_file, scaffolds):
    scaffolds = read_scaffolds(scaffolds)
    start_stop = {}
    depth = {}
    samfile = pysam.AlignmentFile(sam_file, "r")
    for read in samfile.fetch():

        # Skip unmapped, secondary, and supplementary reads
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        # Get the chromosome name
        chrom = samfile.getrname(read.reference_id)

        # Skip scaffolds not in the list
        if chrom not in scaffolds:
            continue

        # Get the start and end positions of where the read maps
        start = read.reference_start - 1
        end = read.reference_end - 1
        
        # Loop through each of the positions in the read and add 1 to the depth
        if chrom in depth:
            # If the read starts before the start of the chromosome, add 0s to the start of the depth list  
            if start < start_stop[chrom][0]:
                depth[chrom] = [0] * (start_stop[chrom][0] - start) + depth[chrom]
                start_stop[chrom][0] = start
                
            if end > start_stop[chrom][1]:
                depth[chrom] += [0] * (end - start_stop[chrom][1])
                start_stop[chrom][1] = end

            # Add 1 to the depth for each position in the read
            for i in range(start - start_stop[chrom][0], end - start_stop[chrom][0]):
                depth[chrom][i] += 1
        else:
            start_stop[chrom] = [start, end]
            depth[chrom] = [0] * (start_stop[chrom][1] - start_stop[chrom][0])
            for i in range(start - start_stop[chrom][0], end - start_stop[chrom][0]):
                depth[chrom][i] = 1
          
    samfile.close()
    return depth

def average_depth(depth):
    depths = {}
    for chrom in depth:
        total_bases = 0
        covered_bases = 0
        bases = []
        for i in range(len(depth[chrom])):
            if depth[chrom][i] > 0:
                bases.append(depth[chrom][i])
                total_bases += 1
                covered_bases += depth[chrom][i]
        depths[chrom] = covered_bases / total_bases
    return depths

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python average_depth_pysam.py <sam_file>")
        sys.exit(1)
    sam_file = sys.argv[1]
    scaffolds = sys.argv[2]
    sex_chrom = sys.argv[3]
    figure = sys.argv[4]
    depth = parse_sam(sam_file, scaffolds)
    average_depths = average_depth(depth)

    depth = 0
    for avg_depth in average_depths:
        print("Average depth for chromosome {}: {}".format(avg_depth, average_depths[avg_depth]))
    
    # Get proportion of coverage for the autosomes
    autosome_props = []
    for chrom in average_depths:
        for chrom1 in average_depths:
            if chrom != chrom1 and chrom != sex_chrom and chrom1 != sex_chrom and chrom != "NC_000024.10":
                proportion = average_depths[chrom] / average_depths[chrom1]
                autosome_props.append(proportion)
    
    # Get proportion of sex chromsome coverage to autosomes 
    sex_chrom_props = []
    for chromosome in average_depths:
        if chromosome != sex_chrom:
            prop = average_depths[sex_chrom] / average_depths[chromosome]
            sex_chrom_props.append(prop)
    
    y_chrom_props = []
    for chromosome in average_depths:
        if chromosome != "NC_000024.10":
            prop = average_depths["NC_000024.10"] / average_depths[chromosome]
            y_chrom_props.append(prop)

    # Plot the histograms
    plt.hist(autosome_props, bins=10, color="blue", alpha=0.5, label="Autosomes", density=True)
    plt.hist(sex_chrom_props, bins=10, color="orange", alpha=0.5, label="Sex Chromosome", density=True)
    plt.hist(y_chrom_props, bins=10, color="green", alpha=0.5, label="Y Chromosome", density=True)

    test = mannwhitneyu(sex_chrom_props, autosome_props)

    # Add the p-value for mann-whitney u test
    plt.text(s=f"test_stat: {str(test[0])} p-value: {str(test[1])}", fontsize=6, color="black", transform=plt.gca().transAxes, y=0.95, x=0.05, style="italic", weight="bold")

    plt.legend(loc='upper right', fontsize=8)
    plt.savefig(figure, dpi=300, bbox_inches="tight")

    print(f"test_stat: {str(test[0])} p-value: {str(test[1])}")


            


