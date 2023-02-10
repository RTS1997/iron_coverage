//------------------------------------------------------------------------------------------------------------
//     ___                                                                                      __ _          
//    |_ _|     _ _    ___    _ _       o O O   __      ___    __ __    ___      _ _   __ _    / _` |   ___   
//     | |     | '_|  / _ \  | ' \     o       / _|    / _ \   \ V /   / -_)    | '_| / _` |   \__, |  / -_)  
//    |___|   _|_|_   \___/  |_||_|   TS__[O]  \__|_   \___/   _\_/_   \___|   _|_|_  \__,_|   |___/   \___|  
//  _|"""""|_|"""""|_|"""""|_|"""""| {======|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""| 
//  "`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'./o--000'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-' 
//
//------------------------------------------------------------------------------------------------------------


// Script to calculate coverage from an assembled file and its reads

// Script based on:
// i. https://github.com/mmolari/morbidostat-genome-analysis
// ii. https://www.biostars.org/p/5165/

// Tools used:
// - minimap2
// - samtools depth


// run name
// folder containing the input files
params.input_fld = "test_dataset"

// input files directory
input_dir = file(params.input_fld.replaceAll('/$',''))
assert input_dir.isDirectory()

// results directory
output_dir = "$baseDir/results/${input_dir.getName()}"

// Process to map the reads against the assembled genome using minimap2
process map_reads {

    label 'q30m'
    conda "conda_envs/read_map.yml"

    input:
        path(genome_fa), path(genome_gbk), path(reads)

    output:
        path("reads.sam")

    script:
        """
        minimap2 -a -x map-ont -t ${task.cpus} $genome_fa $reads > reads.sam
        """
}

// Process to sort the mapped reads and create a BAM file
process sort_mapped_reads {

    label 'q30m'
    conda "conda_envs/read_map.yml"

    input:
        path("reads.sam")

    output:
        path("reads.sorted.bam")

    script:
        """
        samtools sort -@ ${task.cpus} reads.sam > reads.sorted.bam
        """
}

// Process to calculate the average read depth and standard deviation.
process average_read_depth {

    label 'q30m'
    conda "conda_envs/read_map.yml"


    publishDir "$output_dir/vial_${vial}/time_${timepoint}/", mode: 'copy'

    input:
        path("reads.sorted.bam")

    output:
        path("reads.sorted.bam")

    script:
        """
        samtools depth reads.sorted.bam  |  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}'
        """
}
