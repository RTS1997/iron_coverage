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
output_dir = "$baseDir/${input_dir.getName()}"

// Process to map the reads against the assembled genome using minimap2
process map_reads {

    conda "conda_envs/read_map.yml"

    input:
        path(genome_fa), path(genome_gbk), path(reads)

    output:
        path("reads.sam")

    script:
        """
        minimap2 -a -x map-ont $genome_fa $reads > reads.sam
        """
}

// Process to sort the mapped reads and create a BAM file
process sort_mapped_reads {

    conda "conda_envs/read_map.yml"

    input:
        path("reads.sam")

    output:
        path("reads.sorted.bam")

    script:
        """
        samtools sort reads.sam > reads.sorted.bam
        """
}

// Process to calculate the average read depth and standard deviation.
process average_read_depth {

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

// Main workflow for the program, executes the required processes
workflow pileup_workflow {
    main:

        // reads input channel. Has items [reads]
        reads = Channel.fromPath("${input_dir}/reads.fastq.gz")

        // filtered so that only the first timepoint is kept
        assembled_genomes = Channel.fromFilePairs("${input_dir}/*.{fna,gbk}")

        // combine genomes with reads
        combined = assembled_genomes.cross(reads)

        // create symlinks of reference genomes and reads
        create_symlinks(combined)
        
        // map and sort reads
        sorted_reads = map_reads(combined) | sort_mapped_reads

        // create index for sorted reads
        depth_reads = average_read_depth(sorted_reads)

}
