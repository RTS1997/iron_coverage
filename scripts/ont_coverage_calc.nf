//------------------------------------------------------------------------------------------------------------
//     ___                                                                                      __ _          
//    |_ _|     _ _    ___    _ _       o O O   __      ___    __ __    ___      _ _   __ _    / _` |   ___   
//     | |     | '_|  / _ \  | ' \     o       / _|    / _ \   \ V /   / -_)    | '_| / _` |   \__, |  / -_)  
//    |___|   _|_|_   \___/  |_||_|   TS__[O]  \__|_   \___/   _\_/_   \___|   _|_|_  \__,_|   |___/   \___|  
//  _|"""""|_|"""""|_|"""""|_|"""""| {======|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""|_|"""""| 
//  "`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'./o--000'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-' 
//
//------------------------------------------------------------------------------------------------------------

// Declare syntax version
nextflow.enable.dsl=2

// Script to calculate coverage from an assembled file and its reads

// Script based on:
// i. https://github.com/mmolari/morbidostat-genome-analysis
// ii. https://www.biostars.org/p/5165/

// Tools used:
// - minimap2
// - samtools depth

params.reads_file = "reads"
params.gen_file = "genome"
params.output_fld = "output"

// Process to map the reads against the assembled genome using minimap2
process map_reads {
    conda "conda_envs/read_map.yml"

    output:
        path("reads.sam")

    script:
        """
        minimap2 -a -x map-ont ${params.gen_file} ${params.reads_file} > reads.sam
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
        samtools sort reads.sam -o reads.sorted.bam
        """
}

// create an index for the bam file
process index_sorted_reads {
    conda "conda_envs/read_map.yml"

    input:
        path("reads.sorted.bam")

    output:
        path("reads.sorted.bam")

    script:
        """
        samtools index reads.sorted.bam
        """
}

// Process to calculate the average read depth and standard deviation.
process average_read_depth {
    conda "conda_envs/read_map.yml"

    publishDir(
        path: "${params.output_fld}/",
        mode: 'copy',
    )

    input:
        path("reads.sorted.bam")

    output:
        path("results.txt")

    script:
        """
        samtools coverage -A -w 32 reads.sorted.bam > results.txt
        """
}

workflow {
   map_reads | sort_mapped_reads | index_sorted_reads | average_read_depth | view
}
