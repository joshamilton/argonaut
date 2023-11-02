process COVERAGE_LR {
    label 'process_low'

    input:
    val genome_size
    path total_bases

    output:
    path("estimatedCoverageLR.txt")        , emit: coverage_est

    script: 
    """
    printf "%.2f\n" \$((10**2 * $total_bases/$genome_size))e-2 > estimatedCoverage.txt
    """
}