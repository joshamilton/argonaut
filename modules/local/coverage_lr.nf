process COVERAGE_LR {
    label 'process_low'

    input:
    val genome_size
    path total_bases

    output:
    path("estimatedCoverageLR.txt")        , emit: coverage_est

    script: 
    def est_size
    """
    est_size = \$(grep -o '[0-9]\\+' $genome_size)
    est_bases = \$(grep -o '[0-9]\\+' $total_bases)

    printf "%.2f\n" \$((10**2 * \$est_bases/\$est_size))e-2 > estimatedCoverage.txt
    """
}