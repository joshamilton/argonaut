process KMER_FREQ {
    tag "$meta.id"
    label 'process_low'

    input:
    path reads

    output:
    path "genome_estimate", emit: gce2.log

    script:
    """
    /core/projects/EBP/conservation/software/kmerfreq/kmerfreq -k 17 -t 10 read_files.lib

    #calculate kmer number
    #less read_files.lib.kmer.freq.stat | grep "#Kmer indivdual number"
    less read_files.lib.kmer.freq.stat | perl -ne 'next if(/^#/ || /^\s/); print; ' | awk '{print $1"\t"$2}' > read_files.lib.kmer.freq.stat.2colum 

    #homozygous mode
    #/core/projects/EBP/conservation/software/GCE/gce-1.0.2/gce -g 59620023369 -f read_files.lib.kmer.freq.stat.2colum >gce.table 2>gce.log
    #heterozygous mode 
    /core/projects/EBP/conservation/software/GCE/gce-1.0.2/gce -g 59858472737 -f read_files.lib.kmer.freq.stat.2colum -c 75 -H 1 >gce2.table 2>gce2.log
    """
}