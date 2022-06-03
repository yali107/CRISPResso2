from pydantic import BaseModel

class MainInput(BaseModel):


    class Config:
        frozen = True


# parser.add_argument('--version', action='version', version='%(prog)s '+__version__)
# parser.add_argument('-r1', '--fastq_r1', type=str,  help='First fastq file', default='', required='fastq_r1' in requiredParams)
# parser.add_argument('-r2', '--fastq_r2', type=str,  help='Second fastq file for paired end reads', default='')
#
# parser.add_argument('-a', '--amplicon_seq', type=str,  help='Amplicon Sequence (can be comma-separated list of multiple sequences)', required='amplicon_seq' in requiredParams)
#
# parser.add_argument('-an', '--amplicon_name', type=str,  help='Amplicon Name (can be comma-separated list of multiple names, corresponding to amplicon sequences given in --amplicon_seq', default='Reference')
# parser.add_argument('-amas', '--amplicon_min_alignment_score', type=str,  help='Amplicon Minimum Alignment Score; score between 0 and 100; sequences must have at least this homology score with the amplicon to be aligned (can be comma-separated list of multiple scores, corresponding to amplicon sequences given in --amplicon_seq)', default="")
# parser.add_argument('--default_min_aln_score', '--min_identity_score',  type=int, help='Default minimum homology score for a read to align to a reference amplicon', default=60)
# parser.add_argument('--expand_ambiguous_alignments', help='If more than one reference amplicon is given, reads that align to multiple reference amplicons will count equally toward each amplicon. Default behavior is to exclude ambiguous alignments.', action='store_true')
# parser.add_argument('--assign_ambiguous_alignments_to_first_reference', help='If more than one reference amplicon is given, ambiguous reads that align with the same score to multiple amplicons will be assigned to the first amplicon. Default behavior is to exclude ambiguous alignments.', action='store_true')
# parser.add_argument('-g', '--guide_seq', '--sgRNA', help="sgRNA sequence, if more than one, please separate by commas. Note that the sgRNA needs to be input as the guide RNA sequence (usually 20 nt) immediately adjacent to but not including the PAM sequence (5' of NGG for SpCas9). If the PAM is found on the opposite strand with respect to the Amplicon Sequence, ensure the sgRNA sequence is also found on the opposite strand. The CRISPResso convention is to depict the expected cleavage position using the value of the parameter '--quantification_window_center' nucleotides from the 3' end of the guide. In addition, the use of alternate nucleases besides SpCas9 is supported. For example, if using the Cpf1 system, enter the sequence (usually 20 nt) immediately 3' of the PAM sequence and explicitly set the '--cleavage_offset' parameter to 1, since the default setting of -3 is suitable only for SpCas9.", default='')
# parser.add_argument('-gn', '--guide_name', help="sgRNA names, if more than one, please separate by commas.", default='')
# parser.add_argument('-fg', '--flexiguide_seq', help="sgRNA sequence (flexible) (can be comma-separated list of multiple flexiguides). The flexiguide sequence will be aligned to the amplicon sequence(s), as long as the guide sequence has homology as set by --flexiguide_homology.")
# parser.add_argument('-fh', '--flexiguide_homology', type=int, help="flexiguides will yield guides in amplicons with at least this homology to the flexiguide sequence.", default=80)
# parser.add_argument('-fgn', '--flexiguide_name', help="flexiguide name", default='')
# parser.add_argument('--discard_guide_positions_overhanging_amplicon_edge', help="If set, for guides that align to multiple positions, guide positions will be discarded if plotting around those regions would included bp that extend beyond the end of the amplicon. ", action='store_true')
# parser.add_argument('-e', '--expected_hdr_amplicon_seq', help='Amplicon sequence expected after HDR', default='')
# parser.add_argument('-c', '--coding_seq',  help='Subsequence/s of the amplicon sequence covering one or more coding sequences for frameshift analysis. If more than one (for example, split by intron/s), please separate by commas.', default='')
#
# #quality filtering options
# parser.add_argument('-q', '--min_average_read_quality', type=int, help='Minimum average quality score (phred33) to keep a read', default=0)
# parser.add_argument('-s', '--min_single_bp_quality', type=int, help='Minimum single bp score (phred33) to keep a read', default=0)
# parser.add_argument('--min_bp_quality_or_N', type=int, help='Bases with a quality score (phred33) less than this value will be set to "N"', default=0)
#
# #output options
# parser.add_argument('--file_prefix',  help='File prefix for output plots and tables', default='')
# parser.add_argument('-n', '--name',  help='Output name of the report (default: the name is obtained from the filename of the fastq file/s used in input)', default='')
# parser.add_argument('-o', '--output_folder',  help='Output folder to use for the analysis (default: current folder)', default='')
#
# ## read preprocessing params
# parser.add_argument('--split_interleaved_input', '--split_paired_end', help='Splits a single fastq file containing paired end reads into two files before running CRISPResso', action='store_true')
# parser.add_argument('--trim_sequences', help='Enable the trimming of Illumina adapters with Trimmomatic', action='store_true')
# parser.add_argument('--trimmomatic_command', type=str, help='Command to run trimmomatic', default='trimmomatic')
# parser.add_argument('--trimmomatic_options_string', type=str, help='Override options for Trimmomatic, e.g. "ILLUMINACLIP:/data/NexteraPE-PE.fa:0:90:10:0:true"', default='')
# parser.add_argument('--flash_command', type=str, help='Command to run flash', default='flash')
# parser.add_argument('--min_paired_end_reads_overlap',  type=int, help='Parameter for the FLASH read merging step. Minimum required overlap length between two reads to provide a confident overlap. ', default=10)
# parser.add_argument('--max_paired_end_reads_overlap',  type=int, help='Parameter for the FLASH merging step.  Maximum overlap length expected in approximately 90%% of read pairs. Please see the FLASH manual for more information.', default=100)
# parser.add_argument('--stringent_flash_merging', help='Use stringent parameters for flash merging. In the case where flash could merge R1 and R2 reads ambiguously, the expected overlap is calculated as 2*average_read_length - amplicon_length. The flash parameters for --min-overlap and --max-overlap will be set to prefer merged reads with length within 10bp of the expected overlap. These values override the --min_paired_end_reads_overlap or --max_paired_end_reads_overlap CRISPResso parameters.', action='store_true')
# parser.add_argument('--force_merge_pairs', help=argparse.SUPPRESS, action='store_true')#help=Force-merges R1 and R2 if they cannot be merged using flash (use with caution -- may create non-biological apparent indels at the joining site)
#
# #quantification window params
# parser.add_argument('-w', '--quantification_window_size', '--window_around_sgrna', type=str, help='Defines the size (in bp) of the quantification window extending from the position specified by the "--cleavage_offset" or "--quantification_window_center" parameter in relation to the provided guide RNA sequence(s) (--sgRNA). Mutations within this number of bp from the quantification window center are used in classifying reads as modified or unmodified. A value of 0 disables this window and indels in the entire amplicon are considered. Default is 1, 1bp on each side of the cleavage position for a total length of 2bp. Multiple quantification window sizes (corresponding to each guide specified by --guide_seq) can be specified with a comma-separated list.', default='1')
# parser.add_argument('-wc', '--quantification_window_center', '--cleavage_offset', type=str, help="Center of quantification window to use within respect to the 3' end of the provided sgRNA sequence. Remember that the sgRNA sequence must be entered without the PAM. For cleaving nucleases, this is the predicted cleavage position. The default is -3 and is suitable for the Cas9 system. For alternate nucleases, other cleavage offsets may be appropriate, for example, if using Cpf1 this parameter would be set to 1. For base editors, this could be set to -17 to only include mutations near the 5' end of the sgRNA. Multiple quantification window centers (corresponding to each guide specified by --guide_seq) can be specified with a comma-separated list.", default='-3')
# #    parser.add_argument('--cleavage_offset', type=str, help="Predicted cleavage position for cleaving nucleases with respect to the 3' end of the provided sgRNA sequence. Remember that the sgRNA sequence must be entered without the PAM. The default value of -3 is suitable for the Cas9 system. For alternate nucleases, other cleavage offsets may be appropriate, for example, if using Cpf1 this parameter would be set to 1. To suppress the cleavage offset, enter 'N'.", default=-3)
# parser.add_argument('--exclude_bp_from_left', type=int, help='Exclude bp from the left side of the amplicon sequence for the quantification of the indels', default=15)
# parser.add_argument('--exclude_bp_from_right', type=int, help='Exclude bp from the right side of the amplicon sequence for the quantification of the indels', default=15)
# parser.add_argument('--use_legacy_insertion_quantification', help='If set, the legacy insertion quantification method will be used (i.e. with a 1bp quantification window, indels at the cut site and 1bp away from the cut site would be quantified). By default (if this parameter is not set) with a 1bp quantification window, only insertions at the cut site will be quantified.', action='store_true')
#
# parser.add_argument('--ignore_substitutions', help='Ignore substitutions events for the quantification and visualization', action='store_true')
# parser.add_argument('--ignore_insertions', help='Ignore insertions events for the quantification and visualization', action='store_true')
# parser.add_argument('--ignore_deletions', help='Ignore deletions events for the quantification and visualization', action='store_true')
# parser.add_argument('--discard_indel_reads', help='Discard reads with indels in the quantification window from analysis', action='store_true')
#
# # alignment parameters
# parser.add_argument('--needleman_wunsch_gap_open', type=int, help='Gap open option for Needleman-Wunsch alignment', default=-20)
# parser.add_argument('--needleman_wunsch_gap_extend', type=int, help='Gap extend option for Needleman-Wunsch alignment', default=-2)
# parser.add_argument('--needleman_wunsch_gap_incentive', type=int, help='Gap incentive value for inserting indels at cut sites', default=1)
# parser.add_argument('--needleman_wunsch_aln_matrix_loc', type=str, help='Location of the matrix specifying substitution scores in the NCBI format (see ftp://ftp.ncbi.nih.gov/blast/matrices/)', default='EDNAFULL')
# parser.add_argument('--aln_seed_count', type=int, default=5, help=argparse.SUPPRESS)#help='Number of seeds to test whether read is forward or reverse',default=5)
# parser.add_argument('--aln_seed_len', type=int, default=10, help=argparse.SUPPRESS)#help='Length of seeds to test whether read is forward or reverse',default=10)
# parser.add_argument('--aln_seed_min', type=int, default=2, help=argparse.SUPPRESS)#help='number of seeds that must match to call the read forward/reverse',default=2)
#
# #plotting parameters
# parser.add_argument('--plot_histogram_outliers', help="If set, all values will be shown on histograms. By default (if unset), histogram ranges are limited to plotting data within the 99 percentile.", action='store_true')
#
# #allele plot parameters
# parser.add_argument('--plot_window_size', '--offset_around_cut_to_plot',  type=int, help='Defines the size of the window extending from the quantification window center to plot. Nucleotides within plot_window_size of the quantification_window_center for each guide are plotted.', default=20)
# parser.add_argument('--min_frequency_alleles_around_cut_to_plot', type=float, help='Minimum %% reads required to report an allele in the alleles table plot.', default=0.2)
# parser.add_argument('--expand_allele_plots_by_quantification', help='If set, alleles with different modifications in the quantification window (but not necessarily in the plotting window (e.g. for another sgRNA)) are plotted on separate lines, even though they may have the same apparent sequence. To force the allele plot and the allele table to be the same, set this parameter. If unset, all alleles with the same sequence will be collapsed into one row.', action='store_true')
# parser.add_argument('--allele_plot_pcts_only_for_assigned_reference', help='If set, in the allele plots, the percentages will show the percentage as a percent of reads aligned to the assigned reference. Default behavior is to show percentage as a percent of all reads.', action='store_true')
# parser.add_argument('-qwc', '--quantification_window_coordinates', type=str, help='Bp positions in the amplicon sequence specifying the quantification window. This parameter overrides values of the "--quantification_window_center", "--cleavage_offset", "--window_around_sgrna" or "--window_around_sgrna" values. Any indels/substitutions outside this window are excluded. Indexes are 0-based, meaning that the first nucleotide is position 0. Ranges are separted by the dash sign (e.g. "start-stop"), and multiple ranges can be separated by the underscore (_). ' +
#                                                                                   'A value of 0 disables this filter. (can be comma-separated list of values, corresponding to amplicon sequences given in --amplicon_seq e.g. 5-10,5-10_20-30 would specify the 5th-10th bp in the first reference and the 5th-10th and 20th-30th bp in the second reference)', default=None)
# parser.add_argument('--annotate_wildtype_allele', type=str, help='Wildtype alleles in the allele table plots will be marked with this string (e.g. **).', default='')
#
# #verbosity parameters
# parser.add_argument('--keep_intermediate', help='Keep all the  intermediate files', action='store_true')
# parser.add_argument('--dump', help='Dump numpy arrays and pandas dataframes to file for debugging purposes', action='store_true')
# parser.add_argument('--write_detailed_allele_table', help='If set, a detailed allele table will be written including alignment scores for each read sequence.', action='store_true')
# parser.add_argument('--fastq_output', help='If set, a fastq file with annotations for each read will be produced.', action='store_true')
#
# #report style parameters
# parser.add_argument('--max_rows_alleles_around_cut_to_plot',  type=int, help='Maximum number of rows to report in the alleles table plot.', default=50)
# parser.add_argument('--suppress_report',  help='Suppress output report', action='store_true')
# parser.add_argument('--place_report_in_output_folder',  help='If true, report will be written inside the CRISPResso output folder. By default, the report will be written one directory up from the report output.', action='store_true')
# parser.add_argument('--suppress_plots',  help='Suppress output plots', action='store_true')
# parser.add_argument('--write_cleaned_report', action='store_true', help=argparse.SUPPRESS)#trims working directories from output in report (for web access)
#
# #base editor parameters
# parser.add_argument('--base_editor_output', help='Outputs plots and tables to aid in analysis of base editor studies.', action='store_true')
# parser.add_argument('--conversion_nuc_from',  help='For base editor plots, this is the nucleotide targeted by the base editor', default='C')
# parser.add_argument('--conversion_nuc_to',  help='For base editor plots, this is the nucleotide produced by the base editor', default='T')
#
# #prime editing parameters
# parser.add_argument('--prime_editing_pegRNA_spacer_seq', type=str, help="pegRNA spacer sgRNA sequence used in prime editing. The spacer should not include the PAM sequence. The sequence should be given in the RNA 5'->3' order, so for Cas9, the PAM would be on the right side of the given sequence.", default='')
# parser.add_argument('--prime_editing_pegRNA_extension_seq', type=str, help="Extension sequence used in prime editing. The sequence should be given in the RNA 5'->3' order, such that the sequence starts with the RT template including the edit, followed by the Primer-binding site (PBS).", default='')
# parser.add_argument('--prime_editing_pegRNA_extension_quantification_window_size', type=int, help="Quantification window size (in bp) at flap site for measuring modifications anchored at the right side of the extension sequence. Similar to the --quantification_window parameter, the total length of the quantification window will be 2x this parameter. Default: 5bp (10bp total window size)", default=5)
# parser.add_argument('--prime_editing_pegRNA_scaffold_seq', type=str, help="If given, reads containing any of this scaffold sequence before extension sequence (provided by --prime_editing_extension_seq) will be classified as 'Scaffold-incorporated'. The sequence should be given in the 5'->3' order such that the RT template directly follows this sequence. A common value is 'GGCACCGAGUCGGUGC'.", default='')
# parser.add_argument('--prime_editing_pegRNA_scaffold_min_match_length', type=int, help="Minimum number of bases matching scaffold sequence for the read to be counted as 'Scaffold-incorporated'. If the scaffold sequence matches the reference sequence at the incorporation site, the minimum number of bases to match will be minimally increased (beyond this parameter) to disambiguate between prime-edited and scaffold-incorporated sequences.", default=1)
# parser.add_argument('--prime_editing_nicking_guide_seq', type=str, help="Nicking sgRNA sequence used in prime editing. The sgRNA should not include the PAM sequence. The sequence should be given in the RNA 5'->3' order, so for Cas9, the PAM would be on the right side of the sequence", default='')
# parser.add_argument('--prime_editing_override_prime_edited_ref_seq', type=str, help="If given, this sequence will be used as the prime-edited reference sequence. This may be useful if the prime-edited reference sequence has large indels or the algorithm cannot otherwise infer the correct reference sequence.", default='')
#
# #special running modes
# parser.add_argument('--crispresso1_mode', help='Parameter usage as in CRISPResso 1', action='store_true')
# parser.add_argument('--dsODN', help='Label reads with the dsODN sequence provided', default='')
# parser.add_argument('--auto', help='Infer amplicon sequence from most common reads', action='store_true')
# parser.add_argument('--debug', help='Show debug messages', action='store_true')
# parser.add_argument('--no_rerun', help="Don't rerun CRISPResso2 if a run using the same parameters has already been finished.", action='store_true')
# parser.add_argument('-p', '--n_processes', type=str, help='Specify the number of processes to use for analysis.\
#     Please use with caution since increasing this parameter will significantly increase the memory required to run CRISPResso. Can be set to \'max\'.', default='1')
#
# #processing of aligned bam files
# parser.add_argument('--bam_input', type=str,  help='Aligned reads for processing in bam format', default='')
# parser.add_argument('--bam_chr_loc', type=str,  help='Chromosome location in bam for reads to process. For example: "chr1:50-100" or "chrX".', default='')
#
# #deprecated params
# parser.add_argument('--save_also_png', default=False, help=argparse.SUPPRESS) #help='Save also .png images in addition to .pdf files') #depreciated -- now pngs are automatically created. Pngs can be suppressed by '--suppress_report'
