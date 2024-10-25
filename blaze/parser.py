import textwrap
import logging
import argparse
import numpy as np
import multiprocessing as mp
import sys

import blaze.helper as helper
from blaze.config import *

# setup logging
LOG_FORMAT = \
'(%(asctime)s) %(message)s'
DATE_FORMATE = '%d/%m/%Y %H:%M:%S' #'%a, %d %b %Y %H:%M:%S'
logging.basicConfig(format=LOG_FORMAT, datefmt=DATE_FORMATE)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def parse_arg(arg_string=None):
    parser = argparse.ArgumentParser(
    description=textwrap.dedent(helper.bold_text(
    '''BLAZE2 is a tool for demultiplexing 10X single cell long-read RNA-seq data. It takes fastq files as 
    input and output a whitelist of barcodes and a fastq with demultiplexed reads.\n
    ''')),
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    usage='%(prog)s [options] <input fastq filename/directory>')

    # Type functions
    def existing_file(fn):
        if fn is not None:
            helper.check_files_exist(fn)
        return fn


    def get_files_from_dir(fastq_dir):
        helper.check_files_exist(fastq_dir)
        if os.path.isdir(fastq_dir):
            fastq_fns = helper.get_files_by_suffix(
                search_dir = fastq_dir, 
                suffix = ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz'],
                recursive=True)
        elif os.path.isfile(fastq_dir):
            fastq_fns = [fastq_dir]
        else:
            helper.err_msg(f"File type of input file/dir {fastq_fns} is not supported.", printit=True)
            sys.exit(1)
        return fastq_fns

    parser.add_argument('fastq_fns', metavar='<input fastq filename/directory>', 
                        help='Filename of input fastq files. Can be a directory or a single file.',
                        type = get_files_from_dir)  # input files


    parser.add_argument('--output-prefix', type=str, default=DEFAULT_PREFIX,
        help= textwrap.dedent(
        """
        Filename of output files. Note that the output can be directed to a different directory 
        by specifying the path in the prefix. E.g., --output-prefix /path/to/output/prefix
        """))

    parser.add_argument('--output-fastq', '-o', type=str, default=DEFAULT_GRB_OUT_FASTQ,
        help='Filename of output fastq file name. Note that the filename has to end \nwith .fastq, .fq, .fastq.gz or .fq.gz.')
    
    parser.add_argument('--threads', type=int, default=mp.cpu_count()-1,
        help='Number of threads to used.')

    parser.add_argument('--batch-size', type=int, default=1000,
        help='Number of reads process together as a batch. Note that if the specified number \n'
            'larger than the number of reads in each fastq files, all reads in the file will be processed as a batch.')

    parser.add_argument('--overwrite', action='store_true',
        help='By default, BLAZE will skip the steps generating the existing file(s). \n'
        'Specify this option to rerun the steps those steps and overwrite the existing output files.')

    parser.add_argument('--minimal_stdout', action='store_true', 
        help='Minimise the command-line printing.')


    # required name argment for whitelist
    whitelist_arg_required = parser.add_argument_group(helper.bold_text(
        "One of these argument is required to generate the whitelist.Note "
        "that if multiple options are specified, the priority is as follows:"
         " --no-whitelisting > --force-cells > --count-threshold > --expect-cells"))
    whitelist_arg_required.add_argument('--expect-cells',type=int,
        help='Expected number of cells.')
    whitelist_arg_required.add_argument('--count-threshold', type=int,
        help='Output the whitelist base of the count threshold of high-quality putative barcodes') 
    whitelist_arg_required.add_argument('--force-cells', type=int,
        help='Force the number of cells to be the specified number. This option is useful when \nthe expected number of cells is known.')
    whitelist_arg_required.add_argument('--no-whitelisting', action='store_false', dest='do_whitelisting',
        help='Do not perform whitelisting, if dumultiplexing is enabled, the reads will be assigned to know list of barcodes specified by --known-bc-list.')


    # Optional argument for whitelist
    whitelist_arg_opt = parser.add_argument_group(
        helper.bold_text("\nOptional arguments for generating the whitelist")
        )
    whitelist_arg_opt.add_argument('--minQ', type=int, default=DEFAULT_GRB_MIN_SCORE,
        help='Minimum phred score for all bases in a putative BC to define a "high quality putative barcode".')
    whitelist_arg_opt.add_argument('--max-edit-distance', type=int, default=DEFAULT_ASSIGNMENT_ED,
        help='Maximum edit distance allowed between a putative barcode and a barcode \nfor a read to be assigned to the barcode.')                    
    whitelist_arg_opt.add_argument('--10x-kit-version', '--kit-version', dest="kit_version", choices=['3v4', '3v3', '3v2', '5v3', '5v2'], default=DEFAULT_GRB_KIT,
        help='Choose from 10X Single Cell 3ʹ gene expression v4, v3, v2 (3v4, 3v3, 3v2) or 5ʹ gene expression v3, v2 (5v3, 5v2). If using other protocols, \n'
            'please do not specify this option and specify --full-bc-whitelist and --umi-len instead.') 
    whitelist_arg_opt.add_argument('--umi-len', dest="umi_len", type=int, default=DEFAULT_UMI_SIZE,
        help='UMI length, will only be used when --kit-version is not specified.')
    whitelist_arg_opt.add_argument('--full-bc-whitelist', 
        type=lambda x: x if helper.check_files_exist(x) else None, 
        default=None,
        help='Filename of the full barcode whitelist. If not specified, the corresponding version of 10X whitelist will be used.')
    whitelist_arg_opt.add_argument('--high-sensitivity-mode', action='store_true', 
        help='Turn on the sensitivity mode, which increases the sensitivity of barcode \n'
            'detections but potentially increase the number false/uninformative BC in \n'
            'the whitelist. Identification of empty droplets are recommended in downstream\n')


    # Demultiplexing argument
    demux_option_opt = parser.add_argument_group(helper.bold_text("Demultiplexing options"))
    demux_option_opt.add_argument('--no-demultiplexing', action='store_false', dest='do_demultiplexing',
        help='Do not perform the demultiplexing step.')
    demux_option_opt.add_argument('--known-bc-list', type=existing_file, default=None,
        help='A file specifies a list of barcodes for demultiplexing. If not specified, the barcodes will be assigned to the whitelist from the whitelisting step.')
    demux_option_opt.add_argument('--no-restrand', dest="restrand", action='store_false',
        help='By default, blaze2 re-strands all reads to transcript strand: \n'
        'reads from the reverse strand (those with ployT instead of polyA) will be reverse complemented \n'
        'the their quality scores will be reversed. This option will disable the re-stranding.')
    ###############################
    ####### checking the argument:
    ###############################
    if arg_string is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(arg_string.split())

    # config the output filenames
    args.output_fastq = args.output_prefix + args.output_fastq
    args.out_raw_bc_fn = args.output_prefix + DEFAULT_GRB_OUT_RAW_BC
    args.out_whitelist_fn = args.output_prefix + DEFAULT_GRB_OUT_WHITELIST
    args.out_emptydrop_fn = args.output_prefix + DEFAULT_EMPTY_DROP_FN
    args.out_plot_fn = args.output_prefix + DEFAULT_KNEE_PLOT_FN
    args.summary_fn = args.output_prefix + DEFAULT_BC_STAT_FN
    # validate output filename 
    if not helper.check_suffix(args.output_fastq, ['.fastq', '.fq', '.fastq.gz', '.fq.gz']):
        helper.err_msg(
            f"Error in filename configuration:"
            f"Filename '{args.output_fastq}' should end with '.fastq', '.fq', '.fastq.gz' or '.fq.gz'. Please check the config.py file in BLAZE.", printit=True)
        sys.exit(1)
    if not helper.check_suffix(args.out_raw_bc_fn, '.csv'):
        helper.err_msg(
            f"Error in filename configuration:"
            f"Filename '{args.out_raw_bc_fn}' should end with '.csv'. Please check the config.py file in BLAZE.",  printit=True)
        sys.exit(1)
    if not helper.check_suffix(args.out_whitelist_fn, '.csv'):
        helper.err_msg(
            f"Error in filename configuration:"
            f"Filename '{args.out_whitelist_fn}' should end with '.csv'. Please check the config.py file in BLAZE.", printit=True)
        sys.exit(1)
    if not helper.check_suffix(args.out_emptydrop_fn, '.csv'):
        helper.err_msg(
            f"Error in filename configuration:"
            f"Filename '{args.out_emptydrop_fn}' should end with '.csv'. Please check the config.py file in BLAZE.", printit=True)
        sys.exit(1)
    if not helper.check_suffix(args.out_plot_fn, '.png'):
        helper.err_msg(
            f"Error in filename configuration:"
            f"Filename '{args.out_plot_fn}' should end with '.png'. Please check the config.py file in BLAZE.", printit=True)
        sys.exit(1)
    
    # create output directory if not exist
    out_dir = os.path.dirname(args.output_fastq)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # config the full whitelist
    if args.full_bc_whitelist:
        helper.warning_msg(textwrap.dedent(
            f'You are using {os.path.basename(args.full_bc_whitelist)} as the full barcode '\
            'whitelist. Note that the barcodes not listed in the file will never be found.'))
    elif args.kit_version == '3v4':
        args.full_bc_whitelist = DEFAULT_GRB_WHITELIST_3V4
        args.umi_len = DEFAULT_UMI_SIZE
    elif args.kit_version == '3v3':
        args.full_bc_whitelist = DEFAULT_GRB_WHITELIST_3V3
        args.umi_len = DEFAULT_UMI_SIZE
    elif args.kit_version == '5v3':
        args.full_bc_whitelist = DEFAULT_GRB_WHITELIST_5V3
        args.umi_len = DEFAULT_UMI_SIZE
    elif args.kit_version == '5v2' or args.kit_version == '3v2':
        args.full_bc_whitelist = DEFAULT_GRB_WHITELIST_V2
        args.umi_len = V2_UMI_SIZE
    else:
        helper.err_msg("Error: Invalid value of --kit-version, please choose from v3 or v2 or specify --full-bc-whitelist.", printit=True) 
        sys.exit(1)

    # Update pipeline args to determine which step to run
    args, pipeline_summary = update_pipeline_args(args)
    if not any([args.do_whitelisting, args.do_bc_search, args.do_demultiplexing]):
        helper.err_msg("No pipeline step will be run. Please check the summary below: \n", printit=True)
        print(pipeline_summary)
        sys.exit(1)

    print(pipeline_summary)


    # check whitelisting argument
    if args.do_whitelisting:
        if not args.expect_cells and not args.count_threshold and not args.force_cells:
            helper.err_msg("One of the following arguments is required: --expect-cells, --count-threshold --no-whitelisting or --force-cells.", printit=True),
            sys.exit(1)


    # hard coded arugument:
    args.empty_max_count=np.inf# count to define the maximum count for empty droplets

    return args


def update_pipeline_args(args):
    """
    Update the pipeline arguments based on the existing output files and the specified arguments
    The argument to check/update:
        1. do_bc_search
        2. do_whitelisting
        3. do_demultiplexing
    Output:
        args: updated args
        pipeline_summary: a summary of the pipeline steps
    """

    # summary of pipeline steps
    pipeline_summary = helper.bold_text('#'*10 + 'Pipeline steps to run:' + '#'*10 + '\n')

    # run bc search?
    if not os.path.exists(args.out_raw_bc_fn) or args.overwrite:
        args.do_bc_search = True
        pipeline_summary += helper.green_msg('Search barcode in reads: Yes\n', printit=False)
        if os.path.exists(args.out_raw_bc_fn) and args.overwrite:
            pipeline_summary += helper.warning_msg(
                f"\tWarning: The output putative barcodes table `{args.out_raw_bc_fn}` exist. It will be overwritten.\n",
                 printit = False)
    
    else:
        args.do_bc_search = False
        pipeline_summary += textwrap.dedent(
        """
        Search barcode in reads: Skipped
            *Reason: The output putative barcodes table exists and --overwrite was not specified.
        """)

    # run whitelisting?
    if not args.do_whitelisting:
        args.out_whitelist_fn = args.known_bc_list
        pipeline_summary += textwrap.dedent(
        f"""
        Generate the whitelist and knee plot: Skipped
           *Reason: --no-whitelisting specified.
        """)
        
    else: 
        assert args.do_whitelisting == True

        whitelist_output_fn = [args.out_emptydrop_fn, args.out_plot_fn, args.out_whitelist_fn]
        output_exists = all([os.path.exists(x) for x in whitelist_output_fn])
        output_time_correct = \
            all([os.path.getmtime(x)>os.path.getmtime(args.out_raw_bc_fn) for x in whitelist_output_fn]) if output_exists else True
        
        if not args.overwrite and output_exists:
            args.do_whitelisting = False

            pipeline_summary += textwrap.dedent(
            f"""
            Generate the whitelist and knee plot: Skipped
                *Reason: All output files from whitelisting step (`{args.out_whitelist_fn}`, `{args.out_emptydrop_fn}` and `{args.out_plot_fn}`) exist and --overwrite was not specified.
            """)
            
            if not output_time_correct:
                pipeline_summary += helper.warning_msg(
                    f"\tWarning: some of these files are older than the putative barcode file {args.out_raw_bc_fn}.\n"
                )

        else:
            args.do_whitelisting = True
            pipeline_summary += helper.green_msg('Generate the whitelist and knee plot: Yes\n', printit=False)
            if args.force_cells is not None:
                pipeline_summary += f"\t*Forced number of cells: {args.force_cells}\n"
                args.count_threshold = None
                args.expect_cells = None
            elif args.count_threshold is not None:
                args.expect_cells = None

            
    # run demultiplexing?
    if not args.do_demultiplexing:
        pipeline_summary += textwrap.dedent(
            f"""
            Demultiplexing (assign reads to cells): Skipped
                *Reason: --no-demultiplexing specified.
            """)

    else:
        args.do_demultiplexing = True
        output_exists = os.path.exists(args.output_fastq) and os.path.exists(args.out_whitelist_fn)
        output_time_correct = \
            os.path.getmtime(args.output_fastq) > os.path.getmtime(args.out_whitelist_fn) if output_exists else True



        if args.overwrite and output_exists:
            args.do_demultiplexing = True
            pipeline_summary += helper.green_msg('Demultiplexing (assign reads to cells): Yes\n', printit=False)
            pipeline_summary += helper.warning_msg(f"\t*NOTE:  {args.output_fastq} will be overwritten.",printit = False)

        elif not args.overwrite and output_exists:
            args.do_demultiplexing = False
            pipeline_summary += textwrap.dedent(
                f"""
                Demultiplexing (assign reads to cells): Skipped
                    *Reason: Output file from the demultiplexing step (`{args.output_fastq}`) exists and --overwrite was not specified.
                """)
            
            if not output_time_correct:
                pipeline_summary += helper.warning_msg(
                    f"\tWarning: some of these files are older than the whitelist {args.out_whitelist_fn}.\n"
                )
        elif not args.out_whitelist_fn:
            args.do_demultiplexing = False
            pipeline_summary += textwrap.dedent(
                f"""
                Demultiplexing (assign reads to cells): Skipped
                    *Reason: No valid whitelist to use.
                """)
        else:
            args.do_demultiplexing = True
            pipeline_summary += helper.green_msg('Demultiplexing (assign reads to cells): Yes\n', printit=False)
            if args.known_bc_list:
                pipeline_summary += f"\t*Barcode list for demultiplexing: {args.out_whitelist_fn}"
        
    pipeline_summary += "\n" + '#'*40 + "\n"
    return args, pipeline_summary

