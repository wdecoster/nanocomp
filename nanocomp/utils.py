import logging
import sys
import os
import pandas as pd
import numpy as np
from datetime import datetime as dt
from time import time
import textwrap as _textwrap
from .version import __version__
from argparse import ArgumentParser, FileType, HelpFormatter


def make_output_dir(path):
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except IOError:
        sys.exit("ERROR: No writing permission to the output directory.")


def chunks(values, chunks):
    if values:
        chunksize = int(len(values) / chunks)
        return([' '.join(values[i:i + chunksize]) for i in range(0, len(values), chunksize)])
    else:
        return [" "] * chunks


def stats2html(statsf):
    '''for legacy stats output files'''
    df = pd.read_csv(statsf, sep=':', header=None, names=['feature', 'value'])
    values = df["value"].str.strip().str.replace('\t', ' ').str.split().replace(np.nan, '')
    num = len(values[0]) or 1
    v = [chunks(i, num) for i in values]
    return pd.DataFrame(v, index=df["feature"]).to_html(header=False)


def init_logs(args, tool="NanoComp"):
    """Initiate log file and log arguments."""
    start_time = dt.fromtimestamp(time()).strftime('%Y%m%d_%H%M')
    logname = os.path.join(args.outdir, args.prefix + tool + "_" + start_time + ".log")
    handlers = [logging.FileHandler(logname)]
    if args.verbose:
        handlers.append(logging.StreamHandler())
    logging.basicConfig(
        format='%(asctime)s %(message)s',
        handlers=handlers,
        level=logging.INFO)
    logging.info('{} {} started with arguments {}'.format(tool, __version__, args))
    logging.info('Python version is: {}'.format(sys.version.replace('\n', ' ')))
    return logname


def validate_split_runs_file(split_runs_file):
    """Check if structure of file is as expected and return dictionary linking names to run_IDs."""
    try:
        content = [l.strip() for l in split_runs_file.readlines()]
        if content[0].upper().split('\t') == ['NAME', 'RUN_ID']:
            return {c.split('\t')[1]: c.split('\t')[0] for c in content[1:] if c}
        else:
            sys.exit("ERROR: Mandatory header of --split_runs tsv file not found: 'NAME', 'RUN_ID'")
            logging.error("Mandatory header of --split_runs tsv file not found: 'NAME', 'RUN_ID'")
    except IndexError:
        sys.exit("ERROR: Format of --split_runs tab separated file not as expected")
        logging.error("ERROR: Format of --split_runs tab separated file not as expected")


def change_identifiers(datadf, split_dict):
    """Change the dataset identifiers based on the names in the dictionary."""
    for rid, name in split_dict.items():
        datadf.loc[datadf["runIDs"] == rid, "dataset"] = name


class CustomHelpFormatter(HelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string

    def _fill_text(self, text, width, indent):
        return ''.join(indent + line for line in text.splitlines(keepends=True))

    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        return _textwrap.wrap(text, 80)


def custom_formatter(prog):
    return CustomHelpFormatter(prog)


def get_args():
    epilog = """EXAMPLES:
    NanoComp --bam alignment1.bam alignment2.bam --outdir compare-runs
    NanoComp --fastq reads1.fastq.gz reads2.fastq.gz reads3.fastq.gz  --names run1 run2 run3
    """
    parser = ArgumentParser(
        description="Compares long read sequencing datasets.",
        epilog=epilog,
        formatter_class=custom_formatter,
        add_help=False)
    general = parser.add_argument_group(
        title='General options')
    general.add_argument("-h", "--help",
                         action="help",
                         help="show the help and exit")
    general.add_argument("-v", "--version",
                         help="Print version and exit.",
                         action="version",
                         version='NanoComp {}'.format(__version__))
    general.add_argument("-t", "--threads",
                         help="Set the allowed number of threads to be used by the script",
                         default=4,
                         type=int)
    general.add_argument("-o", "--outdir",
                         help="Specify directory in which output has to be created.",
                         default=".")
    general.add_argument("-p", "--prefix",
                         help="Specify an optional prefix to be used for the output files.",
                         default="",
                         type=str)
    general.add_argument("--verbose",
                         help="Write log messages also to terminal.",
                         action="store_true")
    general.add_argument("--raw",
                         help="Store the extracted data in tab separated file.",
                         action="store_true")
    general.add_argument("--store",
                         help="Store the extracted data in a pickle file for future plotting.",
                         action="store_true")
    general.add_argument("--tsv_stats",
                         help="Output the stats file as a properly formatted TSV.",
                         action='store_true')
    filtering = parser.add_argument_group(
        title='Options for filtering or transforming input prior to plotting')
    filtering.add_argument("--readtype",
                           help="Which read type to extract information about from summary. \
                             Options are 1D, 2D, 1D2",
                           default="1D",
                           choices=['1D', '2D', '1D2'])
    filtering.add_argument("--maxlength",
                           help="Drop reads longer than length specified.",
                           type=int,
                           metavar='N')
    filtering.add_argument("--minlength",
                           help="Drop reads shorter than length specified.",
                           type=int,
                           metavar='N')
    filtering.add_argument("--barcoded",
                           help="Barcoded experiment in summary format, splitting per barcode.",
                           action="store_true")
    filtering.add_argument("--split_runs",
                           help="File: Split the summary on run IDs and use names in tsv file. "
                                "Mandatory header fields are 'NAME' and 'RUN_ID'.",
                           default=False,
                           type=FileType('r'),
                           metavar="TSV_FILE")
    visual = parser.add_argument_group(
        title='Options for customizing the plots created')
    visual.add_argument("-f", "--format",
                        help="Specify the output format of the plots, which are in addition to the html files",
                        default="png",
                        type=str,
                        choices=['png', 'jpg', 'jpeg', 'webp', 'svg', 'pdf', 'eps', 'json'])
    visual.add_argument("-n", "--names",
                        help="Specify the names to be used for the datasets",
                        nargs="+",
                        metavar="names")
    visual.add_argument("-c", "--colors",
                        help="Specify the colors to be used for the datasets",
                        nargs="+",
                        metavar="colors")
    visual.add_argument("--plot",
                        help="Which plot type to use: "
                             "'box', 'violin' (default), 'ridge' (joyplot) or 'false' (no plots)",
                        type=str,
                        choices=['violin', 'box', 'ridge', 'false'],
                        default="violin")
    visual.add_argument("--title",
                        help="Add a title to all plots, requires quoting if using spaces",
                        type=str,
                        default=None)
    visual.add_argument("--dpi",
                        help="Set the dpi for saving images (deprecated)",
                        type=int,
                        default=100)
    target = parser.add_argument_group(
        title="Input data sources, one of these is required.")
    mtarget = target.add_mutually_exclusive_group(
        required=True)
    mtarget.add_argument("--fasta",
                         help="Data is in (compressed) fasta format.",
                         nargs='+',
                         metavar="file")
    mtarget.add_argument("--fastq",
                         help="Data is in (compressed) fastq format.",
                         nargs='+',
                         metavar="files")
    mtarget.add_argument("--summary",
                         help="Data is in (compressed) summary files generated by \
                               albacore or guppy.",
                         nargs='+',
                         metavar="files")
    mtarget.add_argument("--bam",
                         help="Data is in sorted bam files.",
                         nargs='+',
                         metavar="files")
    mtarget.add_argument("--ubam",
                         help="Data is in one or more unmapped bam file(s).",
                         nargs='+',
                         metavar="file")
    mtarget.add_argument("--cram",
                         help="Data is in one or more sorted cram file(s).",
                         nargs='+',
                         metavar="file")
    mtarget.add_argument("--pickle",
                         help="Data is in one or more pickle file(s) from using NanoComp/NanoPlot.",
                         nargs='+',
                         metavar="file")
    mtarget.add_argument("--feather",
                         help="Data is in one or more feather file(s).",
                         nargs='+',
                         metavar="file")
    args = parser.parse_args()
    sources = [args.fastq, args.summary, args.bam, args.fasta,
               args.ubam, args.cram, args.pickle, args.feather]
    if args.names:
        if not len(args.names) == [len(i) for i in sources if i][0]:
            sys.exit("ERROR: Number of names (-n) should be same as number of files specified!")
        if len(args.names) != len(set(args.names)):
            sys.stderr.write("\nWarning: duplicate values in -n/--names detected. ")
            sys.stderr.write("Datasets with the same name will be merged.\n\n")
    if args.colors:
        if not len(args.colors) == [len(i) for i in sources if i][0]:
            sys.exit("ERROR: Number of colors (-c) should be same as number of files specified!")
    settings = vars(args)
    settings["path"] = os.path.join(args.outdir, args.prefix)
    return settings, args
