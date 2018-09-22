import sys
import nanoget
from os import path
from argparse import ArgumentParser, FileType
from nanoplot import utils
from nanoplot.filteroptions import filter_and_transform_data
import nanoplotter
import nanoplotter.compplots as compplots
import numpy as np
import logging
from .version import __version__
from nanomath import write_stats


def main():
    '''
    Organization function
    -setups logging
    -gets inputdata
    -calls plotting function
    '''
    args = get_args()
    try:
        utils.make_output_dir(args.outdir)
        utils.init_logs(args, tool="NanoComp")
        args.format = nanoplotter.check_valid_format(args.format)
        settings = vars(args)
        settings["path"] = path.join(args.outdir, args.prefix)
        sources = [args.fastq, args.bam, args.summary, args.fasta]
        sourcename = ["fastq", "bam", "summary", "fasta"]
        if args.split_runs:
            split_dict = validate_split_runs_file(args.split_runs)
        datadf = nanoget.get_input(
            source=[n for n, s in zip(sourcename, sources) if s][0],
            files=[f for f in sources if f][0],
            threads=args.threads,
            readtype=args.readtype,
            names=args.names,
            barcoded=args.barcoded,
            combine="track")
        datadf, settings = filter_and_transform_data(datadf, vars(args))
        if args.raw:
            datadf.to_csv("NanoComp-data.tsv.gz", sep="\t", index=False, compression="gzip")
        if args.split_runs:
            change_identifiers(datadf, split_dict)
        if args.barcoded:
            datadf["dataset"] = datadf["barcode"]
        identifiers = list(datadf["dataset"].unique())
        write_stats(
            datadfs=[datadf[datadf["dataset"] == i] for i in identifiers],
            outputfile=settings["path"] + "NanoStats.txt",
            names=identifiers)
        if args.plot in ['violin', 'box']:
            plots = make_plots(datadf, settings)
            make_report(plots, path.join(args.outdir, args.prefix))
        logging.info("Succesfully processed all input.")
    except Exception as e:
        logging.error(e, exc_info=True)
        raise


def get_args():
    epilog = """EXAMPLES:
    NanoComp --bam alignment1.bam alignment2.bam --outdir compare-runs
    NanoComp --fastq reads1.fastq.gz reads2.fastq.gz reads3.fastq.gz  --names run1 run2 run3
    """
    parser = ArgumentParser(
        description="Compares long read sequencing datasets.",
        epilog=epilog,
        formatter_class=utils.custom_formatter,
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
                        help="Specify the output format of the plots.",
                        default="png",
                        type=str,
                        choices=['eps', 'jpeg', 'jpg', 'pdf', 'pgf', 'png', 'ps',
                                 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff'])
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
                             "'boxplot', 'violinplot' (default) or 'false' (no plots)",
                        type=str,
                        choices=['violin', 'box', 'false'],
                        default="violin")
    visual.add_argument("--title",
                        help="Add a title to all plots, requires quoting if using spaces",
                        type=str,
                        default=None)
    visual.add_argument("--dpi",
                        help="Set the dpi for saving images",
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
    args = parser.parse_args()
    if args.names:
        if not len(args.names) == [len(i) for i in
                                   [args.fastq, args.summary, args.bam, args.fasta] if i][0]:
            sys.exit("ERROR: Number of names (-n) should be same as number of files specified!")
        if len(args.names) != len(set(args.names)):
            sys.stderr.write("\nWarning: duplicate values in -n/--names detected. ")
            sys.stderr.write("Datasets with the same name will be merged.\n\n")
    if args.colors:
        if not len(args.colors) == [len(i) for i in [args.fastq, args.summary, args.bam] if i][0]:
            sys.exit("ERROR: Number of colors (-c) should be same as number of files specified!")
    return args


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


def make_plots(df, settings):
    nanoplotter.plot_settings(dict(), dpi=settings["dpi"])
    df["log length"] = np.log10(df["lengths"])
    if settings["plot"] == "violin":
        violin = True
    else:
        violin = False
    plots = []
    plots.extend(
        compplots.output_barplot(
            df=df,
            figformat=settings["format"],
            path=settings["path"],
            title=settings["title"],
            palette=settings["colors"])
    )
    plots.extend(
        compplots.violin_or_box_plot(
            df=df,
            y="lengths",
            figformat=settings["format"],
            path=settings["path"],
            y_name="Read length",
            violin=violin,
            title=settings["title"],
            palette=settings["colors"])
    )
    plots.extend(
        compplots.violin_or_box_plot(
            df=df,
            y="log length",
            figformat=settings["format"],
            path=settings["path"],
            y_name="Log-transformed read length",
            violin=violin,
            log=True,
            title=settings["title"],
            palette=settings["colors"])
    )
    if "quals" in df:
        plots.extend(
            compplots.violin_or_box_plot(
                df=df,
                y="quals",
                figformat=settings["format"],
                path=settings["path"],
                y_name="Average base call quality score",
                violin=violin,
                title=settings["title"],
                palette=settings["colors"])
        )
    if "percentIdentity" in df:
        plots.extend(
            compplots.violin_or_box_plot(
                df=df[df["percentIdentity"] > np.percentile(df["percentIdentity"], 1)],
                y="percentIdentity",
                figformat=settings["format"],
                path=settings["path"],
                y_name="Percent reference identity",
                violin=violin,
                title=settings["title"],
                palette=settings["colors"])
        )
    if "start_time" in df:
        plots.extend(
            compplots.compare_cumulative_yields(
                df=df,
                path=settings["path"],
                title=settings["title"],
                palette=settings["colors"])
        )
    plots.extend(
        compplots.overlay_histogram(
            df=df,
            path=settings["path"],
            palette=settings["colors"]
        )
    )
    return plots


def make_report(plots, path):
    '''
    Creates a fat html report based on the previously created files
    plots is a list of Plot objects defined by a path and title
    statsfile is the file to which the stats have been saved,
    which is parsed to a table (rather dodgy)
    '''
    logging.info("Writing html report.")
    html_head = """<!DOCTYPE html>
    <html>
        <head>
        <meta charset="UTF-8">
            <style>
            table, th, td {
                text-align: left;
                padding: 2px;
                /* border: 1px solid black;
                border-collapse: collapse; */
            }
            h2 {
                line-height: 0pt;
            }
            </style>
            <title>NanoComp Report</title>
        </head>"""
    html_content = ["\n<body>\n<h1>NanoComp report</h1>"]
    html_content.append("<h2>Summary statistics</h2>")
    with open(path + "NanoStats.txt") as stats:
        html_content.append('\n<table>')
        for line in stats:
            html_content.append('')
            linesplit = line.strip().split('\t')
            if line.startswith('Data'):
                html_content.append('\n<tr></tr>\n<tr>\n\t<td colspan="2">' +
                                    line.strip() + '</td>\n</tr>')
                break
            if len(linesplit) > 1:
                data = ''.join(["<td>" + e + "</td>" for e in linesplit])
                html_content.append("<tr>\n\t" + data + "\n</tr>")
            else:
                html_content.append('\n<tr></tr>\n<tr>\n\t<td colspan="2"><b>' +
                                    line.strip() + '</b></td>\n</tr>')
        for line in stats:
            html_content.append('\n<tr>\n\t<td colspan="2">' +
                                line.strip() + '</td>\n</tr>')
        html_content.append('</table>')
    html_content.append('\n<br>\n<br>\n<br>\n<br>')
    html_content.append("<h2>Plots</h2>")
    for plot in plots:
        html_content.append("\n<h3>" + plot.title + "</h3>\n" + plot.encode())
        html_content.append('\n<br>\n<br>\n<br>\n<br>')
    html_body = '\n'.join(html_content) + "</body></html>"
    html_str = html_head + html_body
    with open(path + "NanoComp-report.html", "w") as html_file:
        html_file.write(html_str)
    return path + "NanoComp-report.html"


if __name__ == '__main__':
    main()
