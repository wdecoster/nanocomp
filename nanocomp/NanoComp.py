import pickle
import nanoget
from nanoplot.filteroptions import filter_and_transform_data
import nanocomp.utils as utils
import nanocomp.compplots as compplots
import numpy as np
import logging
from nanomath import write_stats
from nanoplot.utils import subsample_datasets


def main():
    '''
    Organization function
    -setups logging
    -gets inputdata
    -calls plotting function
    '''
    settings, args = utils.get_args()
    try:
        utils.make_output_dir(args.outdir)
        utils.init_logs(args)
        args.format = utils.check_valid_format(args.format)
        sources = {
            "fastq": args.fastq,
            "bam": args.bam,
            "cram": args.cram,
            "summary": args.summary,
            "fasta": args.fasta,
            "ubam": args.ubam,
        }
        if args.split_runs:
            split_dict = utils.validate_split_runs_file(args.split_runs)
        if args.pickle:
            from nanoget import combine_dfs
            datadf = combine_dfs(dfs=[pickle.load(open(p, 'rb')) for p in args.pickle],
                                 names=args.names,
                                 method="track")
        elif args.feather:
            from nanoget import combine_dfs
            from pandas import read_feather
            datadf = combine_dfs([read_feather(p) for p in args.feather],
                                 names=args.names or args.feather,
                                 method="track")
        else:
            datadf = nanoget.get_input(
                source=[n for n, s in sources.items() if s][0],
                files=[f for f in sources.values() if f][0],
                threads=args.threads,
                readtype=args.readtype,
                names=args.names,
                barcoded=args.barcoded,
                combine="track")
        datadf, settings = filter_and_transform_data(datadf, vars(args))
        if args.raw:
            datadf.to_csv(settings["path"] + "NanoComp-data.tsv.gz",
                          sep="\t",
                          index=False,
                          compression="gzip")
        if args.store:
            pickle.dump(
                obj=datadf,
                file=open(settings["path"] + "NanoComp-data.pickle", 'wb'))
        if args.split_runs:
            utils.change_identifiers(datadf, split_dict)
        if args.barcoded:
            datadf["dataset"] = datadf["barcode"]
        identifiers = list(datadf["dataset"].unique())
        stats_df = write_stats(
            datadfs=[datadf[datadf["dataset"] == i] for i in identifiers],
            outputfile=settings["path"] + "NanoStats.txt",
            names=identifiers,
            as_tsv=args.tsv_stats)
        if args.plot != 'false':
            plots = make_plots(datadf, settings)
            make_report(plots, settings["path"], stats_df=stats_df)
        logging.info("Succesfully processed all input.")
    except Exception as e:
        logging.error(e, exc_info=True)
        raise


def make_plots(df, settings):
    sub_df = subsample_datasets(df)
    utils.plot_settings(dict(), dpi=settings["dpi"])
    df["log length"] = np.log10(df["lengths"])
    sub_df["log length"] = np.log10(sub_df["lengths"])
    plots = []
    plots.extend(
        compplots.output_barplot(
            df=df,
            path=settings["path"],
            title=settings["title"])
    )
    plots.extend(
        compplots.n50_barplot(
            df=sub_df,
            path=settings["path"],
            title=settings["title"])
    )
    plots.extend(
        compplots.violin_or_box_plot(
            df=sub_df[sub_df["length_filter"]],
            y="lengths",
            path=settings["path"],
            y_name="Read length",
            plot=settings["plot"],
            title=settings["title"])
    )
    plots.extend(
        compplots.violin_or_box_plot(
            df=sub_df[sub_df["length_filter"]],
            y="log length",
            path=settings["path"],
            y_name="Log-transformed read length",
            plot=settings["plot"],
            log=True,
            title=settings["title"])
    )
    if "quals" in df:
        plots.extend(
            compplots.violin_or_box_plot(
                df=sub_df,
                y="quals",
                path=settings["path"],
                y_name="Average base call quality score",
                plot=settings["plot"],
                title=settings["title"])
        )
    if "duration" in df:
        plots.extend(
            compplots.compare_sequencing_speed(
                df=sub_df,
                path=settings["path"],
                title=settings["title"])
        )
    if "percentIdentity" in df:
        plots.extend(
            compplots.violin_or_box_plot(
                df=sub_df[sub_df["percentIdentity"] > np.percentile(sub_df["percentIdentity"], 1)],
                y="percentIdentity",
                path=settings["path"],
                y_name="Percent reference identity",
                plot=settings["plot"],
                title=settings["title"])
        )
        plots.append(
            compplots.overlay_histogram_identity(
                df=sub_df[sub_df["percentIdentity"] > np.percentile(sub_df["percentIdentity"], 1)],
                path=settings["path"],
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
    if "channelIDs" in df:
        plots.append(
            compplots.active_pores_over_time(
                df=df,
                path=settings["path"],
                palette=settings["colors"],
                title=settings["title"]
            )
        )
    plots.extend(
        compplots.overlay_histogram(
            df=sub_df,
            path=settings["path"],
            palette=settings["colors"]
        )
    )
    return plots


def make_report(plots, path, stats_df):
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
    if stats_df is not None:
        html_content.append(stats_df.to_html())
    else:
        html_content.append(utils.stats2html(path + "NanoStats.txt"))
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
