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
    df["log length"] = np.log10(df["lengths"])
    sub_df["log length"] = np.log10(sub_df["lengths"])
    plots = []
    plots.extend(
        compplots.output_barplot(
            df=df,
            path=settings["path"],
            title=settings["title"],
            figformat=settings["format"])
    )
    plots.extend(
        compplots.n50_barplot(
            df=sub_df,
            path=settings["path"],
            title=settings["title"],
            figformat=settings["format"])
    )
    plots.extend(
        compplots.violin_or_box_plot(
            df=sub_df[sub_df["length_filter"]],
            y="lengths",
            path=settings["path"],
            y_name="Read length",
            plot=settings["plot"],
            title=settings["title"],
            figformat=settings["format"])
    )
    plots.extend(
        compplots.violin_or_box_plot(
            df=sub_df[sub_df["length_filter"]],
            y="log length",
            path=settings["path"],
            y_name="Log-transformed read length",
            plot=settings["plot"],
            log=True,
            title=settings["title"],
            figformat=settings["format"])
    )
    if "quals" in df:
        plots.extend(
            compplots.violin_or_box_plot(
                df=sub_df,
                y="quals",
                path=settings["path"],
                y_name="Average base call quality score",
                plot=settings["plot"],
                title=settings["title"],
                figformat=settings["format"])
        )
    if "duration" in df:
        plots.extend(
            compplots.compare_sequencing_speed(
                df=sub_df,
                path=settings["path"],
                title=settings["title"],
                figformat=settings["format"])
        )
    if "percentIdentity" in df:
        plots.extend(
            compplots.violin_or_box_plot(
                df=sub_df[sub_df["percentIdentity"] > np.percentile(sub_df["percentIdentity"], 1)],
                y="percentIdentity",
                path=settings["path"],
                y_name="Percent reference identity",
                plot=settings["plot"],
                title=settings["title"],
                figformat=settings["format"])
        )
        plots.append(
            compplots.overlay_histogram_identity(
                df=sub_df[sub_df["percentIdentity"] > np.percentile(sub_df["percentIdentity"], 1)],
                path=settings["path"],
                palette=settings["colors"],
                figformat=settings["format"])
        )

        plots.append(
            compplots.overlay_histogram_phred(
                df=sub_df[sub_df["percentIdentity"] > np.percentile(sub_df["percentIdentity"], 1)],
                path=settings["path"],
                palette=settings["colors"],
                figformat=settings["format"])
        )

    if "start_time" in df:
        plots.extend(
            compplots.compare_cumulative_yields(
                df=df,
                path=settings["path"],
                title=settings["title"],
                palette=settings["colors"],
                figformat=settings["format"])
        )
    if "channelIDs" in df:
        plots.append(
            compplots.active_pores_over_time(
                df=df,
                path=settings["path"],
                palette=settings["colors"],
                title=settings["title"],
                figformat=settings["format"]
            )
        )
    plots.extend(
        compplots.overlay_histogram(
            df=sub_df,
            path=settings["path"],
            palette=settings["colors"],
            figformat=settings["format"]
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
            body {margin:0}
.hiddentitle { /* hides titles that are not necessary for content, but are for outline */
  position: absolute;
  width: 1px;
  height: 1px;
  overflow: hidden;
  left: -10000px;
}

h1 { color: #111; font-family: 'Helvetica Neue', sans-serif; font-size: 60px; font-weight: bold; letter-spacing: -1px; line-height: 1; text-align: center; }

h2 { color: #111; font-family: 'Open Sans', sans-serif; font-size: 25px; font-weight: 300; line-height: 32px; text-align: center; padding-bottom: 0;}

h3 { color: #111; font-family: 'Helvetica Neue', sans-serif; font-size: 16px; font-weight: 150; margin: 0 0 0 0; text-align: left; padding:20px 0px 20px 0px;}

table {
  border: none;
  table-layout: auto;
  empty-cells: hide;
  padding: 5px;
  font-family: Arial, Helvetica, sans-serif;
  border-collapse: separate;
  margin-left: auto;
  margin-right: auto;
  overflow-x: auto;
  white-space: nowrap;
  width: 100%;
}

table td, table th {
  border: 1px solid #ddd;
  padding: 8px;
}

table tr:nth-child(even){background-color: #f2f2f2;}

table tr:hover {background-color: #ddd;}

/* Style the button that is used to open and close the collapsible content */
.collapsible {
  background-color: #39CCCC;
  color: white;
  cursor: pointer;
  padding: 18px;
  width: 100%;
  border: none;
  text-align: left;
  outline: none;
  font-size: 15px;
}

/* Add a background color to the button if it is clicked on (add the .active class with JS), and when you move the mouse over it (hover) */
.active, .collapsible:hover {
    color:white;
  background-color: #001f3f;
}

/* Style the collapsible content. Note: hidden by default */
.collapsible-content {
  padding: 0 18px;
  display: block;
  overflow: hidden;
  background-color: #FFFFFF;
}

.collapsible:after {
  content: '-';
  font-size: 20px;
    font-weight: bold;
  float: right;
    color:white;
  margin-left: 5px;
}

.active:after {
  content: '+'; /* Unicode character for "minus" sign (-) */
      color: white;

}

.hiddentitle { /* hides titles that are not necessary for content, but are for outline */
  position: absolute;
  width: 1px;
  height: 1px;
  overflow: hidden;
  left: -10000px;
}

li a, .submenubutton {
    display: inline-block; /* display the list items inline block so the items are vertically displayed */
    color: white;
    text-align: center;
    padding: 14px 16px;
    text-decoration: none; /* removes the underline that comes with the a tag */
}

li a:hover, .submenu:hover .submenubutton { /* when you hover over a submenu item the bkgrnd color is gray */
    background-color: #39CCCC;
}

.submenu {
    display: inline-block; /* idem to above, list items are displayed underneath each other */
}

.submenu-items { /* hides the ul */
    display: none;
    position: absolute;
    background-color: #f9f9f9;
    min-width: 160px;
    z-index: 1;
}

.submenu-items li {
    display: block;
    float: none;
    overflow: hidden;
}

.submenu-items li a { /* styling of the links in the submenu */
    color: black;
    padding: 12px 16px;
    text-decoration: none;
    display: block;
    text-align: left;
}

.submenu-items a:hover {
    background-color: #f1f1f1;
}

.submenu:hover .submenu-items {
    display: block;
    float: bottom;
    overflow: hidden;
}

nav {
    text-align: center;
}

ul {
    border-bottom: 1px solid white;
    font-family: "Trebuchet MS", sans-serif;
    list-style-type: none; /* remove dot symbols from list */
    margin: 0;
    padding: 0;
    overflow: hidden; /* contains the overflow of the element if it goes 'out of bounds' */
    background-color: #001f3f;
    font-size: 1.6em;
}

ul > li > ul {
    font-size: 1em;
}

li {
    float: left; /* floats the list items to the left side of the page */
}

.issue-btn {
  border-right: none;
  float: right;
}

.tablewrapper {
  width: 100%;
  overflow: auto;
}
</style>
            <title>NanoComp Report</title>
        </head>"""

    html_content = []
    html_content.append('<body><nav><ul><li><a href="#stats">Summary Statistics</a></li>')
    html_content.append('<li class="submenu"><a href="#plots" class="submenubtn">Plots</a>')
    html_content.append('<ul class="submenu-items">')
    html_content.extend(['<li><a href="#'
                         + p.title.replace(' ', '_') + '">' + p.title + '</a></li>' for p in plots])
    html_content.append('</ul>')
    html_content.append('</li>')
    html_content.append(
        '<li class="issue-btn"><a href="https://github.com/wdecoster/nanocomp/issues" target="_blank"  class="reporting">Report issue on Github</a></li>')
    html_content.append('</ul></nav>')
    html_content.append("<h1>NanoComp report</h1>")
    html_content.append("<h2 id='stats'>Summary statistics</h2><div class='tablewrapper'>")
    if stats_df is not None:
        html_content.append(stats_df.to_html())
    else:
        html_content.append(utils.stats2html(path + "NanoStats.txt"))
    # html_content.append('\n<br>\n<br>\n<br>\n<br>')
    html_content.append("</div><h2 id='plots'>Plots</h2>")

    for plot in plots:
        html_content.append('<button class="collapsible">' + plot.title + '</button>')
        html_content.append('<section class="collapsible-content"><h4 class="hiddentitle" id="' +
                            plot.title.replace(' ', '_') + '">' + plot.title + '</h4>')
        html_content.append(plot.encode())
        html_content.append('</section>')

    html_content.append(
        '<script>var coll = document.getElementsByClassName("collapsible");var i;for (i = 0; i < coll.length; i++) {'
        'coll[i].addEventListener("click", function() {this.classList.toggle("active");var content = '
        'this.nextElementSibling;if (content.style.display === "none") {content.style.display = "block";} else {'
        'content.style.display = "none";}});}</script>')

    html_body = '\n'.join(html_content) + "</body></html>"
    html_str = html_head + html_body
    with open(path + "NanoComp-report.html", "w") as html_file:
        html_file.write(html_str)
    return path + "NanoComp-report.html"


if __name__ == '__main__':
    main()
