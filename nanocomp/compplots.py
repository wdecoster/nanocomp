from nanoplotter.plot import Plot
from nanoplotter.timeplots import check_valid_time_and_sort
from nanomath import get_N50
import logging
import numpy as np
import plotly
import plotly.graph_objs as go
import sys


def violin_or_box_plot(df, y, path, y_name, settings, title=None, plot="violin", log=False):
    """Create a violin/boxplot/ridge from the received DataFrame.

    The x-axis should be divided based on the 'dataset' column,
    the y-axis is specified in the arguments
    """
    comp = Plot(
        path=f"{path}NanoComp_{y.replace(' ', '_')}_{plot}.html",
        title=f"Comparing {y_name.lower()}",
    )

    if plot == "violin":
        logging.info(f"NanoComp: Creating violin plot for {y}.")

        fig = go.Figure()

        for dataset in df["dataset"].unique():
            fig.add_trace(
                go.Violin(
                    x=df["dataset"][df["dataset"] == dataset],
                    y=df[y][df["dataset"] == dataset],
                    points=False,
                    name=dataset,
                )
            )

        process_violin_and_box(
            fig,
            log=log,
            plot_obj=comp,
            title=title,
            y_name=y_name,
            ymax=np.amax(df[y]),
            settings=settings,
        )

    elif plot == "box":
        logging.info("NanoComp: Creating box plot for {y}.")

        fig = go.Figure()

        for dataset in df["dataset"].unique():
            fig.add_trace(
                go.Box(
                    x=df["dataset"][df["dataset"] == dataset],
                    y=df[y][df["dataset"] == dataset],
                    name=dataset,
                )
            )

        process_violin_and_box(
            fig,
            log=log,
            plot_obj=comp,
            title=title,
            y_name=y_name,
            ymax=np.amax(df[y]),
            settings=settings,
        )

    elif plot == "ridge":
        logging.info("NanoComp: Creating ridges plot for {y}.")

        fig = go.Figure()

        for d in df["dataset"].unique():
            fig.add_trace(go.Violin(x=df[y][df["dataset"] == d], name=d))

        fig.update_traces(orientation="h", side="positive", width=3, points=False)
        fig.update_layout(title=title or comp.title, title_x=0.5)

        comp.fig = fig
        comp.html = comp.fig.to_html(full_html=False, include_plotlyjs="cdn")
        comp.save(settings)

    else:
        logging.error(f"Unknown comp plot type {plot}")
        sys.exit(f"Unknown comp plot type {plot}")

    return [comp]


def process_violin_and_box(fig, log, plot_obj, title, y_name, ymax, settings):
    if log:
        ticks = [10**i for i in range(10) if not 10**i > 10 * (10**ymax)]
        fig.update_layout(
            yaxis=dict(
                tickmode="array",
                tickvals=np.log10(ticks),
                ticktext=ticks,
            )
        )

    fig.update_layout(
        title_text=title or plot_obj.title,
        title_x=0.5,
        yaxis_title=y_name,
    )

    plot_obj.fig = fig
    plot_obj.html = plot_obj.fig.to_html(full_html=False, include_plotlyjs="cdn")
    plot_obj.save(settings)


def output_barplot(df, path, settings, title=None):
    """Create barplots based on number of reads and total sum of nucleotides sequenced."""
    logging.info("NanoComp: Creating barplots for number of reads and total throughput.")
    read_count = Plot(
        path=path + "NanoComp_number_of_reads.html", title="Comparing number of reads"
    )

    counts = df["dataset"].value_counts(sort=False)

    read_count.fig = go.Figure()
    for idx, count in zip(counts.index, counts):
        read_count.fig.add_trace(go.Bar(x=[idx], y=[count], name=idx))

    read_count.fig.update_layout(
        title_text=title or read_count.title,
        title_x=0.5,
        yaxis_title="Number of reads",
    )

    read_count.html = read_count.fig.to_html(full_html=False, include_plotlyjs="cdn")
    read_count.save(settings)

    throughput_bases = Plot(
        path=path + "NanoComp_total_throughput.html",
        title="Comparing throughput in bases",
    )
    if "aligned_lengths" in df:
        throughput = df.groupby("dataset", sort=False)["aligned_lengths"].sum()
        ylabel = "Total bases aligned"
    else:
        throughput = df.groupby("dataset", sort=False)["lengths"].sum()
        ylabel = "Total bases sequenced"

    throughput_bases.fig = go.Figure()
    for idx, sum_dataset in zip(throughput.index, throughput):
        throughput_bases.fig.add_trace(go.Bar(x=[idx], y=[sum_dataset], name=idx))

    throughput_bases.fig.update_layout(
        title=title or throughput_bases.title,
        title_x=0.5,
        yaxis_title=ylabel,
    )

    throughput_bases.html = throughput_bases.fig.to_html(full_html=False, include_plotlyjs="cdn")
    throughput_bases.save(settings)

    return read_count, throughput_bases


def n50_barplot(df, path, settings, title=None):
    """
    Returns Plot object and creates figure(format specified)/html
    containing bar chart of total gb aligned/sequenced read length n50
    """
    n50_bar = Plot(path=path + "NanoComp_N50.html", title="Comparing read length N50")
    datasets = df["dataset"].unique()
    length_column = "aligned_lengths" if "aligned_lengths" in df else "lengths"
    ylabel = "Aligned read length N50" if "aligned_lengths" in df else "Sequenced read length N50"

    n50s = [get_N50(np.sort(df.loc[df["dataset"] == d, length_column])) for d in datasets]
    n50_bar.fig = go.Figure()

    for idx, n50 in zip(datasets, n50s):
        n50_bar.fig.add_trace(go.Bar(x=[idx], y=[n50], name=idx))

    n50_bar.fig.update_layout(
        title=title or n50_bar.title,
        title_x=0.5,
        yaxis_title=ylabel,
    )

    n50_bar.html = n50_bar.fig.to_html(full_html=False, include_plotlyjs="cdn")
    n50_bar.save(settings)
    return [n50_bar]


def compare_sequencing_speed(df, path, settings, title=None):
    logging.info("NanoComp: creating comparison of sequencing speed over time.")
    seq_speed = Plot(
        path=path + "NanoComp_sequencing_speed_over_time.html",
        title="Sequencing speed over time",
    )

    dfs = check_valid_time_and_sort(df, "start_time").set_index("start_time")
    dfs = dfs.loc[dfs["duration"] > 0]

    palette = plotly.colors.DEFAULT_PLOTLY_COLORS * 5

    data = []
    for sample, color in zip(df["dataset"].unique(), palette):
        seqspeed = (
            (
                dfs.loc[dfs["dataset"] == sample, "lengths"]
                / dfs.loc[dfs["dataset"] == sample, "duration"]
            )
            .resample("30T")
            .median()
        )
        data.append(
            go.Scatter(
                x=seqspeed.index.total_seconds() / 3600,
                y=seqspeed,
                opacity=0.75,
                name=sample,
                mode="lines",
                marker=dict(color=color),
            )
        )

    seq_speed.fig = go.Figure({"data": data})

    seq_speed.fig.update_layout(
        title=title or seq_speed.title,
        title_x=0.5,
        xaxis_title="Interval (hours)",
        yaxis_title="Sequencing speed (nucleotides/second)",
    )

    seq_speed.html = seq_speed.fig.to_html(full_html=False, include_plotlyjs="cdn")
    seq_speed.save(settings)
    return [seq_speed]


def compare_cumulative_yields(df, path, settings, palette=None, title=None):
    if palette is None:
        palette = plotly.colors.DEFAULT_PLOTLY_COLORS * 5
    dfs = check_valid_time_and_sort(df, "start_time").set_index("start_time")

    logging.info("NanoComp: Creating cumulative yield plots using {len(dfs)} reads.")
    cum_yield_gb = Plot(
        path=path + "NanoComp_CumulativeYieldPlot_Gigabases.html",
        title="Cumulative yield",
    )
    data = []
    annotations = []
    for sample, color in zip(df["dataset"].unique(), palette):
        cumsum = dfs.loc[dfs["dataset"] == sample, "lengths"].cumsum().resample("10T").max() / 1e9
        data.append(
            go.Scatter(
                x=cumsum.index.total_seconds() / 3600,
                y=cumsum,
                opacity=0.75,
                name=sample,
                marker=dict(color=color),
            )
        )
        annotations.append(
            dict(
                xref="paper",
                x=0.99,
                y=cumsum[-1],
                xanchor="left",
                yanchor="middle",
                text="{}Gb".format(round(cumsum[-1])),
                showarrow=False,
            )
        )

    cum_yield_gb.fig = go.Figure(
        {
            "data": data,
            "layout": go.Layout(
                barmode="overlay",
                title=title or cum_yield_gb.title,
                xaxis=dict(title="Time (hours)"),
                yaxis=dict(title="Yield (gigabase)"),
                annotations=annotations,
            ),
        }
    )
    cum_yield_gb.html = cum_yield_gb.fig.to_html(full_html=False, include_plotlyjs="cdn")
    cum_yield_gb.save(settings)
    return [cum_yield_gb]


def overlay_histogram(df, path, settings, palette=None):
    """
    Use plotly to create an overlay of length histograms
    Return html code, but also save as figure (format specified)

    Only has 10 colors, which get recycled up to 5 times.
    """
    if palette is None:
        palette = plotly.colors.DEFAULT_PLOTLY_COLORS * 5

    hist = Plot(path=path + "NanoComp_OverlayHistogram.html", title="Histogram of read lengths")
    hist.html, hist.fig = plot_overlay_histogram(df, palette, column="lengths", title=hist.title)
    hist.save(settings)

    hist_norm = Plot(
        path=path + "NanoComp_OverlayHistogram_Normalized.html",
        title="Normalized histogram of read lengths",
    )
    hist_norm.html, hist_norm.fig = plot_overlay_histogram(
        df, palette, column="lengths", title=hist_norm.title, density=True
    )
    hist_norm.save(settings)

    log_hist = Plot(
        path=path + "NanoComp_OverlayLogHistogram.html",
        title="Histogram of log transformed read lengths",
    )
    log_hist.html, log_hist.fig = plot_log_histogram(df, palette, title=log_hist.title)
    log_hist.save(settings)

    log_hist_norm = Plot(
        path=path + "NanoComp_OverlayLogHistogram_Normalized.html",
        title="Normalized histogram of log transformed read lengths",
    )
    log_hist_norm.html, log_hist_norm.fig = plot_log_histogram(
        df, palette, title=log_hist_norm.title, density=True
    )
    log_hist_norm.save(settings)

    return [hist, hist_norm, log_hist, log_hist_norm]


def overlay_histogram_identity(df, path, settings, palette=None):
    if palette is None:
        palette = plotly.colors.DEFAULT_PLOTLY_COLORS * 5
    hist_pid = Plot(
        path=path + "NanoComp_OverlayHistogram_Identity.html",
        title="Histogram of percent reference identity",
    )
    hist_pid.html, hist_pid.fig = plot_overlay_histogram(
        df, palette, "percentIdentity", hist_pid.title, density=True
    )
    hist_pid.save(settings)

    return hist_pid


def overlay_histogram_phred(df, path, settings, palette=None):
    """
    Reads with a perfect alignment and thus a percentIdentity of 100
    get a phred score of Inf
    Which is not cool
    So these are set to 60, a very high phred score
    """
    df["phredIdentity"] = -10 * np.log10(1 - (df["percentIdentity"] / 100))
    df["phredIdentity"][np.isinf(df["phredIdentity"])] = 60

    if palette is None:
        palette = plotly.colors.DEFAULT_PLOTLY_COLORS * 5

    hist_phred = Plot(
        path=path + "NanoComp_OverlayHistogram_PhredScore.html",
        title="Histogram of Phred scores",
    )

    hist_phred.html, hist_phred.fig = plot_overlay_histogram(
        df, palette, "phredIdentity", hist_phred.title, bins=20, density=True
    )

    hist_phred.save(settings)

    return hist_phred


def plot_overlay_histogram(df, palette, column, title, bins=None, density=False):
    data = []
    if not bins:
        bins = max(round(int(np.amax(df.loc[:, column])) / 500), 10)

    for d, c in zip(df["dataset"].unique(), palette):
        counts, bins = np.histogram(df.loc[df["dataset"] == d, column], bins=bins, density=density)
        data.append(
            go.Bar(
                x=bins[1:],
                y=counts,
                opacity=0.4,
                name=d,
                text=bins[1:],
                hoverinfo="text",
                hovertemplate=None,
                marker=dict(color=c),
            )
        )

    fig = go.Figure({"data": data, "layout": go.Layout(barmode="overlay", title=title, bargap=0)})

    fig.update_layout(title_x=0.5, yaxis_title="Density" if density else "Number of reads")

    return fig.to_html(full_html=False, include_plotlyjs="cdn"), fig


def plot_log_histogram(df, palette, title, density=False):
    """
    Plot overlaying histograms with log transformation of length
    Return both html and figure
    """
    data = []
    bins = max(round(int(np.amax(df.loc[:, "lengths"])) / 500), 10)
    for d, c in zip(df["dataset"].unique(), palette):
        counts, bins = np.histogram(
            np.log10(df.loc[df["dataset"] == d, "lengths"]), bins=bins, density=density
        )
        data.append(
            go.Bar(
                x=bins[1:],
                y=counts,
                opacity=0.4,
                name=d,
                text=[10**i for i in bins[1:]],
                hoverinfo="text",
                hovertemplate=None,
                marker=dict(color=c),
            )
        )
    xtickvals = [10**i for i in range(10) if not 10**i > 10 * np.amax(df["lengths"])]

    fig = go.Figure(
        {
            "data": data,
            "layout": go.Layout(
                barmode="overlay",
                title=title,
                xaxis=dict(tickvals=np.log10(xtickvals), ticktext=xtickvals),
                bargap=0,
            ),
        }
    )

    fig.update_layout(title_x=0.5, yaxis_title="Density" if density else "Number of reads")

    return fig.to_html(full_html=False, include_plotlyjs="cdn"), fig


def active_pores_over_time(df, path, settings, palette=None, title=None):
    if palette is None:
        palette = plotly.colors.DEFAULT_PLOTLY_COLORS * 5
    dfs = check_valid_time_and_sort(df, "start_time").set_index("start_time")

    logging.info("NanoComp: Creating active pores plot using {len(dfs)} reads.")
    active_pores = Plot(
        path=path + "NanoComp_ActivePoresOverTime.html", title="Active pores over time"
    )
    data = []
    for sample, color in zip(df["dataset"].unique(), palette):
        pores = dfs.loc[dfs["dataset"] == sample, "channelIDs"].resample("10T").nunique()
        data.append(
            go.Scatter(
                x=pores.index.total_seconds() / 3600,
                y=pores,
                opacity=0.75,
                name=sample,
                marker=dict(color=color),
            )
        )

    active_pores.fig = go.Figure(
        {
            "data": data,
            "layout": go.Layout(
                barmode="overlay",
                title=title or active_pores.title,
                xaxis=dict(title="Time (hours)"),
                yaxis=dict(title="Active pores (per 10 minutes)"),
            ),
        }
    )

    active_pores.fig.update_layout(title_x=0.5)

    active_pores.html = active_pores.fig.to_html(full_html=False, include_plotlyjs="cdn")
    active_pores.save(settings)

    return active_pores
