from nanoplotter.plot import Plot
from nanoplotter.timeplots import check_valid_time_and_sort, add_time_bins
from nanomath import get_N50
import logging
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import plotly
import plotly.graph_objs as go
import sys
import joypy


def violin_or_box_plot(df, y, figformat, path, y_name,
                       title=None, plot="violin", log=False, palette=None):
    """Create a violin or boxplot from the received DataFrame.

    The x-axis should be divided based on the 'dataset' column,
    the y-axis is specified in the arguments
    """
    comp = Plot(path=path + "NanoComp_" + y.replace(' ', '_') + '.' + figformat,
                title="Comparing {}".format(y_name.lower()))

    if plot == 'violin':
        logging.info("NanoComp: Creating violin plot for {}.".format(y))
        process_violin_and_box(ax=sns.violinplot(x="dataset",
                                                 y=y,
                                                 data=df,
                                                 inner=None,
                                                 cut=0,
                                                 palette=palette,
                                                 linewidth=0),
                               log=log,
                               plot_obj=comp,
                               title=title,
                               y_name=y_name,
                               figformat=figformat,
                               ymax=np.amax(df[y]))
    elif plot == 'box':
        logging.info("NanoComp: Creating box plot for {}.".format(y))
        process_violin_and_box(ax=sns.boxplot(x="dataset",
                                              y=y,
                                              data=df,
                                              palette=palette),
                               log=log,
                               plot_obj=comp,
                               title=title,
                               y_name=y_name,
                               figformat=figformat,
                               ymax=np.amax(df[y]))
    elif plot == 'ridge':
        logging.info("NanoComp: Creating ridges plot for {}.".format(y))
        comp.fig, axes = joypy.joyplot(df,
                                       by="dataset",
                                       column=y,
                                       title=title or comp.title,
                                       x_range=[-0.05, np.amax(df[y])])
        if log:
            xticks = [float(i.get_text()) for i in axes[-1].get_xticklabels()]
            axes[-1].set_xticklabels([10**i for i in xticks])
        axes[-1].set_xticklabels(axes[-1].get_xticklabels(), rotation=30, ha='center')
        comp.save(format=figformat)
    else:
        logging.error("Unknown comp plot type {}".format(plot))
        sys.exit("Unknown comp plot type {}".format(plot))
    plt.close("all")
    return [comp]


def process_violin_and_box(ax, log, plot_obj, title, y_name, figformat, ymax):
    if log:
        ticks = [10**i for i in range(10) if not 10**i > 10 * (10**ymax)]
        ax.set(yticks=np.log10(ticks),
               yticklabels=ticks)
    ax.set(title=title or plot_obj.title,
           ylabel=y_name)
    plt.xticks(rotation=30, ha='center')
    plot_obj.fig = ax.get_figure()
    plot_obj.save(format=figformat)


def output_barplot(df, figformat, path, title=None, palette=None):
    """Create barplots based on number of reads and total sum of nucleotides sequenced."""
    logging.info("NanoComp: Creating barplots for number of reads and total throughput.")
    read_count = Plot(path=path + "NanoComp_number_of_reads." + figformat,
                      title="Comparing number of reads")
    ax = sns.countplot(x="dataset",
                       data=df,
                       palette=palette)
    ax.set(ylabel='Number of reads',
           title=title or read_count.title)
    plt.xticks(rotation=30, ha='center')
    read_count.fig = ax.get_figure()
    read_count.save(format=figformat)
    plt.close("all")

    throughput_bases = Plot(path=path + "NanoComp_total_throughput." + figformat,
                            title="Comparing throughput in gigabases")
    if "aligned_lengths" in df:
        throughput = df.groupby('dataset')['aligned_lengths'].sum()
        ylabel = 'Total gigabase aligned'
    else:
        throughput = df.groupby('dataset')['lengths'].sum()
        ylabel = 'Total gigabase sequenced'
    ax = sns.barplot(x=list(throughput.index),
                     y=throughput / 1e9,
                     palette=palette,
                     order=df["dataset"].unique())
    ax.set(ylabel=ylabel,
           title=title or throughput_bases.title)
    plt.xticks(rotation=30, ha='center')
    throughput_bases.fig = ax.get_figure()
    throughput_bases.save(format=figformat)
    plt.close("all")
    return read_count, throughput_bases


def n50_barplot(df, figformat, path, title=None, palette=None):
    n50_bar = Plot(path=path + "NanoComp_N50." + figformat,
                   title="Comparing read length N50")
    if "aligned_lengths" in df:
        n50s = [get_N50(np.sort(df.loc[df["dataset"] == d, "aligned_lengths"]))
                for d in df["dataset"].unique()]
        ylabel = 'Total gigabase aligned'
    else:
        n50s = [get_N50(np.sort(df.loc[df["dataset"] == d, "lengths"]))
                for d in df["dataset"].unique()]
        ylabel = 'Sequenced read length N50'
    ax = sns.barplot(x=list(df["dataset"].unique()),
                     y=n50s,
                     palette=palette,
                     order=df["dataset"].unique())
    ax.set(ylabel=ylabel,
           title=title or n50_bar.title)
    plt.xticks(rotation=30, ha='center')
    n50_bar.fig = ax.get_figure()
    n50_bar.save(format=figformat)
    plt.close("all")
    return [n50_bar]


def compare_sequencing_speed(df, figformat, path, title=None, palette=None):
    logging.info("NanoComp: creating comparison of sequencing speed over time.")
    seq_speed = Plot(path=path + "NanoComp_sequencing_speed_over_time." + figformat,
                     title="Sequencing speed over time")
    dfs = check_valid_time_and_sort(df, "start_time")
    dfs['timebin'] = add_time_bins(dfs)
    dfs = dfs.loc[dfs["duration"] > 0]
    ax = sns.violinplot(x=dfs["timebin"],
                        y=dfs["lengths"] / dfs["duration"],
                        hue=dfs["dataset"],
                        inner=None,
                        cut=0,
                        linewidth=0)
    ax.set(xlabel='Interval (hours)',
           ylabel="Sequencing speed (nucleotides/second)")
    plt.xticks(rotation=45, ha='center', fontsize=8)
    seq_speed.fig = ax.get_figure()
    seq_speed.save(format=figformat)
    plt.close("all")
    return [seq_speed]


def compare_cumulative_yields(df, path, palette=None, title=None):
    if palette is None:
        palette = plotly.colors.DEFAULT_PLOTLY_COLORS * 5
    dfs = check_valid_time_and_sort(df, "start_time").set_index("start_time")

    logging.info("NanoComp: Creating cumulative yield plots using {} reads.".format(len(dfs)))
    cum_yield_gb = Plot(path=path + "NanoComp_CumulativeYieldPlot_Gigabases.html",
                        title="Cumulative yield")
    data = []
    annotations = []
    for sample, color in zip(df["dataset"].unique(), palette):
        cumsum = dfs.loc[dfs["dataset"] == sample, "lengths"].cumsum().resample('10T').max() / 1e9
        data.append(go.Scatter(x=cumsum.index.total_seconds() / 3600,
                               y=cumsum,
                               opacity=0.75,
                               name=sample,
                               marker=dict(color=color))
                    )
        annotations.append(dict(xref='paper',
                                x=0.99,
                                y=cumsum[-1],
                                xanchor='left',
                                yanchor='middle',
                                text='{}Gb'.format(round(cumsum[-1])),
                                showarrow=False)
                           )

    cum_yield_gb.html = plotly.offline.plot({
        "data": data,
        "layout": go.Layout(barmode='overlay',
                            title=title or cum_yield_gb.title,
                            xaxis=dict(title="Time (hours)"),
                            yaxis=dict(title="Yield (gigabase)"),
                            annotations=annotations
                            )},
        output_type="div",
        show_link=False)

    cum_yield_gb.fig = go.Figure({
        "data": data,
        "layout": go.Layout(barmode='overlay',
                            title=title or cum_yield_gb.title,
                            xaxis=dict(title="Time (hours)"),
                            yaxis=dict(title="Yield (gigabase)"),
                            annotations=annotations
                            )})
    cum_yield_gb.save()
    return [cum_yield_gb]


def overlay_histogram(df, path, palette=None):
    """
    Use plotly to create an overlay of length histograms
    Return html code, but also save as png

    Only has 10 colors, which get recycled up to 5 times.
    """
    if palette is None:
        palette = plotly.colors.DEFAULT_PLOTLY_COLORS * 5

    hist = Plot(path=path + "NanoComp_OverlayHistogram.html",
                title="Histogram of read lengths")
    hist.html, hist.fig = plot_overlay_histogram(df, palette, title=hist.title)
    hist.save()

    hist_norm = Plot(path=path + "NanoComp_OverlayHistogram_Normalized.html",
                     title="Normalized histogram of read lengths")
    hist_norm.html, hist_norm.fig = plot_overlay_histogram(
        df, palette, title=hist_norm.title, histnorm="probability density")
    hist_norm.save()

    log_hist = Plot(path=path + "NanoComp_OverlayLogHistogram.html",
                    title="Histogram of log transformed read lengths")
    log_hist.html, log_hist.fig = plot_log_histogram(df, palette, title=log_hist.title)
    log_hist.save()

    log_hist_norm = Plot(path=path + "NanoComp_OverlayLogHistogram_Normalized.html",
                         title="Normalized histogram of log transformed read lengths")
    log_hist_norm.html, log_hist_norm.fig = plot_log_histogram(
        df, palette, title=log_hist_norm.title, histnorm="probability density")
    log_hist_norm.save()

    return [hist, hist_norm, log_hist, log_hist_norm]


def plot_overlay_histogram(df, palette, title, histnorm=""):
    data = []
    for d, c in zip(df["dataset"].unique(), palette):
        numreads = len(df.loc[df["dataset"] == d])
        data.append(
            go.Histogram(
                x=df.loc[df["dataset"] == d, "lengths"].sample(min(numreads, 10000)),
                opacity=0.4,
                name=d,
                histnorm=histnorm,
                marker=dict(color=c))
        )

    html = plotly.offline.plot(
        {"data": data,
         "layout": go.Layout(barmode='overlay',
                             title=title)},
        output_type="div",
        show_link=False)
    fig = go.Figure(
        {"data": data,
         "layout": go.Layout(barmode='overlay',
                             title=title)})
    return html, fig


def plot_log_histogram(df, palette, title, histnorm=""):
    """
    Plot overlaying histograms with log transformation of length
    Return both html and fig for png
    """
    data = []
    for d, c in zip(df["dataset"].unique(), palette):
        numreads = len(df.loc[df["dataset"] == d])
        data.append(
            go.Histogram(
                x=np.log10(df.loc[df["dataset"] == d, "lengths"].sample(min(numreads, 10000))),
                opacity=0.4,
                name=d,
                histnorm=histnorm,
                marker=dict(color=c))
        )
    xtickvals = [10**i for i in range(10) if not 10**i > 10 * np.amax(df["lengths"])]
    html = plotly.offline.plot(
        {"data": data,
         "layout": go.Layout(barmode='overlay',
                             title=title,
                             xaxis=dict(tickvals=np.log10(xtickvals),
                                        ticktext=xtickvals))},
        output_type="div",
        show_link=False)
    fig = go.Figure(
        {"data": data,
         "layout": go.Layout(barmode='overlay',
                             title=title,
                             xaxis=dict(tickvals=np.log10(xtickvals),
                                        ticktext=xtickvals))})
    return html, fig


def active_pores_over_time(df, path, palette=None, title=None):
    if palette is None:
        palette = plotly.colors.DEFAULT_PLOTLY_COLORS * 5
    dfs = check_valid_time_and_sort(df, "start_time").set_index("start_time")

    logging.info("NanoComp: Creating active pores plot using {} reads.".format(len(dfs)))
    active_pores = Plot(path=path + "NanoComp_ActivePoresOverTime.html",
                        title="Active pores over time")
    data = []
    for sample, color in zip(df["dataset"].unique(), palette):
        pores = dfs.loc[dfs["dataset"] == sample, "channelIDs"].resample('10T').nunique()
        data.append(go.Scatter(x=pores.index.total_seconds() / 3600,
                               y=pores,
                               opacity=0.75,
                               name=sample,
                               marker=dict(color=color))
                    )

    active_pores.html = plotly.offline.plot({
        "data": data,
        "layout": go.Layout(barmode='overlay',
                            title=title or active_pores.title,
                            xaxis=dict(title="Time (hours)"),
                            yaxis=dict(title="Active pores (per 10 minutes)"),
                            )},
        output_type="div",
        show_link=False)

    active_pores.fig = go.Figure({
        "data": data,
        "layout": go.Layout(barmode='overlay',
                            title=title or active_pores.title,
                            xaxis=dict(title="Time (hours)"),
                            yaxis=dict(title="Active pores (per 10 minutes)"),
                            )})
    active_pores.save()
    return active_pores
