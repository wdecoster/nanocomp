import numpy as np
from argparse import ArgumentParser
import plotly.graph_objects as go
from os.path import basename
from pandas import read_feather


def main():
    args = get_args()
    readlengths = [read_feather(p)["lengths"] for p in args.input]
    names = args.names if args.names else [basename(f).replace('.feather', '') for f in args.input]
    yield_bars = plot_yields(readlengths, names)
    readplot, logreadplot = plot_read_length(readlengths, names, maxlength=args.maxlength)
    with open(args.output, 'w') as output:
        output.write(yield_bars)
        output.write(readplot)
        output.write(logreadplot)


def get_args():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", help="input feather files", nargs='+')
    parser.add_argument("-n", "--names", help="names of datasets in input", nargs='+')
    parser.add_argument("-m", "--maxlength", help="maximum read length to show", default=100000)
    parser.add_argument("-o", "--output", help="output html file", default='read_lengths.html')
    return parser.parse_args()


def plot_read_length(lengths, names, maxlength=100000):
    sampled_lengths = [np.random.choice(a, size=10000) for a in lengths]

    fig = go.Figure()
    bins = round(maxlength / 500)
    for dataset, name in zip(sampled_lengths, names):
        counts, bins = np.histogram([l for l in dataset if not l > maxlength],
                                    bins=bins,
                                    density=True)
        fig.add_trace(go.Bar(x=bins[1:],
                             y=counts,
                             opacity=0.4,
                             name=name,
                             text=bins[1:],
                             hoverinfo="text",
                             hovertemplate=None))
    fig.update_layout(barmode='overlay',
                      title_text='Aligned read length histogram',
                      hovermode="x unified",
                      bargap=0)
    fig["layout"]["xaxis"].update(title_text="Read length")
    length = fig.to_html(full_html=False, include_plotlyjs='cdn')

    fig = go.Figure()
    bins = round(maxlength / 500)
    for dataset, name in zip(sampled_lengths, names):
        counts, bins = np.histogram(np.log10([l for l in dataset if not l > maxlength]),
                                    bins=bins,
                                    density=True)
        fig.add_trace(go.Bar(x=bins[1:],
                             y=counts,
                             opacity=0.4,
                             name=name,
                             text=[10**i for i in bins[1:]],
                             hoverinfo="text",
                             hovertemplate=None))
    fig.update_layout(barmode='overlay',
                      title_text='Aligned log-transformed read length histogram',
                      hovermode="x unified",
                      bargap=0)
    maxl = max(np.amax(l) for l in sampled_lengths)
    xtickvals = [10**i for i in range(10) if not 10**i > 10 * maxl]
    fig.update_layout(xaxis=dict(tickvals=np.log10(xtickvals), ticktext=xtickvals))
    fig["layout"]["xaxis"].update(title_text="Read length")
    loglength = fig.to_html(full_html=False, include_plotlyjs='cdn')

    return length, loglength


def plot_yields(lengths, names):
    yields = [np.sum(l)/1e9 for l in lengths]
    for n, y in zip(names, yields):
        print(f"{n}: {y}")
    fig = go.Figure(go.Bar(x=names, y=yields))
    fig["layout"]["yaxis"].update(title_text="Yield (Gb)")
    fig.update_layout(title_text='Yield per individual')
    return fig.to_html(full_html=True, include_plotlyjs='cdn')


if __name__ == '__main__':
    main()
