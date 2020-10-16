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


def plot_read_length(lengths, names, maxlength=2000):
    sampled_lengths = [np.random.choice(a, size=10000) for a in lengths]

    fig = go.Figure()
    for dataset, name in zip(sampled_lengths, names):
        fig.add_trace(go.Histogram(x=[l for l in dataset if not l > maxlength],
                                   opacity=0.4,
                                   name=name,
                                   histnorm="probability density"))
    fig.update_layout(barmode='overlay', title_text='Aligned read length histogram')
    fig["layout"]["xaxis"].update(title_text="Read length")
    lenght = fig.to_html(full_html=False, include_plotlyjs='cdn')

    fig = go.Figure()
    for dataset, name in zip(sampled_lengths, names):
        fig.add_trace(go.Histogram(x=np.log10([l for l in dataset if not l > maxlength]),
                                   opacity=0.4,
                                   name=name,
                                   histnorm="probability density"))
    maxl = max(np.amax(l) for l in sampled_lengths)
    xtickvals = [10**i for i in range(10) if not 10**i > 10 * maxl]
    fig.update_layout(xaxis=dict(tickvals=np.log10(xtickvals), ticktext=xtickvals))
    fig.update_layout(barmode='overlay', title_text='Aligned log-transformed read length histogram')
    loglength = fig.to_html(full_html=False, include_plotlyjs='cdn')
    return lenght, loglength


def plot_yields(lengths, names):
    yields = [np.sum(l)/1e9 for l in lengths]
    for n, y in zip(names, yields):
        print(f"{n}: {y}")
    fig = go.Figure(go.Bar(x=names, y=yields))
    fig["layout"]["yaxis"].update(title_text="Yield (Gb)")
    return fig.to_html(full_html=True, include_plotlyjs='cdn')


if __name__ == '__main__':
    main()
