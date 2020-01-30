import pandas as pd
import plotly.graph_objs as go
import plotly
from argparse import ArgumentParser


def main():
    args = get_args()

    df = pd.concat([pd.read_csv(f, sep="\t")["end_reason"].value_counts()
                    for f in args.summaries],
                   axis='columns')
    df.columns = args.names
    with open("End_reason_comparison.html", 'w') as output:
        output.write(plot(trace=[go.Bar(name=t, x=args.names, y=df.loc[t])
                                 for t in df.index],
                          layout=dict(barmode='stack',
                                      title="End reason per dataset [absolute]",
                                      yaxis_title="Number of reads")))
        output.write(plot(trace=[go.Bar(name=t, x=args.names, y=df.loc[t] / df.sum())
                                 for t in df.index],
                          layout=dict(barmode='stack',
                                      title="End reason per dataset [relative]",
                                      yaxis_title="Fraction of reads")))


def get_args():
    parser = ArgumentParser()
    parser.add_argument("-s", "--summaries",
                        help="Unbasecalled summary files to compare", nargs='+', required=True)
    parser.add_argument("-n", "--names", help="Names of dataset in --summaries",
                        nargs='+', required=True)
    return parser.parse_args()


def plot(trace, layout):
    return plotly.offline.plot(dict(data=trace, layout=layout),
                               output_type="div",
                               show_link=False,
                               include_plotlyjs='cdn')


if __name__ == '__main__':
    main()
