import pandas as pd
import plotly.graph_objects as go
from argparse import ArgumentParser
from plotly import subplots


def main():
    args = get_args()
    states = ['single_pore', 'saturated', 'multiple', 'other', 'zero', 'unavailable']
    colors = ['#4feb34', '#000000', '#f77d02', '#8c8c8b', '#326ba1', '#399bf7']
    showlegends = [True] + [False] * len(args.mux_scan_files)
    fig = subplots.make_subplots(
        rows=len(args.mux_scan_files),
        cols=1,
        shared_xaxes=True,
        specs=[[{}] for i in range(len(args.mux_scan_files))],
        print_grid=False,
        vertical_spacing=0.05)
    for index, (muxfile, showlegend) in enumerate(zip(args.mux_scan_files, showlegends), start=1):
        df = pd.read_csv(muxfile, usecols=['mux_scan_assessment', 'repeat']) \
            .groupby("repeat")['mux_scan_assessment'] \
            .value_counts() \
            .unstack()

        for trace in [go.Bar(name=s, x=df.index, y=df[s],
                             marker_color=c, legendgroup=s, showlegend=showlegend)
                      for s, c in zip(states, colors)]:
            fig.append_trace(trace=trace, row=index, col=1)
    fig.update_layout(barmode='stack')
    for index, name in enumerate(args.names or args.mux_scan_files, start=1):
        fig.update_yaxes(title_text=name, row=index, col=1, range=[0, 12000])
    with open(args.output, 'w') as output:
        output.write(fig.to_html(full_html=True, include_plotlyjs='cdn'))


def get_args():
    parser = ArgumentParser(description="compare mux_scan_data files")
    parser.add_argument("mux_scan_files",
                        help="mux_scan_data csv files, optionally compressed",
                        nargs='*')
    parser.add_argument('-o', '--output', help="output html file name", default='poor_state.html')
    parser.add_argument('-n', '--names', help="names of datasets in mux_scan_files",
                        default=None, nargs='*')
    return parser.parse_args()


if __name__ == '__main__':
    main()
