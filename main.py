import os.path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from statsmodels.distributions.empirical_distribution import ECDF


filename = 'Protein lengths.xlsx'
bacteria = pd.read_excel(filename, sheet_name='Бактерии')
eukaryotes = pd.read_excel(filename, sheet_name='Эукариоты')

def processing(df, series):
    df = df.rename(columns={'Name': 'Name.0', 'Length': 'Length.0'})
    frames = dict()
    for column in df.columns.tolist():
        if 'Name' not in column: continue
        name, n = column.split('.')
        split = df[column].str.split('_', expand=True)
        frame = pd.DataFrame({
            'Name': split[0],
            'Mnemonics': split[1],
            'Length': df[f'Length.{n}'],
        }).dropna(axis=0)
        index = frame.pop('Mnemonics').unique()
        assert len(index) == 1
        name = index[0]
        frame['Organism'] = series[name]
        frames[name] = frame.set_index(['Organism', 'Name'])
    return frames

bacteria_mnemonics = pd.read_csv('csv/bacteria.csv', index_col='Мнемоника', squeeze=True)
eukaryotes_mnemonics = pd.read_csv('csv/eukaryotes.csv', index_col='Мнемоника', squeeze=True)
bacteria_ = processing(bacteria, bacteria_mnemonics)
eukaryotes_ = processing(eukaryotes, eukaryotes_mnemonics)

dist = [bacteria_['ECOLI'].copy(), eukaryotes_['YEAST'].copy()]
log_dist = [np.log10(d.values) for d in dist]
legend = [d.index[0][0] for d in dist]

def figure():
    fig = plt.figure(figsize=(12, 9))
    return fig, fig.gca()


fig, ax = figure()
ax.hist(log_dist, bins=50, alpha=0.75)
ax.legend(legend, loc='upper left')
ax.set_title('a. Joint Histogram')
ax.set_xlabel('Protein Length, log10')
ax.set_ylabel('Probability Density')
fig.savefig('img/a.png')


fig, ax = figure()
ax.boxplot(log_dist, widths=0.4)
ax.set_xticklabels(legend)
ax.set_title('b. Box Plots')
ax.set_ylabel('Protein Length, log10')
fig.savefig('img/b.png')


ecdf = [ECDF(d.squeeze()) for d in log_dist]
fig, ax = figure()
ax.plot(ecdf[0].x, ecdf[0].y, alpha=0.75)
ax.plot(ecdf[1].x, ecdf[1].y, alpha=0.75)
ax.legend(legend, loc='upper left')
ax.set_title('c. Empirical Distribution')
ax.set_xlabel('Protein Length, log10')
fig.savefig('img/c.png')
