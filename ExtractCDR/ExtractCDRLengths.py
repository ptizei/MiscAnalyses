import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re

seqList = []
idList = []
CDR3List = []
lenCDR3List = []
i = 0


MotifUpstream = 'AVYYC'
MotifDownstream = 'WGQGT'
FileToOpen = 'TestSeqs.fasta'

with open(FileToOpen, 'r') as seqsFile:
    Lines = seqsFile.readlines()
    for line in Lines:
        if i % 2 == 0:
            idList.append(line[1:].rstrip('\n'))
        else:
            seqList.append(line.rstrip('\n'))
        i+=1

ZipIdSeqs = zip(idList,seqList)
dfSeqs = pd.DataFrame(ZipIdSeqs, columns = ['id', 'seq'])


for seq in dfSeqs['seq']:
    CDR3 = re.search(f'{MotifUpstream}(.+?){MotifDownstream}', seq)
    CDR3List.append(CDR3.group(1))
    lenCDR3List.append(len(CDR3.group(1)))

dfSeqs['CDR3Seq'] = CDR3List
dfSeqs['CDR3Length'] = lenCDR3List
NumBins = max(dfSeqs.CDR3Length) - min(dfSeqs.CDR3Length) +1
plot = sns.distplot(dfSeqs.CDR3Length, bins = NumBins)
fig = plot.get_figure()
fig.savefig('CDR3LengthHistogram.png')

dfSeqs.to_csv('CDR3LengthTable.csv')
