import csv
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--intraspecies", help="flag to determine to include intraspecies corr scores, set to anything")
parser.add_argument("-f", "--filename", help="filename of output, default is data.tsv")
parser.add_argument("-m", "--max", help="produce a list of the best matching corr score only, set to anything")
parser.add_argument("-q", "--input", help="input file name")
args = parser.parse_args()

filename = args.filename if args.filename is not None else "data.tsv"

with open(filename, 'w') as csvfile:
  filewriter = csv.writer(csvfile, delimiter='\t')
  filewriter.writerow(['SpeciesA', 'SpeciesB', 'Corr_value', 'Seq_similarity'])
  with open(args.input, 'r') as f:
    reader = csv.reader(f)

    # read file row by row
    numOfMtces = 0
    geneList = []
    for i, row in enumerate(reader):
      if row[0] is "": # starter of a new corr-matrix 
        if len(row) is 1 : continue # get rid of empty corr-matrices
        numOfMtces += 1
        geneList = row # i.e. ['', 'AT3G48930', 'AT3G50020', 'BCRA1', 'EGF']
        # print(geneList)
        continue  

      maxScore = -1.01
      maxGene = ""
      for gIdx, gene in enumerate(geneList):
        if gene is "": continue # skip placeholder
        if gene == row[0]: continue # skip diagonals
        if (gene[0:2] == row[0][0:2]) and (args.intraspecies is None):
          continue # quick and easy way to filter intraspecies, NOTE: change for your species
        
        if args.max is None: # i.e. print all CC scores
          filewriter.writerow([row[0], gene, row[gIdx]])
        else: # write 'best' expressologs
          if row[gIdx] == 'NA': continue
          print(row[0])
          if float(row[gIdx]) > float(maxScore):
            maxScore = row[gIdx]
            maxGene = gene
          print(gIdx, len(geneList) - 1)
          if gIdx is len(geneList) - 1: # Before the loop finishes, write the largest CC score for that gene
            print(gIdx, len(geneList) - 1)
            print(row[0])
            filewriter.writerow([row[0], maxGene, maxScore])

    print(numOfMtces, "number of matrices processed!")

print(args.intraspecies)