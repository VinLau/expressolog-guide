import csv
from subprocess import PIPE, run

with open('final.tsv', 'w') as csvfile:
  filewriter = csv.writer(csvfile, delimiter='\t')
  with open('data.tsv', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for i, row in enumerate(reader):
      if i is 0: continue # skip header line
      # below command checks if order is in place
      result = run(['grep', row[0]+'.*;.*'+row[1], "all.bpo"], stdout=PIPE, stderr=PIPE, universal_newlines=True)
      splitStr = result.stdout.split(';')
      if len(splitStr) >= 7:
        filewriter.writerow([row[0], row[1], row[2], splitStr[6]])
      else:
        resultRevGrep = run(['grep', row[1]+'.*;.*'+row[0], "all.bpo"], stdout=PIPE, stderr=PIPE, universal_newlines=True)
        splitStrRevGrep = resultRevGrep.stdout.split(';')
        if len(splitStrRevGrep) >= 7:
          filewriter.writerow([row[0], row[1], row[2], splitStrRevGrep[6]])
        else:
          filewriter.writerow([row[0], row[1], row[2], "NA"])