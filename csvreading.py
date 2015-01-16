import csv
import StringIO
filename = open('RungeKuttaOutput.csv','rb')
reader = csv.reader(filename,delimiter=',')
x1list = [row[1] for row in reader]
x1list.pop(0)
x1listout=[float(x1list[i]) for i in range(len(x1list))]
