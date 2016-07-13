import sys
import math

args = sys.argv
args.pop(0) #ignore self

output_file = args.pop()
data_file = args.pop()

sample_int = int(args.pop())
start_time = int(args.pop())

writer = open(output_file, "w")

with open(data_file, "r") as f:
    count = 0
    avg = 0.0
    avgSq = 0.0

    for line in f:
        if (not line.startswith("#")):
            data = line.strip().split()
            if (data != []): #ignore any lines start with \n
                time = int(data[0])
                value = float(data[1])
                dt = time - start_time
                if (dt >= 0 and dt % sample_int == 0):
                    print str(time) + ' ' + str(value)
                    avg += value
                    avgSq += value * value
                    count += 1
    avgSq /= float(count)
    avg /= float(count)
    std = avgSq - avg * avg
    error = std / math.sqrt(count) 
    writer.write(str(avg) + " " + str(error))

writer.close()
    
        
    
