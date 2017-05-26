import os
import time
import re
import sys

total = 0.0
N = 100 # number of cases (seeds) to test

r = re.compile("\d+.\d+")

for seed in range(1,N):

    vis_command = "java AbstractWarsVis -exec \"/home/tadek/Coding/Topcoder/AbstractWars/./" + sys.argv[1] + "\" -seed "
    vis_command = vis_command + str(seed) + " -novis"

    start_time = time.time()
    output = os.popen(vis_command).readlines()
    finish_time = time.time()
    time_elapsed = finish_time - start_time

    #total = total + int(output[-1])

    f = re.search(r, output[-1]).group()
    total = total + float(f)

    #print(output)
    print("Case " + str(seed) + " time: " + str(time_elapsed) + " score: " + output[-1], end="")


print("Mean: " + str(total/N))
