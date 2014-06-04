# benchmarks: usage
#

# Run the benchmarks.py script on each platform
# types and number of tests are set in the global variables 
#   at the beginning of the file
platform="platform_name"
python benchmarks.py 20 | tee benchmarks-${platform}.log

# Once all the benchmarks are recorded in their own log files, run the graph
# script (if needed, titles and parameters can be changed at lines 200-225)
python benchmarks_graphs.py benchmarks-platform1.log benchmarks-platform2.log

