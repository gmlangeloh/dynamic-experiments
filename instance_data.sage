import glob

load("benchmarks.sage")

instances = glob.glob('./instances/*.ideal')
number_of_vars = {}

for instance in instances:
    b = Benchmark(instance)
    n = b.ideal.ring().ngens()
    if n not in number_of_vars:
        number_of_vars[n] = 1
    else:
        number_of_vars[n] += 1

total = 0
for key, value in sorted(number_of_vars.iteritems()):
    total += value
    print key, value, total
