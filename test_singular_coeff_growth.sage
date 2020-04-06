import time
load("benchmarks.sage")
load("dynamicgb.pyx")
c7h = Benchmark("./instances-char0/cyclicnh7.ideal")
t = time.time()
gb = list(c7h.ideal.groebner_basis())
t = time.time() - t
print(t)
print(bitlength_stats(gb))

gb = dynamic_gb(c7h.ideal.gens(), algorithm="caboara_perry", print_results=True, print_coefficients=True)
gb = dynamic_gb(c7h.ideal.gens(), algorithm="perturbation", print_results=True, print_coefficients=True)
gb = dynamic_gb(c7h.ideal.gens(), algorithm="random", print_results=True, print_coefficients=True)
