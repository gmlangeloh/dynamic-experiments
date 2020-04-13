'''
Computes some Gröbner Fans for my AAECC paper.
'''

load("benchmarks.sage")
load("dynamicgb.pyx") #I can use the bitlength evaluation tools in stats.pxi

#These were selected because the GS algorithm finishes in less than 1 hour
instances = [ "cyclicn4",
              "cyclicnh4",
              "cyclicn5",
              "katsuran4",
              "katsuranh4",
              "econ4",
              "econh4",
              "econ5" ]

def get_ideal(name, char0):
  prefix = "./instances/" if not char0 else "./instances-char0/"
  suffix = ".ideal"
  return Benchmark(prefix + name + suffix).ideal

def gfan_data(I, outfile, char0):
  gf = I.groebner_fan()
  bases = gf.reduced_groebner_bases()

  def basis_len(G):
    return len(G)

  def basis_deg(G):
    return max([g.total_degree() for g in G])

  #minlen = min([ basis_len(G) for G in bases ])
  #mindeg = min([ basis_deg(G) for G in bases ])

  f.write("polynomials degree\n")
  for G in bases:
    f.write(str(basis_len(G)) + " " + str(basis_deg(G)))
    if char0:
      tot, avg, mx = bitlength_stats(list(G))
      f.write(str(tot) + " " + str(avg) + " " + str(mx) + "\n")
    else:
      f.write("\n")

char0 = False
if len(sys.argv) > 1:
  if sys.argv[1] == "0":
    char0 = True

for instance in instances:
  I = get_ideal(instance, char0)

  out_prefix = "./results/"
  out_suffix = ".gf"
  filename = out_prefix + instance + out_suffix
  with open(filename, "w") as f:
    gfan_data(I, f, char0)
