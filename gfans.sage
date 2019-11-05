'''
Computes some Gröbner Fans for my AAECC paper.
'''

load("benchmarks.sage")

#These were selected because the GS algorithm finishes in less than 1 hour
instances = [ "cyclicn4",
              "cyclicnh4",
              "katsuran4",
              "katsuranh4",
              "econ4",
              "econh4",
              "econ5" ]

def get_ideal(name):
  prefix = "./instances/"
  suffix = ".ideal"
  return Benchmark(prefix + name + suffix).ideal

def gfan_data(I, outfile):
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
    f.write(str(basis_len(G)) + " " + str(basis_deg(G)) + "\n")

for instance in instances:
  I = get_ideal(instance)

  out_prefix = "./results/"
  out_suffix = ".gf"
  filename = out_prefix + instance + out_suffix
  with open(filename, "w") as f:
    gfan_data(I, f)
