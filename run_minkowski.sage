import sys

load("minkowski.sage")

if len(sys.argv) == 2 and sys.argv[1] == "min_hilbert":
    run_all_min()
else:
    run_all()
