import sys

def read_polynomial(line, variables, R):
    words = line.split()
    polynomial = R(0)
    monomial = R(1)
    last_variable = ""
    is_exponent = False
    for word in words:
        if word in variables:
            monomial *= R(word)
            last_variable = word
            is_exponent = False
        elif word == "+":
            polynomial += monomial
            monomial = R(1)
            is_exponent = False
        elif word == "-":
            polynomial += monomial
            monomial = R(-1)
            is_exponent = False
        elif word == "^":
            is_exponent = True
        elif "^" in word:
            parts = word.split("^")
            variable = parts[0]
            exponent = parts[1]
            monomial *= R(variable) ** int(exponent)
        else:
            try:
                coef = int(word)
                if is_exponent:
                    monomial *= last_variable ** (coef - 1)
                else:
                    monomial *= coef
            except:
                pass
    polynomial += monomial
    return polynomial

if len(sys.argv) < 2:
    raise Exception("Not enough arguments, requires filename.")

filename = sys.argv[1]

with open(filename, "r") as f:
    lines = f.readlines()
    algorithm = lines[0]
    characteristic = int(lines[1])
    num_vars = int(lines[2])
    dummy = lines[3]
    variables = lines[4].split()
    num_polys = int(lines[5])
    polynomials = []
    R = PolynomialRing(GF(characteristic), names = variables)
    for i in range(6, 6 + num_polys):
        p = read_polynomial(lines[i], variables, R)
        polynomials.append(p)
    print(polynomials)
