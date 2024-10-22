#Converts the coefficients created in polyfit to a copy/pasteable polynomail

def polynomial_to_string(coefficients):
    terms = []
    degree = len(coefficients) - 1
    for i, coeff in enumerate(coefficients):
        if coeff != 0:
            power = degree - i
            if power == 0:
                terms.append(f"{coeff}")
            elif power == 1:
                terms.append(f"{coeff} * x")
            else:
                terms.append(f"{coeff} * x**{power}")
    return " + ".join(terms)