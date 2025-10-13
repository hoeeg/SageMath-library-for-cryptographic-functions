def differential_uniformity(poly, F):
    """
    Compute the differential uniformity of a polynomial over a finite field
    
    :param poly: Polynomial function
    :param F: Finite field GF(2^n)
    :return: Differential uniformity of the polynomial
    """
    poly_values = {x: poly(x) for x in F}

    max_count = 0
    
    for delta_i in F: 
        if delta_i == 0:
            continue

        for delta_o in F:
            count = 0
            for x in F:
                if poly_values[x + delta_i] - poly_values[x] == delta_o:
                    count += 1

            max_count = max(max_count, count)

    return max_count


# Example usage instructions:
# from mycrypto import differential_uniformity
# F.<a> = GF(2^8)                                               # Define the finite field GF(2^8) with generator 'a'
# R.<y> = PolynomialRing(F) OR R.<x> = PolynomialRing(F)        # Define the polynomial ring over the finite field
# poly = ...                                                    # Define your polynomial here
# differential_uniformity(poly, F)

