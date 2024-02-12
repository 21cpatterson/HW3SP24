import math


def t_distribution_probability(degrees_of_freedom, z_value):
    """
    Compute the right-hand side of the t-distribution equation.

    Parameters:
    - degrees_of_freedom: int, degrees of freedom.
    - z_value: float, z value.

    Returns:
    - probability: float, probability from the t-distribution.
    """
    numerator = math.gamma((degrees_of_freedom + 1) / 2)
    denominator = math.sqrt(degrees_of_freedom * math.pi) * math.gamma(degrees_of_freedom / 2)
    coefficient = math.pow(1 + (math.pow(z_value, 2) / degrees_of_freedom), -((degrees_of_freedom + 1) / 2))

    probability = (numerator / denominator) * coefficient
    return probability


def print_table_comparison(degrees_of_freedom, z_values):
    """
    Print the computed probabilities and compare them with Table 25.2.

    Parameters:
    - degrees_of_freedom: int, degrees of freedom.
    - z_values: list of float, z values.
    """
    print(f"\nComparison for {degrees_of_freedom} degrees of freedom:")
    print("z       | Computed Probability | Table 25.2 Probability")
    print("--------|-----------------------|-----------------------")

    for z_value in z_values:
        computed_probability = t_distribution_probability(degrees_of_freedom, z_value)
        table_probability = table_25_2_probability(degrees_of_freedom, z_value)
        print(f"{z_value:.2f}    | {computed_probability:.6f}               | {table_probability:.6f}")


def table_25_2_probability(degrees_of_freedom, z_value):
    """
    Retrieve the probability from Table 25.2.

    Parameters:
    - degrees_of_freedom: int, degrees of freedom.
    - z_value: float, z value.

    Returns:
    - probability: float, probability from Table 25.2.
    """
    # Example values for Table 25.2, replace with actual values from the table
    table_values = {
        (7, 1.0): 0.775,
        (7, 1.5): 0.825,
        (7, 2.0): 0.860,
        (11, 1.0): 0.775,
        (11, 1.5): 0.825,
        (11, 2.0): 0.860,
        (15, 1.0): 0.775,
        (15, 1.5): 0.825,
        (15, 2.0): 0.860,
    }

    return table_values.get((degrees_of_freedom, z_value), 0.0)


def main():
    degrees_of_freedom = int(input("Enter degrees of freedom: "))
    z_values = [float(input(f"Enter z value {i + 1}: ")) for i in range(3)]

    for z_value in z_values:
        computed_probability = t_distribution_probability(degrees_of_freedom, z_value)
        print(
            f"\nComputed probability for {degrees_of_freedom} degrees of freedom and z = {z_value:.2f}: {computed_probability:.6f}")

    print_table_comparison(degrees_of_freedom, z_values)


if __name__ == "__main__":
    main()
