import numpy as np
np.set_printoptions(precision=7, suppress=True, linewidth=100)


def nevilles_method(x_points, y_points, x):
    # creating the matrix
    size = len(x_points)
    matrix = np.zeros((size, size))
    # populating the matrix
    for counter, row in enumerate(matrix):
        row[0] = y_points[counter]

    num_of_points = size


    # for loop to evaluate coefficient of equation
    for i in range(1, num_of_points):
        for j in range(1, i + 1):
            first_multiplication = (x - x_points[i - j]) * matrix[i][j - 1]
            second_multiplication = (x - x_points[i]) * matrix[i - 1][j - 1]
            denominator = x_points[i] - x_points[i - j]

        coefficient = (first_multiplication - second_multiplication) / denominator
        matrix[i][j] = coefficient

    return matrix[i - 1][j - 1]

# function to find diagonal values of the matrix in Newton's forward method
def divided_difference_table(x_points_2, y_points_2):
    # set up the matrix
    size: int = len(x_points_2)
    matrix: np.array = np.zeros((size, size))

    # fill the matrix
    for index, row in enumerate(matrix):
        row[0] = y_points_2[index]

    # populate the matrix (end points are based on matrix size and max operations
    for i in range(1, size):
        for j in range(1, i + 1):
            # the numerator is the immediate left and diagonal left indices...
            numerator = matrix[i][j - 1] - matrix[i - 1][j - 1]
            denominator = x_points_2[i] - x_points_2[i - j]
            operation = numerator / denominator

            # cut it off to view it simpler
            matrix[i][j] = '{0:.7g}'.format(operation)

    # prints only the diagonal values of the matrix
    print("[", matrix[1][1], ",", matrix[2][2], ",", matrix[3][3], "]")
    return matrix

# function for approximating the result of Newton's forward method
def get_approximate_result(matrix, x_points_2, value):
    # p0 is always y0 and, we use a reoccuring x to avoid having to recalculate x
    reoccuring_x_span = 1
    reoccuring_px_result = matrix[0][0]

    # we only need the diagonals...and that starts at the first row...
    size = len(matrix)
    for index in range(1, size):
        polynomial_coefficient = matrix[index][index]
        # we use the previous index for x_points....
        reoccuring_x_span *= (value - x_points_2[index - 1])

        # get a_of_x * the x_span
        mult_operation = polynomial_coefficient * reoccuring_x_span
        # add the reoccuring px result
        reoccuring_px_result += mult_operation

    # final result
    return reoccuring_px_result



def apply_div_dif(matrix: np.array):
    size = len(matrix)
    for i in range(2, size):
        for j in range(2, i + 2):
            # skipping already filled values of the matrix
            if j >= len(matrix[i]) or matrix[i][j] != 0:
                continue

            # something get left and diagonal left
            left: float = matrix[i][j - 1]
            diagonal_left: float = matrix[i - 1][j - 1]
            numerator: float = left - diagonal_left
            denominator = matrix[i][0] - matrix[i - j + 1][0]

            # something save into matrix
            operation = numerator / denominator
            matrix[i][j] = operation

    print(matrix)

# function for hermite interpolating a matrix of values
def hermite_interpolation():
    # defining x, y, and slopes of the matrix
    x_points = [3.6, 3.8, 3.9]
    y_points = [1.675, 1.436, 1.318]
    slopes = [-1.195, -1.188, -1.182]

    # creating a matrix
    num_of_points = len(x_points)
    matrix = np.zeros((num_of_points * 2, num_of_points * 2))

    # populate x values
    index = 0
    for x in range(0, len(matrix), 2):
        matrix[x][0] = x_points[index]
        matrix[x + 1][0] = x_points[index]
        size = len(matrix)


        if index >= size:
            index = 0

        index += 1

    # prepopulate
    # prepopulate y values
    index = 0
    for x in range(0, len(matrix), 2):
        matrix[x][1] = y_points[index]
        matrix[x + 1][1] = y_points[index]
        index += 1

    # prepopulate with derivatives (every other row)
    index = 0
    for x in range(1, len(matrix), 2):
        matrix[x][2] = slopes[index]
        index += 1

    apply_div_dif(matrix)

# function for determining the matrix A of cubic spline interpolating
def cubic_spline_matrixA(x_points_3):
    # creating a matrix
    size = len(x_points_3)
    matrix: np.array = np.zeros((size, size))

    # values to be input into matrix
    h0 = x_points_3[1] - x_points_3[0]
    h1 = x_points_3[2] - x_points_3[1]
    h2 = x_points_3[3] - x_points_3[2]

    # populating the matrix
    for index, row in enumerate(matrix):
        row[0] = [1, h0, 0, 0][index]
        row[1] = [0, 2 * (h0 + h1), h1, 0][index]
        row[2] = [0, h1, 2 * (h1 + h2), 0][index]
        row[3] = [0, 0, h2, 1][index]

    print(matrix)
    return matrix

# finding vector b for cubic spline interpolating
def cubic_spline_vectorB(x_points_3, y_points_3):
    # values for completing calculation
    h0 = x_points_3[1] - x_points_3[0]
    h1 = x_points_3[2] - x_points_3[1]
    h2 = x_points_3[3] - x_points_3[2]

    a1 = y_points_3[0]
    a2 = y_points_3[1]
    a3 = y_points_3[2]
    a4 = y_points_3[3]

    # values calculated for vector
    value_1 = 0
    value_2 = ((3 / h1) * (a3 - a2)) - ((3 / h0) * (a2 - a1))
    value_3 = ((3 / h2) * (a4 - a3)) - ((3 / h1) * (a3 - a2))
    value_4 = 0

    # creating the vector
    list1 = [value_1, value_2, value_3, value_4]

    vector_b = np.array(list1)

    print(vector_b)
    return vector_b


def main():
    # Question 1 Neville's Method
    x_points = [3.6, 3.8, 3.9]
    y_points = [1.675, 1.436, 1.318]
    x = 3.7
    print(nevilles_method(x_points, y_points, x))
    print("")


    # Question 2 Newton's Forward Method
    x_points_2 = [7.2, 7.4, 7.5, 7.6]
    y_points_2 = [23.5492, 25.3913, 26.8224, 27.4589]
    divided_table = divided_difference_table(x_points_2, y_points_2)
    print("")


    # Question 3 Newton's Forward Method
    approximating_x = 7.3
    final_approximation = get_approximate_result(divided_table, x_points_2, approximating_x)
    print(final_approximation)
    print("")


    # Question 4 Divided Difference Method
    hermite_interpolation()
    print("")

    # Question 5 Cubic Spline Interpolation
    x_points_3 = [2, 5, 8, 10]
    y_points_3 = [3, 5, 7, 9]
    a = cubic_spline_matrixA(x_points_3)
    print("")
    b = cubic_spline_vectorB(x_points_3, y_points_3)
    print("")
    # using numpy library and linear algebra equation to calculate x vector values
    print(np.dot(np.linalg.inv(a), b))

if __name__ == "__main__":
    main()
