try:
    import stats
except:
    pass

def chisquare(matrix):
    sum_all = float(sum([sum(x) for x in matrix]))
    def sum_col(i):
        return sum(matrix[i])
    def sum_row(i):
        return sum([x[i] for x in matrix])
    dim_x = len(matrix)
    dim_y = len(matrix[0])

    matrix_expected = []
    for i in range(0, dim_x):
        col = []
        for j in range(0, dim_y):
            element = sum_col(i)*sum_row(j)/sum_all
            col.append(element)
        matrix_expected.append(col)

    matrix_chi = []
    for i in range(0, dim_x):
        col = []
        for j in range(0, dim_y):
            element = matrix[i][j]-matrix_expected[i][j]
            element *= element
            divide_by = matrix_expected[i][j] if matrix_expected[i][j]!=0 else 1
            element /= divide_by
            col.append(element)
        matrix_chi.append(col)
    chi = sum([sum(x) for x in matrix_chi])
    return chi, stats.chisqprob(chi, (dim_x-1)*(dim_y-1))
