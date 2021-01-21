from scipy.optimize import linprog
import numpy as np


class Farrell:

    def __init__(self, results):
        self._efficiency = []
        self._lambdas = []
        self._slacks = []
        for res in results:
            self._efficiency.append(res.fun)
            self._lambdas.append(res.x[:len(res.x) - 1])
            self._slacks.append(res.slack[:len(res.slack - 1)])

    def eff(self):
        return self._efficiency

    def slack(self):
        return self._slacks

    def lambdas(self):
        return self._lambdas

    def print_attr(self, array, name):
        result = ''
        for i in range(1, len(array[0]) + 1):
            result += f'{name}{i}\t\t'
        result += '\n'
        for sls in array:
            result += '\t' * 2
            for sl in sls:
                if name == 'l' and sum(sls) == 0:
                    continue
                result += (f'%.2f\t' % sl)
            result += '\n'
        return result

    def __str__(self):
        result = ''
        result += 'Efficiency:\n'
        for e in self._efficiency:
            result += '%.2f\t' % e
        result += '\nLambdas: '
        result += self.print_attr(self._lambdas, 'l')
        result += 'Slacks: '
        result += self.print_attr(self._slacks, 'sl')
        return result


def create_equations(X, Y, RTS='vrs'):
    """
    Create equations for input efficiency
    :param X: 2-D array,
    Inputs
    :param Y: 2-D array
    Outputs
    :param RTS: technology for dea
    :return:
        c : 1-D array
            The coefficients of the linear objective function to be minimized.
        A_ub : 2-D array, optional
            The inequality constraint matrix. Each row of ``A_ub`` specifies the
            coefficients of a linear inequality constraint on ``x``.
        b_ub : 1-D array, optional
            The inequality constraint vector. Each element represents an
            upper bound on the corresponding value of ``A_ub @ x``.
        A_eq : 2-D array, optional
            The equality constraint matrix. Each row of ``A_eq`` specifies the
            coefficients of a linear equality constraint on ``x``.
        b_eq : 1-D array, optional
            The equality constraint vector. Each element of ``A_eq @ x`` must equal
            the corresponding element of ``b_eq``.
    details:
        This function creates LPP for input orientated efficiency
    """
    k, m, n = len(X), len(X[0]), len(Y[0])
    c = [0] * k + [1]
    A_ub = []
    b_ub = [0] * (m + n)
    A_eq = None
    b_eq = None
    # Input equations El >= SUM(x l) -> SUM(x l) <= El, let l is the second dmu
    for i in range(0, m):
        a = [0] * (k + 1)
        for j in range(0, k):
            a[j] = X[j][i]
        A_ub.append(a)
    # Output equations y <= SUM(y l) -> -SUM(y l) <= -l
    for i in range(0, n):
        a = [0] * (k + 1)
        for j in range(0, k):
            a[j] = -Y[j][i]
        A_ub.append(a)
    if RTS =='vrs':
        A_eq = []
        b_eq = []
        A_eq.append([1] * k + [0])
        b_eq = [1]
    return (c, A_ub, b_ub, A_eq, b_eq)


def simplex(X, Y, c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, dmu_number=None):
    """
    :param X:2-D array
    Inputs
    :param Y: 2-D array
    Outputs
    :param c:
    :param A_ub:
    :param b_ub:
    :param A_eq:
    :param b_eq:
    :param dmu_number:
    :return: res : OptimizeResult
        A :class:`scipy.optimize.OptimizeResult` consisting of the fields:

        x : 1-D array
            The values of the decision variables that minimizes the
            objective function while satisfying the constraints.
        fun : float
            The optimal value of the objective function ``c @ x``.
        slack : 1-D array
            The (nominally positive) values of the slack variables,
            ``b_ub - A_ub @ x``.
        con : 1-D array
            The (nominally zero) residuals of the equality constraints,
            ``b_eq - A_eq @ x``.
        success : bool
            ``True`` when the algorithm succeeds in finding an optimal
            solution.
        status : int
            An integer representing the exit status of the algorithm.

            ``0`` : Optimization terminated successfully.

            ``1`` : Iteration limit reached.

            ``2`` : Problem appears to be infeasible.

            ``3`` : Problem appears to be unbounded.

            ``4`` : Numerical difficulties encountered.

        nit : int
            The total number of iterations perfrormed in all phases.
        message : str
            A string descriptor of the exit status of the algorithm.

    """
    k, m, n = len(X), len(X[0]), len(Y[0])
    for i in range(0, m):
        A_ub[i][k] = -X[dmu_number][i]
    for i in range(0, n):
        b_ub[i + m] = -Y[dmu_number][i]
    result = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, method='highs-ds')
    return result


def dea(X, Y, RTS='vrs'):
    """
    :param X: 2-D array,
    Inputs
    :param Y: 2-D array
    Outputs
    :param RTS: technology for dea
    :return: 1-D array
        array of Linear problem solutions for each dmu
    details:
        Input orientated efficiency
    """
    results = []
    c, A_ub, b_ub, A_eq, b_eq = create_equations(X, Y, RTS=RTS)
    for dmu_number in range(0, len(X)):
        results.append(simplex(X, Y, c, A_ub, b_ub, A_eq, b_eq, dmu_number))
    farell = Farrell(results)
    return farell
