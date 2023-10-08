import re
from sympy import Matrix, lcm


def input_eq():
    """Prompts the user for the equation"""
    eq = input("Input the chemical equation: ")
    return eq


def format_elem_num(elem_num):
    """Helper function for add_elem. Returns the chemical element and its coefficient: Ca2 -> Ca, 2"""
    for i, char in enumerate(elem_num):
        if char.isdigit():
            return elem_num[:i], int(elem_num[i:])
    return elem_num, 1


def handle_parentheses_group(compound_elements):
    new_compound_elements = []
    for elem_num in compound_elements:
        if elem_num[0] == '(':
            temp = re.split(r'[()]', elem_num)
            group_elem, mult = temp[1], temp[2]
            elems = re.findall(r'[A-Z][a-z]*\d*|\([A-Za-z0-9]*\)[0-9]*', group_elem)
            for element in elems:
                elem, num = format_elem_num(element)
                num *= mult
                new_compound_elements.append(str())
        else:
            new_compound_elements.append(elem_num)


def add_elem(matrix_elements, side, sign, jump_ind):
    """Isolates the chemical elements from the compound and adds the coefficient to a dictionary"""
    for index_compound, compound in enumerate(side):
        compound_elements = re.findall(r'[A-Z][a-z]*\d*|\([A-Za-z0-9]*\)[0-9]*', compound)
        elem, num = None, None
        for elem_num in compound_elements:
            if elem_num[0] == '(':
                temp = re.split(r'[()]', elem_num)
                group_elem, mult = temp[1], int(temp[2])
                elems = re.findall(r'[A-Z][a-z]*\d*|\([A-Za-z0-9]*\)[0-9]*', group_elem)
                for element in elems:
                    elem, num = format_elem_num(element)
                    num *= mult
                    if elem not in matrix_elements:
                        matrix_elements[elem] = []
                    matrix_elements[elem].append((num * sign, index_compound + jump_ind))
            else:
                elem, num = format_elem_num(elem_num)
                if elem not in matrix_elements:
                    matrix_elements[elem] = []
                matrix_elements[elem].append((num * sign, index_compound + jump_ind))


def split_input(eq):
    """Splits the input into compounds on the left side of = and on right side"""
    left, right = re.split(r'\s*=\s*', eq)
    left = re.split(r'\s*\+\s*', left)
    right = re.split(r'\s*\+\s*', right)
    return left, right


def make_matrix(matrix_elements, rows, cols):
    """Completes a matrix with the dictionary data, representing a system of linear equations"""
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    row = 0
    for value in matrix_elements.values():
        for num, col in value:
            matrix[row][col] = num
        row += 1
        if row > len(matrix):
            break
    return matrix


def find_solution(matrix):
    """Returns the solution of nullspace matrix"""
    matrix = Matrix(matrix)
    sol = matrix.nullspace()[0]
    mult = lcm([val.q for val in sol])
    sol = sol * mult
    sol_list = sol.tolist()
    return sol_list


def build_output(solution, left, right):
    """Formats the output from the solution"""
    output = ""
    for i, compound in enumerate(left):
        if solution[i][0] != 1:
            output += str(solution[i][0])
        output += compound
        if i != len(left)-1:
            output += " + "
    output += " = "
    for i, compound in enumerate(right):
        if solution[i + len(left)][0] != 1:
            output += str(solution[i + len(left)][0])
        output += compound
        if i != len(right)-1:
            output += " + "
    return output


def solve_eq(eq):
    """Balances a chemical equation"""
    left, right = split_input(eq)
    matrix_elements = {}
    add_elem(matrix_elements, left, 1, 0)
    add_elem(matrix_elements, right, -1, len(left))
    matrix = make_matrix(matrix_elements, len(matrix_elements), len(left) + len(right))
    solution = find_solution(matrix)
    output = build_output(solution, left, right)
    return output


eq = "LiOH + Fe3(PO4)2 = Fe(OH)2 + Li3PO4"
# eq = input_eq()
print(solve_eq(eq))
