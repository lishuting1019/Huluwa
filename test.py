#分式方程——去分母
import os
import sys
sys.path.append(os.path.join(os.path.expanduser("~"), "git/latex2sympy"))
from process_latex import process_sympy
from sympy import *
from sympy import ratsimp
from sympy.parsing.sympy_parser import parse_expr
def ReductionFrac(expr,*vars):
    expr_1 = expr.args[0]
    expr_2 = expr.args[1]
    n = together(expr_1,vars[0])
    n = parse_expr(str(n), evaluate=False)
    #得到分母m
    if n.func.__name__ == 'Pow' and n.args[1] == -1:
            m = n.args[0]
    if n.func.__name__ == 'Mul':
        for sub_n in n.args:
            if sub_n.func.__name__ == 'Pow' and sub_n.args[1] == -1:
                m = sub_n.args[0]
    expr_lh = Mul(n,m,evaluate=False)#左边乘分母m，不进行化简
    expr_rh = Mul(expr_2,m,evaluate=False)#右边乘分母m，不进行化简
    expr_new = expr.func(expr_lh,expr_rh)
    desc = "分式方程去分母:"
    result = [{'desc': desc,'expr': expr_new}]
    return result