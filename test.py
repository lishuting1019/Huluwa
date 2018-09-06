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

#分式方程——去分母的check
def ReductionFrac_check(expr):
    lh = expr.args[0]
    rh = expr.args[1]
    #分式方程右边是常数，且只有一项
    if rh.is_real and len(rh.args) < 2:
        return True
    else:
        return False

#分式方程——合并同类项
def TogetherFrac(expr):
    """
    此处合并同类项局限于经过主线加法合并
    :param expr:
        Sympy表达式,可能是:
            多项式:
                合并同类项
            关系式:
                左右两边同时合并
    :param x:
        指定变量
    :return:
         result:
            desc:
                描述选项,为None表示经过本函数没有任何变化
            expr:
                新的表达式
    """

    desc = "分式方程合并同类项:"
    # 判断是否为关系式(等式、不等式)
    if expr.is_Relational:
        expr_new = expr.func(apart(simplify(expr.args[0])), simplify(expr.args[1]))
    else:
        expr_new = simplify(expr)
    return [{'desc': desc, 'expr': expr_new}]

#分式方程——合并同类项的check
def TogetherFrac_check(expr):
    """
    # 多项式
            加法化简
            能合并同类项
    # 关系式
            左边、右边分别按照多项式处理
    :param expr:
    :param x:
    :return:
    """
    need_simplify = False
    if expr.is_Relational:
        for sub_expr in expr.args:
            if not is_sum_simplify(sub_expr):
                need_simplify = True
                break
    else:
        if not is_sum_simplify(expr):
            need_simplify = True
    return need_simplify

#分式方程——移项
def TranspositionFrac(expr):
    # 加载需要的Package
    import os
    import sys
    sys.path.append(os.path.join(os.path.expanduser("~"), "git/latex2sympy"))
    from process_latex import process_sympy
    import sympy
    from sympy.parsing.sympy_parser import parse_expr
    expr_new = expr
    desc = None # desc == None 表示经过该函数,没有任何变化
    try:
        # 判断是否为关系式(等式、不等式)
        if expr.is_Relational:
            lh = expr.args[0]
            rh = expr.args[1]
            if rh.is_Add:
                if len(rh.args) == 0:
                    item = rh
                else:
                    item = rh.args[-1]
            else:
                item = rh
            # 相反数,也许可以优化，此处用字符串处理判断
            if expr.is_Equality:
                desc = "移项:等式两边同时"
            else:
                desc = "移项:不等式两边同时"
            if str(item).startswith('-'):
                item = parse_expr(str(item)[1:], evaluate=False)
                desc = ''.join([desc, "加上", str(item)])
            elif str(item).startswith('(-'):
                item = str(item)
                item = ''.join(['(+', item[2:]])
                item = parse_expr(item, evaluate=False)
                desc = ''.join([desc, "加上", str(item)])
            else:
                desc = ''.join([desc, "减去", str(item)])
                item = sympy.Mul(-1, item, evaluate=False)
            expr_new = expr.func(sympy.Add(lh, item, evaluate=False),
                                 sympy.Add(rh, item, evaluate=False))
    except:
        pass
    result = [{'expr': expr_new, 'desc': desc}]
    return result


#分式方程——移项的check
def TranspositionFrac_check(expr):
    """
        关系式右边不为0
    """
    if expr.is_Relational and str(expr.args[1]) != '0':
        return True
    else:
        return False

#分式方程——常数项移项
#形如a/x+b=c的一元一次方程（其中c可以为任何实数，a、b不为0）的常数项移项
import os
import sys
sys.path.append(os.path.join(os.path.expanduser("~"), "git/latex2sympy"))
from process_latex import process_sympy
from sympy import *
import sympy
from sympy.parsing.sympy_parser import parse_expr
def TransferFrac(expr):
    lh = expr.args[0]
    rh = expr.args[1]
    lis1 = 0
    lis2 = 0
    #把常数项移到等式右边，方便两边同乘分母
    for sub in lh.args:
        if sub.is_real:
            lis1 = lis1 + sub
        else:
            lis2 = lis2 + sub
    lh_new = lis2
    rh_new = -lis1
    expr_new = Eq(lh_new,rh_new)
    desc = "分式方程的常数项移项:"
    result = [{'desc': desc,'expr': expr_new}]
    return(result)


#分式方程——常数项移项的check
#左边有常数，右边为0
def TransferFrac_check(expr):
    lh = expr.args[0]
    rh = expr.args[1]
    for sub in lh.args:
        if sub.is_real and rh ==0:
            return True
    else:
        return False


#分式方程——通分
def tongfen(expr):
    lh = expr.args[0]
    rh = expr.args[1]
    lh_1 = radsimp(lh)#左边通分，方便提取分母
    expr_new = expr.func(lh_1,rh)
    desc = "分式方程通分:"
    result = [{'desc': desc,'expr': expr_new}]
    return(result)


#分式方程——通分的check
def tongfen_check(expr):
    lh = expr.args[0]
    rh = expr.args[1]
    #左边能够通分到最简式
    if lh.is_Add and radsimp(lh)!=lh:
        return True
    else:
        return False


import os
import sys
sys.path.append(os.path.join(os.path.expanduser("~"), "git/latex2sympy"))
from process_latex import process_sympy
from sympy import *
def SimplifyFrac(expr):
    expr_1 = expr.args[0]
    expr_2 = expr.args[1]
    expr_lh = expand(expr_1,evaluate = False)#左边去括号，不进行化简
    expr_rh = expand(expr_2,evaluate = False)#右边去括号，不进行化简
    expr_new = expr.func(expr_lh,expr_rh)
    desc = "分式方程去括号:"
    result = [{'desc': desc,'expr': expr_new}]
    return result



#分式方程——去括号的check
def SimplifyFrac_check(expr):
    lh = expr.args[0]
    rh = expr.args[1]
    #可以进行去括号
    if expand(lh) != lh or expand(rh) != rh:
        return True
    else:
        return False


import os
import sys
sys.path.append(os.path.join(os.path.expanduser("~"), "git/latex2sympy"))
from process_latex import process_sympy
from sympy import *
from sympy import ratsimp
def simplifyfraction(expr):
    result  = []
    expr_1 = expr.args[0]
    expr_2 = expr.args[1]
    #去分母后进行化简，保留括号形式
    lh = simplify(expr_1)
    rh = simplify(expr_2)
    expr_new1 = Eq(lh,rh)
    result.append({'desc': '分式方程化简', 'expr': expr_new1,'step':'分式方程化简'})
    return result


#分式方程——化简的check
def simplifyfraction_check(expr):
    lh = expr.args[0]
    rh = expr.args[1]
    #可以进行化简
    if simplify(lh) != lh or simplify(rh) != rh:
        return True
    else:
        return False


# 方程组 - 加法消元法(第一步)
def AddEquations(expr, *vars):
    lh = expr[0]
    rh = expr[1]
    lh_1 = lh.args[0]
    lh_2 = lh.args[1]
    rh_1 = rh.args[0]
    rh_2 = rh.args[1]

    # x的系数
    xl = Poly(lh_1, vars[0]).all_coeffs()
    xr = Poly(rh_1, vars[0]).all_coeffs()
    coeff_xl = xl[0]
    coeff_xr = xr[0]
    # x两个系数的最小公倍数，化为相同的系数
    lcmx = lcm(coeff_xl, coeff_xr)
    x1 = lcmx / coeff_xl
    x2 = lcmx / coeff_xr
    if x1 < 0:
        new_xl1 = -x1 * lh_1
        new_xl2 = -x1 * lh_2
        new_xr1 = x2 * rh_1
        new_xr2 = x2 * rh_2
    elif x2 < 0:
        new_xl1 = x1 * lh_1
        new_xl2 = x1 * lh_2
        new_xr1 = -x2 * rh_1
        new_xr2 = -x2 * rh_2
    else:
        new_xl1 = x1 * lh_1
        new_xl2 = x1 * lh_2
        new_xr1 = x2 * rh_1
        new_xr2 = x2 * rh_2

    # y的系数
    yl = Poly(lh_1, vars[1]).all_coeffs()
    yr = Poly(rh_1, vars[1]).all_coeffs()
    coeff_yl = yl[0]
    coeff_yr = yr[0]
    # y两个系数的最小公倍数，化为相同的系数
    lcmy = lcm(coeff_yl, coeff_yr)
    y1 = lcmy / coeff_yl
    y2 = lcmy / coeff_yr
    if y1 < 0:
        new_yl1 = -y1 * lh_1
        new_yl2 = -y1 * lh_2
        new_yr1 = y2 * rh_1
        new_yr2 = y2 * rh_2
    elif y2 < 0:
        new_yl1 = y1 * lh_1
        new_yl2 = y1 * lh_2
        new_yr1 = -y2 * rh_1
        new_yr2 = -y2 * rh_2
    else:
        new_yl1 = y1 * lh_1
        new_yl2 = y1 * lh_2
        new_yr1 = y2 * rh_1
        new_yr2 = y2 * rh_2
    # x的系数小于y的系数，对x进行加法消元
    if abs(lcmx) <= abs(lcmy) and lcmx < 0:
        result = []
        with evaluate(False):
            expr_x1 = [Eq(new_xl1, new_xl2), Eq(new_xr1, new_xr2)]
        result.append({'desc': '化为相同系数', 'expr': expr_x1, 'step': '化为相同系数', 'con_exprs': [Eq(new_xl1, new_xl2)]})
        with evaluate(False):
            expr_x2 = Eq(new_xl1 + new_xr1, new_xl2 + new_xr2)
        result.append({'desc': '方程相加', 'expr': expr_x2, 'step': '方程相加', 'con_exprs': [Eq(new_xl1, new_xl2)]})
        return result
    # y的系数小于x的系数，对y进行加法消元
    elif abs(lcmy) <= abs(lcmx) and lcmy < 0:
        result = []
        with evaluate(False):
            expr_y1 = [Eq(new_yl1, new_yl2), Eq(new_yr1, new_yr2)]
        result.append({'desc': '化为相同系数', 'expr': expr_y1, 'step': '化为相同系数', 'con_exprs': [Eq(new_yl1, new_yl2)]})
        with evaluate(False):
            expr_y2 = Eq(new_yl1 + new_yr1, new_yl2 + new_yr2)
        result.append({'desc': '方程相加', 'expr': expr_y2, 'step': '方程相加', 'con_exprs': [Eq(new_yl1, new_yl2)]})
        return result




#方程组 - 加法消元法的check函数
def AddEquations_check(expr,*vars):
    #x(y)的两个系数一正一负，比较x，y的最小公倍数，对小的进行化简
    lh = expr[0]
    rh = expr[1]
    lh_new = lh.args[0]
    rh_new = rh.args[0]
    #x的系数
    xl= Poly(lh_new,vars[0]).all_coeffs()
    xr = Poly(rh_new,vars[0]).all_coeffs()
    coeff_xl = xl[0]
    coeff_xr = xr[0]
    lcmx = lcm(coeff_xl,coeff_xr)
    #y的系数
    yl = Poly(lh_new,vars[1]).all_coeffs()
    yr = Poly(rh_new,vars[1]).all_coeffs()
    coeff_yl = yl[0]
    coeff_yr = yr[0]
    lcmy = lcm(coeff_yl,coeff_yr)
    t = (lcmx,lcmy)
    if abs(lcmx) <= abs(lcmy) and lcmx<0:
        return True
    if abs(lcmy) <= abs(lcmx) and lcmy<0:
        return True
    else:
        return False



#加法消元法的回代
#把第一步得到的结果带入第一个方程里
def back3(expr,**supports):
    con_exprs = supports['con_exprs']#储存第一个方程
    result = []
    expr1 = expr
    arr = con_exprs[0].subs(expr.args[0],expr.args[1])
    with evaluate(False):
        expr_new = arr
    result.append({'desc': '回代:'+str(con_exprs[0].args[0])+'='+str(con_exprs[0].args[1])+'得', 'expr': expr_new,'step':'回代得','con_exprs': expr1})
    return result


# 方程组 - 加法消元法(回代)的check函数
def back3_check(expr):
    return True


#加减法得到的两个解x,y，合并输出
def back4(expr,**supports):
    con_exprs = supports['con_exprs']#存的第一个解
    result = []
    expr_new3 = [expr,con_exprs]#expr是第二个解
    result.append({'desc': '最终解', 'expr':expr_new3,'step':'最终解'})
    return result


def back4_check(expr):
    lh = expr.args[0]
    if lh.is_Add:
        return False
    else:
        return True


# 方程组 - 减法消元法(第一步)
def SubEquations(expr, *vars):
    lh = expr[0]
    rh = expr[1]
    lh_1 = lh.args[0]
    lh_2 = lh.args[1]
    rh_1 = rh.args[0]
    rh_2 = rh.args[1]

    # x的系数
    xl = Poly(lh_1, vars[0]).all_coeffs()
    xr = Poly(rh_1, vars[0]).all_coeffs()
    coeff_xl = xl[0]
    coeff_xr = xr[0]
    # x两个系数的最小公倍数，化为相同的系数
    lcmx = lcm(coeff_xl, coeff_xr)
    x1 = lcmx / coeff_xl
    x2 = lcmx / coeff_xr
    if x1 > 0:
        new_xl1 = x1 * lh_1
        new_xl2 = x1 * lh_2
        new_xr1 = x2 * rh_1
        new_xr2 = x2 * rh_2
    elif x2 > 0:
        new_xl1 = x1 * lh_1
        new_xl2 = x1 * lh_2
        new_xr1 = x2 * rh_1
        new_xr2 = x2 * rh_2
    else:
        new_xl1 = x1 * lh_1
        new_xl2 = x1 * lh_2
        new_xr1 = x2 * rh_1
        new_xr2 = x2 * rh_2

    # y的系数
    yl = Poly(lh_1, vars[1]).all_coeffs()
    yr = Poly(rh_1, vars[1]).all_coeffs()
    coeff_yl = yl[0]
    coeff_yr = yr[0]
    # y两个系数的最小公倍数，化为相同的系数
    lcmy = lcm(coeff_yl, coeff_yr)
    y1 = lcmy / coeff_yl
    y2 = lcmy / coeff_yr
    if y1 > 0:
        new_yl1 = y1 * lh_1
        new_yl2 = y1 * lh_2
        new_yr1 = y2 * rh_1
        new_yr2 = y2 * rh_2
    elif y2 > 0:
        new_yl1 = y1 * lh_1
        new_yl2 = y1 * lh_2
        new_yr1 = y2 * rh_1
        new_yr2 = y2 * rh_2
    else:
        new_yl1 = y1 * lh_1
        new_yl2 = y1 * lh_2
        new_yr1 = y2 * rh_1
        new_yr2 = y2 * rh_2
    # x的系数小于y的系数，对x进行减法消元
    if abs(lcmx) <= abs(lcmy) and lcmx > 0:
        result = []
        with evaluate(False):
            expr_x1 = [Eq(new_xl1, new_xl2), Eq(new_xr1, new_xr2)]
        result.append({'desc': '化为相同系数', 'expr': expr_x1, 'step': '化为相同系数', 'con_exprs': [Eq(new_xl1, new_xl2)]})
        with evaluate(False):
            expr_x2 = Eq(new_xl1 - new_xr1, new_xl2 - new_xr2)
        result.append({'desc': '方程相减', 'expr': expr_x2, 'step': '方程相减', 'con_exprs': [Eq(new_xl1, new_xl2)]})
        return result
    # y的系数小于x的系数，对y进行减法消元
    elif abs(lcmy) < abs(lcmx) and lcmy > 0:
        result = []
        with evaluate(False):
            expr_y1 = [Eq(new_yl1, new_yl2), Eq(new_yr1, new_yr2)]
        result.append({'desc': '化为相同系数', 'expr': expr_y1, 'step': '化为相同系数', 'con_exprs': [Eq(new_yl1, new_yl2)]})
        with evaluate(False):
            expr_y2 = Eq(new_yl1 - new_yr1, new_yl2 - new_yr2)
        result.append({'desc': '方程相减', 'expr': expr_y2, 'step': '方程相减', 'con_exprs': [Eq(new_yl1, new_yl2)]})
        return result



#方程组 - 减法消元法的check函数
def SubEquation_check(expr,*vars):
    #x(y)的两个系数都是正的或负的，比较x，y的最小公倍数，对小的进行化简
    lh = expr[0]
    rh = expr[1]
    lh_new = lh.args[0]
    rh_new = rh.args[0]
    #x的系数
    xl= Poly(lh_new,vars[0]).all_coeffs()
    xr = Poly(rh_new,vars[0]).all_coeffs()
    coeff_xl = xl[0]
    coeff_xr = xr[0]
    lcmx = lcm(coeff_xl,coeff_xr)
    #y的系数
    yl = Poly(lh_new,vars[1]).all_coeffs()
    yr = Poly(rh_new,vars[1]).all_coeffs()
    coeff_yl = yl[0]
    coeff_yr = yr[0]
    lcmy = lcm(coeff_yl,coeff_yr)
    t = (lcmx,lcmy)
    if abs(lcmx)<= abs(lcmy) and lcmx>0:
        return True
    if abs(lcmy)<= abs(lcmx) and lcmy>0:
        return True
    else:
        return False



#方程组 - 减法消元法(回代)
#把第一步得到的结果带入第一个方程里
def back2(expr,**supports):
    con_exprs = supports['con_exprs']#储存第一个方程
    result = []
    expr1 = expr
    arr = con_exprs[0].subs(expr.args[0],expr.args[1])
    with evaluate(False):
        expr_new = arr
    result.append({'desc': '回代:'+str(con_exprs[0].args[0])+'='+str(con_exprs[0].args[1])+'得', 'expr': expr_new,'step':'回代得','con_exprs': expr1})
    return result


#一次函数-属性(斜率,截距,交点)(y=kx+b,k≠0)
def functionk(expr,*vars):
    expr1 = expand(expr)
    lh = expr1.args[0]
    rh = expr1.args[1]
    x1= Poly(rh,vars[0]).all_coeffs()
    #斜率
    expr_new = x1[0]
    expr2 = Eq(0,expr.args[1])
    expr3 = solve(expr2)
    #x轴交点
    expr_new2 = (expr3[0],0)
    #增减性
    expr_new4 = 'y随x的增大而增大'
    expr_new5 = 'y随x的增大而减小'
    #截距#y轴交点
    if rh.is_Add:
        results = solve(expr,vars[0])
        expr4 = Eq(vars[0],results[0])
        expr1 = Eq(0,expr4.args[1])
        expr2 = solve(expr1)
        expr_new1 = expr2[0]
        expr_new3 = (0,expr_new1)
    else:
        expr_new1 = 0
        expr_new3 = (0,expr_new1)
    result = []
    result.append({'desc': '函数斜率k=', 'expr': expr_new,'step':'函数斜率k='})
    with evaluate(False):
        result.append({'desc': '函数截距b=', 'expr': expr_new1,'step':'函数截距b='})
    with evaluate(False):
        result.append({'desc': '与x轴交点', 'expr': expr_new2,'step':'与x轴交点','value':expr_new2})
        result.append({'desc': '与y轴交点', 'expr': expr_new3,'step':'与y轴交点','value':expr_new3})
    with evaluate(False):
        if expr_new>0:
            result.append({'desc': '增减性', 'expr': expr_new4, 'step':'增减性'})
        else:
            result.append({'desc': '增减性', 'expr': expr_new5, 'step':'增减性'})
    return result




#二次函数的属性(y=ax**2+bx+c)
def symmetryaxis(expr,*vars):
    rh = expr.args[1]
    coeff = poly(rh,vars[0]).all_coeffs()
    expr2 = Eq(0,rh)
    expr_new2 = solve(expr2)
    expr_new3 = '向上'
    expr_new4 = '向下'
    if 0 not in coeff:
        a = coeff[0]
        b = coeff[1]
        c = coeff[2]
    elif coeff[1] ==0 and coeff.count(0) ==1:
        a = coeff[0]
        b = 0
        c = coeff[2]
    elif coeff[2] ==0 and coeff.count(0) ==1:
        a = coeff[0]
        b = coeff[1]
        c = 0
    elif coeff[1] ==0 and coeff[2] ==0:
        a = coeff[0]
        b = 0
        c = 0
    coff = b*b - 4*a*c
    expr_new = Eq(vars[0],-b/(2*a))
    x2 = Eq(vars[0],-b/(2*a))
    y2 = expr.subs(x2.args[0],x2.args[1])
    expr_new1 = (x2.args[1],y2.args[1])
    result = []
    result.append({'desc': '函数对称轴为直线','expr': expr_new,'step':'函数对称轴为直线'})
    #判断函数与x轴是否有交点
    if coff == 0:
        expr_new7 = (-b/(2*a),0)
        result.append({'desc': '零点为','expr': expr_new7,'step':'零点为','value':expr_new7})
    elif coff < 0:
        expr_new10 = "没有交点"
        result.append({'desc': '零点为','expr': expr_new10,'step':'零点为'})
    elif coff > 0:
        expr_new7 = ((-b+sqrt(coff))/(2*a),0)
        expr_new9 = ((-b-sqrt(coff))/(2*a),0)
        result.append({'desc': '零点1为','expr': expr_new7,'step':'零点1为','value':expr_new7})
        result.append({'desc': '零点2为','expr': expr_new9,'step':'零点2为','value':expr_new9})
    x3 = Eq(vars[0],0)
    y3 = expr.subs(x3.args[0],x3.args[1])
    expr_new8 = (0,y3.args[1])
    result.append({'desc': '与y轴交点', 'expr': expr_new8, 'step':'与y轴交点','value':expr_new8})
    result.append({'desc': '函数顶点坐标','expr': expr_new1,'step':'函数顶点坐标','value':expr_new1})
    if a>0:
        result.append({'desc': '函数最小值点为','expr': expr_new1,'step':'函数最小值点为'})
        result.append({'desc': '开口方向','expr': expr_new3,'step':'开口方向'})
    else:
        result.append({'desc': '函数最大值点为','expr': expr_new1,'step':'函数最大值点为'})
        result.append({'desc': '开口方向','expr': expr_new4,'step':'开口方向'})
    return result


# 反比例函数的属性(对称中心,增减性)
def inversefunction(expr, *vars):
    lh = expr.args[0]
    rh = expr.args[1]
    k = rh * vars[0]
    expr_new = (0, 0)
    expr_new1 = 'y=x'
    expr_new2 = 'y=-x'
    expr_new3 = 'y随x的增大而增大'
    expr_new4 = 'y随x的增大而减小'

    result = []
    result.append({'desc': '对称中心', 'expr': expr_new, 'step': '对称中心'})
    with evaluate(False):
        if k > 0:
            result.append({'desc': '对称轴为直线', 'expr': expr_new2, 'step': '对称中心'})
            result.append({'desc': '增减性', 'expr': expr_new4, 'step': '增减性'})

        else:
            result.append({'desc': '对称轴为直线', 'expr': expr_new1, 'step': '对称中心'})
            result.append({'desc': '增减性', 'expr': expr_new3, 'step': '增减性'})
    return result     