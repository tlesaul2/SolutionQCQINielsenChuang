#!/usr/bin/env python3

import itertools
import functools
import operator
from collections import Counter, defaultdict
import pickle
import subprocess





class Permutation(list):

    def __init__(self, L):
        if set(range(len(L))) != set(L):
            raise ValueError("A Permutation can only be initialized with an iterable of distinct integers from 0 to it's intended length")
        super().__init__(L)


    def __mul__(self, other):
        return Permutation([other[entry] for entry in self])

    def __hash__(self):
        return hash(tuple(self))


# Toffolis indexed by the targeted qubit index
Toffolis = [Permutation([0,1,2,7,4,5,6,3]), \
            Permutation([0,1,2,3,4,7,6,5]), \
            Permutation([0,1,2,3,4,5,7,6])]

# CNOTs indexed by
CNOTs = [Permutation([0,1,3,2,4,5,7,6]), # CNOT(x2,x3)\
         Permutation([0,1,2,3,5,4,7,6]), # CNOT(x1,x3)\
         Permutation([0,1,2,3,6,7,4,5]), # CNOT(x1,x2)\
         Permutation([0,3,2,1,4,7,6,5]), # CNOT(x3,x2)\
         Permutation([0,5,2,7,4,1,6,3]), # CNOT(x3,x1)\
         Permutation([0,1,6,7,4,5,2,3])] # CNOT(x2, x1)

NOTs = [Permutation([4,5,6,7,0,1,2,3]), # NOT(x1)\
        Permutation([2,3,0,1,6,7,4,5]), # NOT(x2)\
        Permutation([1,0,3,2,5,4,7,6])] # NOT(x3)

qpic_lines = {0: b"+a1 a2 a3",
              1: b"a1 +a2 a3",
              2: b"a1 a2 +a3",
              3: b"a2 +a3",
              4: b"a1 +a3",
              5: b"a1 +a2",
              6: b"+a2 a3",
              7: b"+a1 a3",
              8: b"+a1 a2",
              9: b"+a1",
              10:b"+a2",
              11:b"+a3",
             }

def generate_qpic(indices):
    preamble = b"a1 W\na2 W\na3 W\n"
    return preamble+b"\n".join(map(qpic_lines.__getitem__, indices))

def tikz(qpic):
    p = subprocess.Popen(["qpic"], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    stdout,_ = p.communicate(qpic)
    return str(stdout, "ascii")

identity = Permutation(range(8))

Gates = Toffolis+CNOTs+NOTs

#total = 0
#tests = 0
#for i in itertools.product(range(len(Gates)),repeat=2):
    #tests += 1
    #if Gates[i[0]]*Gates[i[1]] in itertools.chain([identity,],Gates):
        #total+=1
        #print(f"{Gates[i[0]]}x{Gates[i[1]]}={Gates[i[0]]*Gates[i[1]]}")

#print(f"Out of {tests} tests, {total} products of 2 gates are the same as single gates, or the identity")


target = Permutation([0,2,3,4,5,6,7,1])

permutations = set(map(lambda x: Permutation([0,]+list(x)), itertools.permutations(range(1,8))))
#permutations = set(map(Permutation, itertools.permutations(range(8))))

gate_counts = Counter()
length_counts = defaultdict(int)

previous = set([identity,])
current = set()
permutations.remove(identity)
length_counts[0] = 1
gate_counts[identity] = 0
length = 0
indices = defaultdict(list)
indices[identity] = [ [], ]
representations = defaultdict(lambda: defaultdict(int))
representations[0][identity] = 1
print(f"Found {length_counts[length]} permutation that can be expressed by {length} Toffolis and/or CNOTs")
while permutations:
    for gate in Gates:
        for perm in previous:
            permutation = perm*gate
            if permutation in permutations:
                gate_counts[permutation] = length+1
                length_counts[length+1] += 1
                permutations.remove(permutation)
            if gate_counts[permutation] == length+1:
                current.add(permutation)
                representations[length+1][permutation] += representations[length][perm]
                for index_list in indices[perm]:
                    indices[permutation].append(index_list+[Gates.index(gate),])
                #if Gates.index(gate) >= 9:
                    #print(f"Found a representation of {permutation} using a NOT. indices: {indices[permutation][-1]}")

    print(f"Found {length_counts[length+1]} permutations that can be expressed by {length+1} Toffolis and/or CNOTs")

    length+=1
    previous=current
    current=set()

min_length = min(length for length in range(max(representations.keys())) if target in representations[length])
print(f"The target permutation {target} can be performed in {min_length} CNOT and Toffoli gates")

exit(0)

for permutation in permutations:
    if len(set(map(len, indices[permutation])))!=1:
        raise Exception(f"There are index lists for {permutation} with different lengths")

def breakdown(index_list):
    return  len(list(filter(lambda x: x >=9, index_list))), len(list(filter(lambda x: 3 <= x < 9, index_list))), len(list(filter(lambda x: x < 3, index_list)))

breakdown_counter = defaultdict(lambda: defaultdict(Counter))
indices_by_breakdown = defaultdict(lambda : defaultdict(list))
target = Permutation([0,1,3,2,5,4,6,7])

if len(gate_counts) != 7*6*5*4*3*2*1:
    raise Exception(f"There aren't 8! permutations in the gate_counts dictionary")


for permutation,length in gate_counts.items():
    for index_list in indices[permutation]:
        break_down = breakdown(index_list)
        breakdown_counter[length][permutation][break_down] += 1
        indices_by_breakdown[permutation][break_down].append(index_list)

length = 3
for permutation in breakdown_counter[length].keys():
    print(f"{permutation} has {representations[permutation][length]} minimal representation of length {length}")
    print(f"\tBy gate qubit-count:",end=" ")
    for break_down, count in sorted(breakdown_counter[length][permutation].items(),reverse=True):

        print(f"{break_down}:{count}")
        if len(set(map(frozenset, indices_by_breakdown[permutation][break_down])))!=1:
            print(f"\tATTENTION: there are distinct collections of gates below")
        for index_list in indices_by_breakdown[permutation][break_down]:
            print(f"\t\t{index_list}")

exit(0)

length = 1


while True:
    found = False
    for i in itertools.product(range(len(Gates)),repeat=length):
        if functools.reduce(operator.__mul__, map(Gates.__getitem__, i))==target:
            print(f"{'x'.join(map(str,list(map(Gates.__getitem__, i))))}={target}")
            print(i)
            break
    else:
        break
    length+=1

permutations = list(map(lambda x: Permutation([0,]+list(x)), itertools.permutations(range(1,8))))
#permutations = list(map(Permutation, itertools.permutations(range(8))))

with open("eight_gate_permutations.pkl", "rb") as F:
    eight_gate_permutations = pickle.load(F)

four_gate_permutations = set()
indices = defaultdict(list)
representations = defaultdict(int)

for i in itertools.product(range(len(Gates)),repeat=4):
    permutation = functools.reduce(operator.__mul__, map(Gates.__getitem__, i))
    four_gate_permutations.add(permutation)
    representations[permutation] += 1
    indices[permutation].append(tuple(i))

for perm1, perm2 in itertools.product(four_gate_permutations, repeat=2):
    permutation = perm1*perm2
    if permutation in eight_gate_permutations:
        representations[permutation] += 1
        for indices1, indices2 in zip(indices[perm1], indices[perm2]):
            indices[permutation].append(tuple(list(indices1)+list(indices2)))

for i,permutation in enumerate(eight_gate_permutations):
    print(f"For eight_gate_permutation #{i}: {permutation}, counted {representations[permutation]} representations as products of two four_gate_permutations and collected {len(indices[permutation])} distinct index lists")

def Bmatrix(permutation):
    def column(elem):
        val = [0,]*8
        val[elem] = 1
        return val
    transpose = [column(elem) for elem in permutation]
    matrix = [ [transpose[j][i] for j in range(8)] for i in range(8)]
    table = r"""$\begin{bmatrix}
"""
    table += """
""".join([' & '.join(list(map(str, row)))+r" \\" for row in matrix])
    table+= """
\end{bmatrix}$"""
    return table

for i,permutation in enumerate(eight_gate_permutations):
    with open(f"eight_gate_permutation.{i}.tex","w") as F:
        print(r"""\documentclass[11pt]{book}
\usepackage[utf8]{inputenc}
\usepackage{amsmath,amsfonts,amssymb,amsthm}
\usepackage[letterpaper,left=0.5in,top=0.65in,bottom=0.65in]{geometry}
\usepackage[dvipdfmx]{graphicx}
\usepackage[english]{babel}

%図の場所をなるべく指定した場所にする
\usepackage{booktabs}
\usepackage{here}

\RequirePackage[l2tabu, orthodox]{nag}
\usepackage[all, warning]{onlyamsmath}

%単位を書くときに使う
\usepackage{siunitx}

\usepackage{CJKutf8}
\usepackage{ascmac} % screen
\usepackage{ulem}
\usepackage{cases}
\usepackage{braket}
\usepackage{dsfont}
\usepackage{ascmac}
\usepackage{url}
\PassOptionsToPackage{hyphens}{url}\usepackage{hyperref} % hyper link
\usepackage{ccicons} % creative commons license icon
\usepackage{fancyhdr} % footer
\usepackage{blkarray}
%\pagestyle{fancy}
%\cfoot[\href{http://creativecommons.org/licenses/by-nc-sa/4.0/}{Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License}.]{}
\usepackage{color}
\usepackage{calc}
\usepackage{stmaryrd}
\usepackage{adjustbox}
\usepackage{lscape}
\usepackage{tocloft}
\usepackage{comment}
\usepackage{changepage}

\usepackage{tikz}
\usepackage{circuitikz}
\usetikzlibrary{decorations.pathreplacing,decorations.pathmorphing}

%\usepackage{xcolor}
%\definecolor{Pagecolor}{RGB}{35,38,41} % Dark gray
%\definecolor{Textcolor}{RGB}{255,255,255} % Light gray
%\pagecolor{Pagecolor}
%\color{Textcolor}

\usepackage[T1]{fontenc}
\usepackage{lmodern}

%\newcommand{\Gate}[1]{{\fontfamily{pcr}\selectfont#1}}
\newcommand{\NAND}{{\fontfamily{pcr}\selectfont NAND}}
\newcommand{\AND}{{\fontfamily{pcr}\selectfont AND}}
\newcommand{\OR}{{\fontfamily{pcr}\selectfont OR}}
\newcommand{\NOR}{{\fontfamily{pcr}\selectfont NOR}}
\newcommand{\NOT}{{\fontfamily{pcr}\selectfont NOT}}
\newcommand{\XOR}{{\fontfamily{pcr}\selectfont XOR}}
\newcommand{\FANOUT}{{\fontfamily{pcr}\selectfont FANOUT}}
\newcommand{\CROSSOVER}{{\fontfamily{pcr}\selectfont CROSSOVER}}
\newcommand{\CNOT}{{\fontfamily{pcr}\selectfont CNOT}}


\usepackage{mathtools}
\DeclarePairedDelimiter\Abs{\lvert}{\rvert}
\DeclarePairedDelimiter\Ceil{\lceil}{\rceil}
\DeclarePairedDelimiter\Floor{\lfloor}{\rfloor}
\newcommand{\abs}[1]{\Abs*{#1}}
\newcommand{\ceil}[1]{\Ceil*{#1}}
\newcommand{\floor}[1]{\Floor*{#1}}
%\newcommand{\abs}[1]{\big|#1\big|}
%\newcommand{\ceil}[1]{\left\lceil#1\right\rceil}
%\newcommand{\floor}[1]{\left\lfloor#1\right\rfloor}

%\renewcommand{\familydefault}{\sfdefault}
%% Only use the math font of mathpazo
%\let\temp\rmdefault
%\usepackage{fourier}
%\let\rmdefault\temp

\usepackage{fancyhdr}
\setlength{\headheight}{15.2pt}
\pagestyle{fancy}
\lhead[\leftmark ]{\thepage}
\rhead[\thepage]{\leftmark}

\cfoot{\footnotesize \textcopyright 2018 goropikari, \textcopyright 2021 tlesaul2 - \href{http://creativecommons.org/licenses/by-nc-sa/4.0/}{Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License}}

% コマンド定義
\DeclareMathOperator{\Tr}{tr}
\newcommand{\norm}[1]{\left\lVert#1\right\rVert} % norm ||x||
\newcommand{\kb}[1]{\ket{#1}\hspace{-1mm} \bra{#1}} % |x><x|
\newcommand{\kbt}[2]{\ket{#1}\hspace{-1mm} \bra{#2}} % |x><y|
\newcommand{\Textbf}[1]{\hspace{3mm}\\\textbf{#1)}}
\newcommand{\Soln}{\textbf{\\Soln: }}
\newcommand{\Mod}[1]{\ (\mathrm{mod}\ #1)}
\newtheorem{thm}{Theorem.}[section]
\newtheorem{prop}{Proposition.}[section]



\title{Select Solutions for ``Quantum Computation and Quantum Information: 10th Anniversary Edition" by Nielsen and Chuang}
\author{Original author: goropikari\\Extended by: tlesaul2}
\date{\today}

\let\tmp\oddsidemargin
\let\oddsidemargin\evensidemargin
\let\evensidemargin\tmp
\reversemarginpar

\let\Re\undefined
\let\Im\undefined
\DeclareMathOperator{\Re}{\mathfrak{Re}}
\DeclareMathOperator{\Im}{\mathfrak{Im}}

\begin{document}
\maketitle
\thispagestyle{empty}
\cleardoublepage
\thispagestyle{empty}
\setcounter{page}{1} % 表紙のページを0ページにする

\section*{Copylight Notice:}
\ccbyncsa\\
    This work is licensed under a \href{http://creativecommons.org/licenses/by-nc-sa/4.0/}{Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License} by the original author.  As such, the second author does the same.


\section*{Repository}
As of November, 2021, the original source \LaTeX \ code, located at \url{https://github.com/goropikari/SolutionForQuantumComputationAndQuantumInformation}\ has not been updated since April 2020.  The extended source \LaTeX \ code is located at \\ \url{https://github.com/tlesaul2/SolutionQCQINielsenChuang}\ .  It may be updated more actively.

\section*{For readers}
This is an unofficial solution manual for "\href{http://www.cambridge.org/jp/academic/subjects/physics/quantum-physics-quantum-information-and-quantum-computation/quantum-computation-and-quantum-information-10th-anniversary-edition?format=HB&isbn=9781107002173#BBFv83H3ofgcgG3A.97}{Quantum Computation and Quantum Information: 10th Anniversary Edition}" (ISBN-13: 978-1107002173) by Michael A. Nielsen and Isaac L. Chuang.\\

\noindent From the original author:\\
\indent I have studied quantum information theory as a hobby.
And I'm not a researcher.
So there is no guarantee that these solutions are correct.
Especially because I'm not good at mathematics, proofs are often wrong.
Don't trust me. Verify yourself!

If you find some mistake or have some comments, please feel free to open an issue or a PR.
\begin{flushright}
\href{https://github.com/goropikari}{goropikari}\\
\end{flushright}

\noindent From the second author:\\
\indent I'm a mathematician relatively new to quantum information theory as of the adoption of this repo, so hope to supplement the original author's work by checking and formalizing the mathematics, overly at times, while I use the task to learn the field.  The original author's sentiments about self-verification are echoed.
\begin{flushright}
\href{https://github.com/tlesaul2}{tlesaul2}
\end{flushright}

\clearpage

{\let\cleardoublepage\clearpage\tableofcontents\thispagestyle{empty}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\frontmatter""", file=F)

        print(r"\begin{tabular}{ccc}",file=F)

        print(f"{Bmatrix(permutation).strip()}", file=F)

        print(r" & = & \begin{tabular}{c} " + f"{tikz(generate_qpic(indices[permutation][0])).strip()}" + r" \\" + "\n"+ f"{tikz(generate_qpic(indices[permutation][1]))}" + r"\end{tabular} \\", end="\n",file=F)

        for i in range(2, 2+11*3+1, 3): #len(indices[permutation]), 3):
            print(" & ".join([f"{tikz(generate_qpic(indices[permutation][j])).strip()}" for j in range(i,min([i+3, len(indices[permutation])]))])+r" \\", end="\n", file=F)
            old_i = i+3

        print(r"\end{tabular}", file=F)
        print(r"\newpage", file=F)

        while old_i < len(indices[permutation]):
            print(r"\begin{tabular}{ccc}", file=F)
            for i in range(old_i, min([old_i+14*3,len(indices[permutation])]), 3):
                print(" & ".join([f"{tikz(generate_qpic(indices[permutation][j])).strip()}" for j in range(i,min([i+3, len(indices[permutation])]))])+r" \\", end="\n", file=F)
                old_i = i+3
            print(r"\end{tabular}",file=F)
            print(r"\newpage", file=F)


        print(r"\end{document}", file=F)




exit(0)

gate_counts = Counter()
length_counts = defaultdict(int)

previous = set([identity,])
current = set()
permutations.remove(identity)
length_counts[0] = 1
gate_counts[identity] = 0
length = 0
indices = defaultdict(list)
representations = defaultdict(lambda: defaultdict(int))
representations[0][identity] = 1
print(f"Found {length_counts[length]} permutation that can be expressed by {length} Toffolis and/or CNOTs")
while permutations:
    for gate in Gates:
        for perm in previous:
            permutation = perm*gate
            if permutation in permutations:
                gate_counts[permutation] = length+1
                length_counts[length+1] += 1
                permutations.remove(permutation)
            if gate_counts[permutation] == length+1:
                current.add(permutation)
                representations[length+1][permutation] += representations[length][perm]
                indices[permutation] = indices[perm]+[Gates.index(gate),]

    print(f"Found {length_counts[length+1]} permutations that can be expressed by {length+1} Toffolis and/or CNOTs")

    length+=1
    previous=current
    current=set()

eight_gate_permutations = list(map(lambda x: x[0], gate_counts.most_common(6)))
for permutation in eight_gate_permutations:
    print(f"{permutation} is a product of {representations[8][permutation]} distinct sequences of 8 gates, one of which is indexed by {indices[permutation]}")

#with open("eight_gate_permutations.pkl", "wb") as F:
    #pickle.dump(eight_gate_permutations, F)




#Found 1 permutations that can be expressed by 0 Toffolis and/or CNOTs
#Found 9 permutations that can be expressed by 1 Toffolis and/or CNOTs
#Found 60 permutations that can be expressed by 2 Toffolis and/or CNOTs
#Found 261 permutations that can be expressed by 3 Toffolis and/or CNOTs
#Found 845 permutations that can be expressed by 4 Toffolis and/or CNOTs
#Found 1784 permutations that can be expressed by 5 Toffolis and/or CNOTs
#Found 1688 permutations that can be expressed by 6 Toffolis and/or CNOTs
#Found 386 permutations that can be expressed by 7 Toffolid
