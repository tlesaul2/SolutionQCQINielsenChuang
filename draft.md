#draft
- 11.6
$$
    H(Y) + H (X, Y, Z) - H(X, Y) - H(Y, Z)\\
        = \sum_{x,y,z} p(x,y,z) \log \left( p(x,y)p(y,z)/p(y)p(x,y,z) \right) \\
        \le \frac{1}{\ln{2}} \sum_{x,y,z} p(x,y,z) \left[1-p(x,y)p(y,z)/p(y)p(x,y,z) \right]\\
        = \frac{1-1}{\ln{2}}
        = 0
$$
The equality occurs if and only if $p(x,y)p(y,z)/p(y)p(x,y,z)=1$, which means a Markov chain condition of $Z \rightarrow Y \rightarrow X$,which is $p(x|y)=p(x|y,z)$

- 12.3
Equality of (12.14) will happen when $\rho_x$ has orthogonal support. It is obvious that n qubits have at most n orthogonal $\rho_x$s, and from (12.6),
$$
H(X:Y) \le \chi(\rho) \le H(X) \le n
$$
So, n qubits can be used to at most n bits of classical information.

- 12.31
Eve makes her qubits entangled with $|\beta_{00}\rangle$, and gets $\rho^E$.
$$
  |ABE\rangle = U|\beta_{00}^{\otimes n}\rangle |0 \rangle _E\\
  \rho^E = tr_{AB} (|ABE \rangle \langle ABE|)
$$
Note that Eve's mutual information with Alice and Bob measurements does not depend on whether Eve measures $\rho^E$ before Alice and Bob's measurement or after.
So we can assume that Eve measures $\rho^E$ after Alice and Bob's measurement.
Alice and Bob measure their Bell state, getting binary string $\vec{k}$ as an outcome.
Let $\rho^E_k$ and $p_k$ are the corresponding Eve's states and probabilities.
Note,
$$
  \rho_E = \sum_k p_k \rho^E_k.
$$
Let $K$ is a variable of $\vec{k}$ and $e$ is an outcom of a measurement of $\rho^E$, and $E$ is its variable.  From Holevo bound,
$$
H(K:E) \le S(\rho^E) - \sum_k p_k \rho^E_k \le S(\rho^E) = S(\rho).
$$
