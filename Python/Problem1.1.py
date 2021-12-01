#!/usr/bin/env python3

# standard library modules
import random
import itertools
import argparse

def balanced_extension_probability(n, k):
    """
    This function calculates the probability that under the randomly selected function,
    the k-th evaluation agrees with the previous k-1.  It calculates P(B | k)
    """
    if k <= 1:
        # There are no previous evaluations to agree with when k <= 0, so all functions
        # satisfy the requirements vacuously
        return 1
    else:
        return (n-2*k+2)/(2*n-2*k+2)

def calculate_theoretical(n, k):
    """
    This function calculates the posterior probability that a function uniformly chosen
    to be constant or balanced is in-fact constant given k evaluations in agreement of
    n possible input values.  It uses Bayesian inference.  This calculates P(C | k).
    """

    if k == 0:
        # the function is recursive.  n is fixed, so we need only a value of k to start
        # from.  k==1 corresponds to the initial evaluation, which does not provide enough
        # information to change the probability that the function is constant vice
        # balanced.  So, the a priori prior probability of 1/2 is still valid.
        return 1/2

    # Bayeseian updating requires knowning the previous probability, so calculate it
    # (recursively)
    prob_given_one_fewer_evaluations = calculate_theoretical(n,k-1)

    # use the previous probability and the probability that the k-th evaluation agrees with the
    # previous k-1 to calculate a new posterior probability that the function is constant given
    # k evaluations in agreement
    return prob_given_one_fewer_evaluations / \
            (prob_given_one_fewer_evaluations + (1-prob_given_one_fewer_evaluations)*
                                                (balanced_extension_probability(n, k)))


def main(args):
    """
    This function runs experiments to empirically determine the probability that chosen
    function equally likely to be constant or balanced is in-fact constant, given k
    evaluations in agreement.  It compares the empirical value to a theoretically
    calculated value.
    """

    # a successful experiment counter, and a set of functions to choose from
    num_experiments_in_which_a_balanced_functions_satisfied_hypotheses = 0
    if not args.conserve_memory:
        combinations = list(itertools.combinations(range(args.domain_size),args.domain_size//2))

    for i in range(args.num_experiments):
        # To be 100% faithful, we should randomly decide whether the function is balanced
        # or constant, and tally the constant functions separately, but this only doubles
        # the size of the loop.  Instead, we use this loop to choose num_experiments
        # balanced functions, and *assume* an equal number of experiments which chose
        # constant functions were also performed.

        # select a balanced function, either by randomly sampling range(domain_size), or
        # from a pre-computed list.
        if args.conserve_memory:
            combination = set(random.sample(range(args.domain_size),args.domain_size//2))
        else:
            combination = set(random.choice(combinations))
        if set(range(args.num_evaluations)) <= combination or \
                (not set(range(args.num_evaluations)) & combination):
            num_experiments_in_which_a_balanced_functions_satisfied_hypotheses+=1

    # implement assumption that an equal number of experiments which chose a constant
    # functions were also performed.
    num_experiments_selecting_a_constant_functions = args.num_experiments

    # calculate empirical probability
    empirical_probability = num_experiments_selecting_a_constant_functions / \
                                (num_experiments_selecting_a_constant_functions +
                                 num_experiments_in_which_a_balanced_functions_satisfied_hypotheses)

    # calculate theoretical probability
    theoretical_probability = calculate_theoretical(args.domain_size, args.num_evaluations)

    # calculate the error rate
    error = abs(empirical_probability-theoretical_probability)/theoretical_probability

    # report
    print("The probability that f is constant, given {} of {} equal evaluations, is {}" \
                        .format(args.num_evaluations, args.domain_size, empirical_probability))
    print("Compare to the theoretical value: {}".format(theoretical_probability))
    print("These values are off by {:.2%}".format(error))
    if error > 0.005:
        print("This error seems high.  Try rerunning with a larger -N/--num-experiments value")


if __name__ == "__main__":
    # set up argument parser
    parser = argparse.ArgumentParser("This script iteratively calculates the theoretical "
                                     "probability that a chosen function equally likely to be "
                                     "constant or balanced is in-fact constant, given k "
                                     "evaluations in agreement.  It also simulates the chose and "
                                     "evaluation to experimentally confirm it's results.  It "
                                     "accompanies problem 1.1 of Mike & Ike's 'Quantum Computation "
                                     "and Quantum Information: 10th Anniversary Edition'")
    parser.add_argument("-n","--domain-size",
                        type=int,
                        required=False,
                        default=20,
                        help="Size of the domain of the function. [Default = %(default)s]. "
                        "CAUTION: By default, this script will construct a *list* of the "
                        "domain_size \choose domain_size/2 balanced bipartitions of "
                        "1,...,domain_size.  If domain_size is too large, this will consume "
                        "exhorbitant amounts of RAM.  To prevent this, use the --conserve-memory "
                        "option, but be warned, it is *significantly* slower")
    parser.add_argument("-k","--num_evaluations",
                        type=int,
                        required=False,
                        default=8,
                        help="Number of evaluations in agreement to assume.  "
                             "[Default = %(default)s]")
    parser.add_argument("-N","--num-experiments",
                        type=int,
                        required=False,
                        default=1_000_000,
                        help="Number of experiments to perform during empirical calculation of "
                             "probability.  [Default = %(default)s].  CAUTION: with larger values "
                             "of num_evaluations, num_experiments should increase.  The results "
                             "depend on there being a significant number of experiments in which a "
                             "balanced function is selected whose evaluation produces seemingly "
                             "constant results.  With larger n_evaluations, there are factorially "
                             "fewer balanced functions to contribute to the total.")
    parser.add_argument("--conserve-memory",
                        action="store_true",
                        required=False,
                        default=False,
                        help="Conserve memory by avoiding the construction of a *list* of the "
                             "domain_size \choose domain_size/2 balanced bipartitions of "
                             "1,...,domain_size.  Instead, randomly sample the balanced bipartion "
                             "from range(domain_size).  CAUTION: this is *significantly* slower, "
                             "but may be necessary for larger values of domain_size.")

    # parse arguments
    args=parser.parse_args()

    # validate arguments
    if args.domain_size < 2:
        print("ERROR: the domain_size must be a possitive event integer")
        exit(-1)
    if args.domain_size % 2 != 0:
        print("ERROR: the domain_size must be even")
        exit(-1)

    if args.num_evaluations < 0:
        print("ERROR: evaluation cannot occur a negative number of times")
        exit(-1)

    if args.num_evaluations <= 1:
        # NOTE: validation is not strictly necessary here.  The fact that the probability is still
        #       1/2 can be verified both empirically and theoretically without issue by allowing
        #       main to be called.  We at least print a message detailing the predetermined value
        print("Having evaluated the function no more than once, there is not enough information to "
              "update the prior probability of 1/2 that the function is constant.  We still "
              "perform the experiments and calculate the probability (0.5) theoretically for "
              "completeness and robustness, but the theoretical value in the case of 0 and 1 "
              "evaluations is essentially hard-coded to be 1/2.  The emprical probability will be "
              "calculated faithfully, however.")
    if args.num_evaluations > args.domain_size//2:
        # NOTE: validation is not strictly necessary here.  The fact that the function must be
        #       constant can be verified both empirically and theoretically without issue by
        #       allowing main to be called.  We at least print a message detailing the predetermined
        #       value
        print("Having found the function to agree on more than half of its inputs, it must be "
              "constant.  We will still perform the experiments and calculate the probability "
              "(1.0) theoretically, for completeness and robustness")

    # execute
    main(args)
