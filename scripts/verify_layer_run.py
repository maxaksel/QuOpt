# -*- coding: utf-8 -*-

"""This module contains a script to ensure that all optimization runs for a
particular set of hyperparameters were performed.
"""
import re
from collections import defaultdict


def main():
    # User input
    print("How many layers did you simulate?")
    num_layers = input("> ")
    num_layers = int(num_layers)
    print("How many initial guess points did you use for each optimization problem?")
    num_guess_points = input("> ")
    num_guess_points = int(num_guess_points)
    print("What is the name of your stdout file?")

    cn_dict = defaultdict(lambda: 0)

    with open(stdout_filename, "r") as stdout_file:
        for line in stdout_file:
            sline = line.split(" ")
            if sline[1] != "circuit":
                continue
            cn = int(sline[2])
            cn_dict[cn] += 1

    for circuit_number in range(3**num_layers - 1):
        if cn_dict[circuit_number] != num_guess_points:
            print("Error! Circuit missed OR circuit optimization not run the proper number of times!")
            print(f"Problem with circuit number {circuit_number}")
            return

    print("All okay.")


if __name__ == '__main__':
    main()
