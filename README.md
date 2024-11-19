# " Design and Implementation of Active Filter using Genetic Algorithm "

This project focuses on automating the design of active filters through genetic algorithms implemented in MATLAB. There is a MATLAB GUI designed for gathering initial values of preference. By analyzing 16 circuit topologies, we found the transfer function of each circuit to formulate the fitness function, which optimizes the values of resistors and capacitors. Practical testing of six unique filters demonstrated the effectiveness of these methods. The outcomes suggest that this automated approach can enhance efficiency in filter design, offering valuable insights for both academic research and industrial applications.

Note: This project is coded and tested in MATLAB R2019b.

Features:

- Gain
- Different initial evaluations
- Number of iterations: This number indicates how many times each stage is designed. The best one will be chosen as a result.
- Op-amp IC selection
- Topology type of circuit: low-pass, band-pass, high-pass
- Method of solving: Butterworth, Chebyshev 1, Chebyshev 2, and elliptic
- Filter type: low-pass, band-pass, high-pass
- Fitness function stop: fitness limitation and maximum generation
- Elements standard: custom selection and E24
