# Task 1

In Task 1, there is given a power network with 6 buses. The main task
was to compute bus voltages by Gauss-Seidel method:

-   a)  after 1 iteration;
-   b)  after 2 iterations;
-   c)  after 2 iterations with slack bus being Bus 6 instead.

The code in Python was developed to solve this task, which is attached
along with this report. Values of bus voltages were initialized at
1.05 + j0 p.u. to match the slack bus. The results given by the code is
given in Figure 1.

*Figure 1: Task 1 code results output*

Based on the results, from Figure 1, the following observations are
made:

-   Voltage values show low convergence after one and two iterations,
    with values significantly deviating from the initial guess values.
-   When the slack bus is changed from Bus 1 to Bus 6, values of the
    voltage show changes in both magnitude and angle values, but the
    overall solution remains consistent. It shows that change os the
    slack bus resulted in the change of the reference in the power
    system that led to differences of the results between (b) and (c).
-   Slack bus plays a role of a reference node in the power system with
    defined value of the voltage. It balances the system power,
    supplying or absorbing the mismatch in active and reactive power. In
    Gauss-Seidel method, slack bus plays a role of being the reference
    to the calculations and its change can affect numerical results.

# Task 2

In Task 2, there is given a power network with 4 buses. The main task
was to compute bus voltages by Newton-Raphson method:

-   a)  after 2 iterations;
-   b)  after 5 iterations;
-   c)  after 5 iterations with addition of two buses (one PV and one
        PQ) and change of line parameters.

The code in Python was developed to solve this task, which is attached
along with this report. to match the slack bus. The results given by the
code is given in Figure 2.

*Figure 2: Task 2 code results output*

Based on the results, from Figure 2, the following observations are
made:

-   Newton-Raphson method shows much faster convergence than
    Gauss-Seidel method with the voltage magnitude values already
    converging after 2 iterations, and more after 5 iterations.
-   Angle values generally show overall increase after iterations as
    seen from the results for (a) and (b).
-   The addition of two other buses (one PV and one PQ) with the
    alternations of the line parameters in (c) led to the significant
    changes of the voltage magnitudes and angles with wide differences
    between buses.
