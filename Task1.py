import numpy as np

"""
TASK 1
Calculation by Gauss-Seidel method
"""
"GIVEN PARAMETERS"
# Bus parameters
V1 = 1.05 + 0j

S2 = 0.6
S3 = 1.2 - 0.6j
S4 = 0.4
S5 = 0.8 - 0.4j
S6 = 0.5

S = np.array([0, S2, S3, S4, S5, S6], dtype=complex)

# Line parameters
Z12 = 0.04 + 0.12j
Z13 = 0.08 + 0.25j
Z15 = 0.1 + 0.3j
Z24 = 0.07 + 0.2j
Z26 = 0.06 + 0.18j
Z35 = 0.05 + 0.15j
Z45 = 0.09 + 0.27j
Z46 = 0.04 + 0.1j
Z56 = 0.03 + 0.09j

Z = {
    (1, 2): Z12,
    (1, 3): Z13,
    (1, 5): Z15,
    (2, 4): Z24,
    (2, 6): Z26,
    (3, 5): Z35,
    (4, 5): Z45,
    (4, 6): Z46,
    (5, 6): Z56,
}



"HELPER FUNCTIONS"
# Compute admittance
def find_admittance(Z):
    return 1 / Z

# Compute Ybus
def Ybus_calc(Y):

    # Determine the size of the matrix
    max_size = 0
    for(i, j) in Y.keys():
        max_size = max(max_size, i, j)

    Ybus = np.zeros((max_size, max_size), dtype = complex)

    # Compute final matrix
    for (i, j), y in Y.items():
        i -= 1
        j -= 1

        Ybus[i, i] += y
        Ybus[j, j] += y

        Ybus[i, j] -= y
        Ybus[j, i] -= y

    return Ybus

# Voltage by Gauss-Seidel method of a single iteration
def gauss_seidel(V, S, Ybus, slack_bus, PV_buses, PQ_buses):

    n = len(V)
    V_new = V.copy()

    for i in range(n):

        # Skip slack bus
        if i == slack_bus:
            continue

        summ = 0
        for j in range(n):
            if j != i:
                if j < i:
                    summ += Ybus[i, j] * V_new[j]
                else:
                    summ += Ybus[i, j] * V[j]

        # PQ BUS
        if i in PQ_buses:
            V_new[i] = (np.conjugate(S[i]) / np.conjugate(V_new[i]) - summ) / Ybus[i, i]

        # PV BUS
        elif i in PV_buses:

            # Step 1: Compute Q from current voltages
            I_i = np.dot(Ybus[i, :], V_new)
            Q_i = -np.imag(np.conjugate(V_new[i]) * I_i)

            # Step 2: Form updated S = P + jQ
            S_i = S[i].real + 1j * Q_i

            # Step 3: Update voltage
            V_new[i] = (np.conjugate(S_i) / np.conjugate(V_new[i]) - summ) / Ybus[i, i]

    return V_new

# Make GS function to do multiple iterations
def gauss_seidel_iterations(V, S, Y, slack_bus, PV_buses, PQ_buses, iter):
    
    V_old = V.copy()

    for i in range(iter):
        V_new = gauss_seidel(V_old, S, Y, slack_bus, PV_buses, PQ_buses)

        V_old = V_new.copy()
    
    return V_new

def print_results(V, slack_bus, iter):
    max_size = len(V)
    for i in range(max_size):
        print(f"V{i+1} value after iteration {iter}: magnitude {abs(V[i])} p.u., angle {np.angle(V[i], deg=True)} degrees")



"""
SOLUTION
"""

"(a)"

# Define slack, PV, and PQ buses, considering the numbering starts with zero
slack_bus = 0
PV_buses = [1, 3, 5]
PQ_buses = [2, 4]

# Compute admittances
Y = {}

for key, z in Z.items():
    Y[key] = find_admittance(z)

# Find Ybus matrix
Ybus = Ybus_calc(Y)

# Define initial values
V2 = V3 = V4 = V5 = V6 = 1.05 + 0j

V = np.array([V1, V2, V3, V4, V5, V6], dtype=complex)

# Calculate voltages after one iteration by Gauss-Seidel method
V_1 = gauss_seidel(V, S, Ybus, slack_bus, PV_buses, PQ_buses)

# Print results
print("\n")
print("Task 1.a results")
print_results(V_1, slack_bus, 1)
print("\n")


"(b)"
# Calculate voltages after two iterations by Gauss-Seidel method
V_2 = gauss_seidel_iterations(V, S, Ybus, slack_bus, PV_buses, PQ_buses, 2)

# Print results
print("Task 1.b results")
print_results(V_2, slack_bus, 2)
print("\n")



"(c)"
"""
As Bus 6 becomes a slack bus
"""
# Bus parameters
V6 = 1.05 + 0j

S = np.array([0, S2, S3, S4, S5, 0], dtype=complex)

# Define slack, PV, and PQ buses, considering the numbering starts with zero
slack_bus = 5
PV_buses = [0, 1, 3]
PQ_buses = [2, 4]

# Compute voltages by GS method after two iterations
V1 = V2 = V3 = V4 = V5 = 1.05 + 0j

V = np.array([V1, V2, V3, V4, V5, V6], dtype=complex)

V_2 = gauss_seidel_iterations(V, S, Ybus, slack_bus, PV_buses, PQ_buses, 2)

# Print results
print("Task 1.c results")
print_results(V_2, slack_bus, 2)
print("\n")