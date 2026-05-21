import numpy as np


"""
TASK 2
Calculation by Newton-Raphson method
"""
"GIVEN PARAMETERS"
# Bus parameters
V1 = 1 + 0j
V2_abs = 1.05

S1 = 0 + 0j
S2 = 1.2
S3 = 0.8 + 0.3j
S4 = 0.6 + 0.2j

S = np.array([S2.real, S3.real, S4.real, S3.imag, S4.imag], dtype=float)

# Ybus matrix
Ybus = np.array([[6.1667 - 18.5j, -2.5 + 7.5j, -2 + 6j, -1.6667 + 5j],
                 [-2.5 + 7.5j, 7.7523 - 24.1743j, -2.7523 + 9.1743j, -2.5 + 7.5j],
                 [-2 + 6j, -2.7523 + 9.1743j, 7.6935 - 26.939j, -2.9412 + 11.7647j],
                 [-1.6667 + 5j, -2.5 + 7.5j, -2.9412 + 11.7647j, 7.1078 - 24.2647j]], dtype=complex)

# Unknown parameters: delta2, delta3, delta4, V3_abs, V4_abs


"HELPER FUNCTIONS"
# Calculate power values
def calculate_power(V, Ybus, slack_bus, PV_buses, PQ_buses):
    n = len(V)

    p_number = len(PV_buses) + len(PQ_buses)
    q_number = len(PQ_buses)

    P_new = np.zeros(p_number)
    Q_new = np.zeros(q_number)

    p_index = 0
    q_index = 0

    for i in range(n):

        # Skip slack bus
        if i == slack_bus:
            continue

        # For PV buses, define equation for only real power
        elif i in PV_buses:
            P = 0
            for j in range(n):
                P += abs(V[i]) * abs(V[j]) * abs(Ybus[i, j]) * np.cos(np.angle(Ybus[i, j]) + np.angle(V[j]) - np.angle(V[i]))
            P_new[p_index] = P
            p_index += 1

        # For PQ buses, define equations for both real and reactive power
        elif i in PQ_buses:
            P = 0
            Q = 0
            for j in range(n):
                P += abs(V[i]) * abs(V[j]) * abs(Ybus[i, j]) * np.cos(np.angle(Ybus[i, j]) + np.angle(V[j]) - np.angle(V[i]))
                Q -= abs(V[i]) * abs(V[j]) * abs(Ybus[i, j]) * np.sin(np.angle(Ybus[i, j]) + np.angle(V[j]) - np.angle(V[i]))
            P_new[p_index] = P
            Q_new[q_index] = Q
            p_index += 1
            q_index += 1

    S_new = np.concatenate((P_new, Q_new))
    return S_new


# Calculate power mismatches
def calculate_mismatches(S, S_new):
    return S - S_new

# Calculate Jacobian matrix
import numpy as np

def jacobian(V, Ybus, PV_buses, PQ_buses):
    n = len(V)

    # counts
    p_number = len(PV_buses) + len(PQ_buses)
    q_number = len(PQ_buses)

    J11 = np.zeros((p_number, p_number))  # dP/d_delta
    J12 = np.zeros((p_number, q_number))  # dP/dV
    J21 = np.zeros((q_number, p_number))  # dQ/d_delta
    J22 = np.zeros((q_number, q_number))  # dQ/dV

    # mapping: angle buses (PV + PQ)
    angle_buses = PV_buses + PQ_buses

    # mapping: voltage buses (PQ only)
    voltage_buses = PQ_buses

    for row_i, i in enumerate(angle_buses):

        Vi = abs(V[i])
        deltai = np.angle(V[i])

        for col_j, j in enumerate(angle_buses):

            Vj = abs(V[j])
            Yij = Ybus[i, j]
            deltaij = np.angle(Yij)

            angle = deltaij + np.angle(V[j]) - deltai

            if i == j:
                # diagonal terms (simplified standard NR form)
                P_i = 0
                Q_i = 0
                for k in range(n):
                    Vk = abs(V[k])
                    Yik = Ybus[i, k]
                    deltaik = np.angle(Yik)
                    angle_k = deltaik + np.angle(V[k]) - deltai

                    P_i += Vi * Vk * abs(Yik) * np.cos(angle_k)
                    Q_i -= Vi * Vk * abs(Yik) * np.sin(angle_k)

                J11[row_i, col_j] = -Q_i - Vi**2 * np.imag(Ybus[i, i])
            else:
                J11[row_i, col_j] = Vi * Vj * abs(Yij) * np.sin(angle)

        # J12
        for col_j, j in enumerate(voltage_buses):
            if i == j:
                # diagonal voltage effect (simplified)
                J12[row_i, col_j] = (Vi * np.real(Ybus[i, i])) + 0
            else:
                Yij = Ybus[i, j]
                deltaij = np.angle(Yij)
                angle = deltaij + np.angle(V[j]) - deltai
                J12[row_i, col_j] = Vi * abs(Yij) * np.cos(angle)

    # Q equations (J21 and J22)
    for row_i, i in enumerate(voltage_buses):

        Vi = abs(V[i])
        deltai = np.angle(V[i])

        # dQ/d_delta
        for col_j, j in enumerate(angle_buses):
            if i == j:
                J21[row_i, col_j] = - (np.real(Ybus[i, i]) * Vi**2)
            else:
                Yij = Ybus[i, j]
                deltaij = np.angle(Yij)
                angle = deltaij + np.angle(V[j]) - deltai
                J21[row_i, col_j] = -Vi * abs(Yij) * np.cos(angle)

        # dQ/dV
        for col_j, j in enumerate(voltage_buses):
            if i == j:
                J22[row_i, col_j] = - (Vi * np.imag(Ybus[i, i]))
            else:
                Yij = Ybus[i, j]
                deltaij = np.angle(Yij)
                angle = deltaij + np.angle(V[j]) - deltai
                J22[row_i, col_j] = -Vi * abs(Yij) * np.sin(angle)

    # final Jacobian
    J = np.block([[J11, J12],
                  [J21, J22]])

    return J

# Calculate voltage corrections
def voltage_corrections(mismatches, J):
    return np.linalg.solve(J, mismatches)

# Update voltage values
def update_voltages(V, corrections, PV_buses, PQ_buses):

    V_new = V.copy()

    angle_buses = PV_buses + PQ_buses
    voltage_buses = PQ_buses

    # split corrections
    d_delta = corrections[:len(angle_buses)]
    d_V = corrections[len(angle_buses):]

    # update angles
    for i, bus in enumerate(angle_buses):
        mag = abs(V_new[bus])
        ang = np.angle(V_new[bus])
        V_new[bus] = mag * np.exp(1j * (ang + d_delta[i]))

    # update magnitudes
    for i, bus in enumerate(voltage_buses):
        ang = np.angle(V_new[bus])
        V_new[bus] = (abs(V_new[bus]) + d_V[i]) * np.exp(1j * ang)

    return V_new


# Unite everything in a single NR function to perform the power flow analysis
# Fore one iteration
def newton_raphson(V, S, Ybus, slack_bus, PV_buses, PQ_buses):
    S_new = calculate_power(V, Ybus, slack_bus, PV_buses, PQ_buses)
    mismatches = calculate_mismatches(S, S_new)
    J = jacobian(V, Ybus, PV_buses, PQ_buses)
    corrections = voltage_corrections(mismatches, J)
    V_new = update_voltages(V, corrections, PV_buses, PQ_buses)

    return V_new


# For multiple iterations
def newton_raphson_iterations(V, S, Ybus, slack_bus, PV_buses, PQ_buses, iter):
    V_old = V.copy()

    for i in range(iter):
        V_new = newton_raphson(V_old, S, Ybus, slack_bus, PV_buses, PQ_buses)

        V_old = V_new.copy()

    return V_new

# print results
def print_results(V, slack_bus, iter):
    max_size = len(V)
    for i in range(max_size):
        print(f"V{i+1} value after iteration {iter}: magnitude {abs(V[i])} p.u., angle {np.angle(V[i], deg=True)} degrees")




"SOLUTION"
"(a)"

# Define slack, PV, and PQ buses, considering the numbering starts with zero
slack_bus = 0
PV_buses = [1]
PQ_buses = [2, 3]

# Define initial values
V2 = V2_abs * np.exp(1j * 0) # Initial guess for PV bus
V3 = V4 = 1 + 0j

V = np.array([V1, V2, V3, V4], dtype=complex)

# Calculate voltages after two iterations by Newton-Raphson method
V_2 = newton_raphson_iterations(V, S, Ybus, slack_bus, PV_buses, PQ_buses, 2)

# Print results
print("\n")
print("Task 2.a results")
print_results(V_2, slack_bus, 2)
print("\n")


"(b)"
# Calculate voltages after five iterations by Newton-Raphson method
V_5 = newton_raphson_iterations(V, S, Ybus, slack_bus, PV_buses, PQ_buses, 5)

# Print results
print("Task 2.b results")
print_results(V_5, slack_bus, 5)
print("\n")


"(c)"

# New given data
S5 = 0.5 + 0.2j
S6 = 0.6

V6 = 1.04 + 0j

S = np.array([S2.real, S3.real, S4.real, S5.real, S6.real, S3.imag, S4.imag, S5.imag], dtype=float)

# Line impedances
Z12 = 0.04 + 0.12j
Z13 = 0.05 + 0.15j
Z14 = 0.06 + 0.18j
Z23 = 0.03 + 0.1j
Z24 = 0.04 + 0.12j
Z34 = 0.02 + 0.08j
Z35 = 0.03 + 0.09j
Z46 = 0.04 + 0.11j
Z56 = 0.02 + 0.07j

Z = {
    (1, 2): Z12,
    (1, 3): Z13,
    (1, 4): Z14,
    (2, 3): Z23,
    (2, 4): Z24,
    (3, 4): Z34,
    (3, 5): Z35,
    (4, 6): Z46,
    (5, 6): Z56,
}

# Helper functions for this part only: calculate admitances and Ybus matrix
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

# Solution
# Define slack, PV, and PQ buses, considering the numbering starts with zero
slack_bus = 0
PV_buses = [1, 5]
PQ_buses = [2, 3, 4]

# Calculate admittances and Ybus matrix
Y = {}

for key, z in Z.items():
    Y[key] = find_admittance(z)

# Find Ybus matrix
Ybus = Ybus_calc(Y)

# Define initial values
V3 = V4 = V5 = 1 + 0j

V = np.array([V1, V2, V3, V4, V5, V6], dtype=complex)

# Calculate voltages after two iterations by Newton-Raphson method
V_5 = newton_raphson_iterations(V, S, Ybus, slack_bus, PV_buses, PQ_buses, 5)

# Print results
print("Task 2.c results")
print_results(V_5, slack_bus, 5)
print("\n")