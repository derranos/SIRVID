import matplotlib.pyplot as plt

def dopri8(f, t0, y0, t_end, h):
    # Коэффициенты метода Дорманда-Принца 8-го порядка
    a = [0] * 13
    b = [0] * 13
    c = [0] * 13
    # Заполнение коэффициентов (значения должны быть взяты из метода DOPRI8)
    
    t = t0
    y = y0
    t_values = [t0]
    y_values = [y0]
    while t < t_end:
        if t + h > t_end:
            h = t_end - t
        k = [0] * 13
        for i in range(13):
            k[i] = h * f(t + c[i] * h, y + sum(a[i][j] * k[j] for j in range(i)))
        y = y + sum(b[i] * k[i] for i in range(13))
        t = t + h
        t_values.append(t)
        y_values.append(y)
    return t_values, y_values

def model(t, y):
    S, I, R, V, D = y
    N = S + I + R + V + D
    dSdt = -beta * I * S / N + sigma * R - alpha * S
    dIdt = beta * I * S / N - gamma * I - delta * I
    dRdt = gamma * I - sigma * R
    dVdt = alpha * S
    dDdt = delta * I
    return [dSdt, dIdt, dRdt, dVdt, dDdt]

# Параметры модели
beta = 0.3
gamma = 0.1
delta = 0.05
alpha = 0.02
sigma = 0.01
N = 1000

# Начальные условия
S0 = 990
I0 = 10
R0 = 0
V0 = 0
D0 = 0
y0 = [S0, I0, R0, V0, D0]

# Временные параметры
t0 = 0
t_end = 100
h = 0.1

# Решение системы ОДУ
t_values, y_values = dopri8(model, t0, y0, t_end, h)

# Извлечение результатов
S_values = [y[0] for y in y_values]
I_values = [y[1] for y in y_values]
R_values = [y[2] for y in y_values]
V_values = [y[3] for y in y_values]
D_values = [y[4] for y in y_values]

# Построение графиков
plt.figure(figsize=(10, 6))
plt.plot(t_values, S_values, label='Susceptible (S)')
plt.plot(t_values, I_values, label='Infected (I)')
plt.plot(t_values, R_values, label='Recovered (R)')
plt.plot(t_values, V_values, label='Vaccinated (V)')
plt.plot(t_values, D_values, label='Deceased (D)')
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Epidemic Model')
plt.legend()
plt.grid(True)
plt.show()
