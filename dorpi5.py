import matplotlib.pyplot as plt

# Параметры модели
beta = 0.3  # Коэффициент передачи
gamma = 0.1  # Коэффициент выздоровления
delta = 0.05  # Коэффициент смертности
alpha = 0.02  # Скорость вакцинации
sigma = 0.01  # Утрата иммунитета
N = 1000  # Общее население

# Начальные условия
S0 = 990
I0 = 10
R0 = 0
V0 = 0
D0 = 0
initial_conditions = [S0, I0, R0, V0, D0]

# Временные параметры
t_start = 0
t_end = 100
dt = 0.1
tol = 1e-6  # Допустимая ошибка

def model(t, y):
    S, I, R, V, D = y
    dSdt = -beta * I * S / N + sigma * R - alpha * S
    dIdt = beta * I * S / N - gamma * I - delta * I
    dRdt = gamma * I - sigma * R
    dVdt = alpha * S
    dDdt = delta * I
    return [dSdt, dIdt, dRdt, dVdt, dDdt]

def vector_add(v1, v2):
    return [v1[i] + v2[i] for i in range(len(v1))]

def vector_scale(v, scalar):
    return [vi * scalar for vi in v]

def vector_dot(v1, v2):
    return sum(v1[i] * v2[i] for i in range(len(v1)))

def dopri5_step(t, y, dt, model):
    # Коэффициенты метода DOPRI5
    a = [0, 1/5, 3/10, 4/5, 8/9, 1, 1]
    b = [
        [0, 0, 0, 0, 0, 0],
        [1/5, 0, 0, 0, 0, 0],
        [3/40, 9/40, 0, 0, 0, 0],
        [44/45, -56/15, 32/9, 0, 0, 0],
        [19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0],
        [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0],
        [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84]
    ]
    c = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]
    c_hat = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40]

    k = [[0] * len(y) for _ in range(7)]
    k[0] = model(t, y)
    for i in range(1, 7):
        y_temp = vector_add(y, vector_scale([vector_dot(b[i][:i], [k[j] for j in range(i)]) for _ in range(len(y))], dt))
        k[i] = model(t + a[i] * dt, y_temp)

    y_new = vector_add(y, vector_scale([vector_dot(c, [k[j][i] for j in range(7)]) for i in range(len(y))], dt))
    y_hat = vector_add(y, vector_scale([vector_dot(c_hat, [k[j][i] for j in range(7)]) for i in range(len(y))], dt))
    error = sum((y_new[i] - y_hat[i]) ** 2 for i in range(len(y))) ** 0.5

    return y_new, error

def dopri5(t_start, t_end, y0, model, dt, tol):
    t = [t_start]
    y = [y0]
    while t[-1] < t_end:
        y_new, error = dopri5_step(t[-1], y[-1], dt, model)
        if error < tol:
            t.append(t[-1] + dt)
            y.append(y_new)
        dt = 0.9 * dt * (tol / error) ** 0.2
    return t, y

# Решение системы ОДУ
t, y = dopri5(t_start, t_end, initial_conditions, model, dt, tol)

# Визуализация результатов
plt.plot(t, [yi[0] for yi in y], label='S(t)')
plt.plot(t, [yi[1] for yi in y], label='I(t)')
plt.plot(t, [yi[2] for yi in y], label='R(t)')
plt.plot(t, [yi[3] for yi in y], label='V(t)')
plt.plot(t, [yi[4] for yi in y], label='D(t)')
plt.xlabel('Time')
plt.ylabel('Population')
plt.legend()
plt.show()
