import matplotlib.pyplot as plt

# Коэффициенты метода DOPRI5
A = [
    [0, 0, 0, 0, 0, 0],
    [1/5, 0, 0, 0, 0, 0],
    [3/40, 9/40, 0, 0, 0, 0],
    [44/45, -56/15, 32/9, 0, 0, 0],
    [19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0],
    [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0]
]

B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84]  # Коэффициенты для основного решения
B_hat = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40]  # Коэффициенты для оценки ошибки
C = [0, 2/10, 3/10, 4/10, 8/10, 8/9, 1]  # Временные точки для стадий

# Параметры модели
beta = 0.2  # Коэффициент передачи
gamma = 0.01  # Коэффициент выздоровления
delta = 0.01  # Коэффициент смертности
alpha = 0.001  # Скорость вакцинации
sigma = 0.05  # Утрата иммунитета
N = 1000  # Общее население

# Система ОДУ
def system(y):
    S, I, R, V, D = y
    dSdt = -beta * I * S / N + sigma * R - alpha * S
    dIdt = beta * I * S / N - gamma * I - delta * I
    dRdt = gamma * I - sigma * R
    dVdt = alpha * S
    dDdt = delta * I
    return [dSdt, dIdt, dRdt, dVdt, dDdt]

# Реализация метода DOPRI5
def dopri5(f, t_span, y0, tol=1e-6, h_min=1e-4, h_max=1.0):
    t, y = t_span[0], y0
    h = h_max  # Начальный шаг
    solution = [(t, y.copy())]

    while t < t_span[1]:
        if t + h > t_span[1]:
            h = t_span[1] - t  # Корректируем шаг, чтобы не выйти за пределы интервала

        # Вычисляем стадии метода Рунге-Кутты
        k = [[0] * len(y0) for _ in range(7)] # явные уравнения коэф Рунге-Куты
        k[0] = f(t, y)
        for i in range(1, 7):
            y_stage = []
            for j in range(len(y0)):
                koef = sum(A[i-1][m] * k[m][j] for m in range(i))
                y_stage.append(y[j] + h * koef)
            k[i] = f(t + C[i] * h, y_stage)

        # Основное решение (5-го порядка)
        y_new = []
        for j in range(len(y0)):
            res = 0
            for m in range(6):
                res += B[m] * k[m][j]
            y_new.append(y[j] + h * res)

        # Оценка ошибки (4-го порядка)
        y_hat = [y[j] + h * sum(B_hat[m] * k[m][j] for m in range(7)) for j in range(len(y0))]
        error = max(abs(y_new[j] - y_hat[j]) for j in range(len(y0)))  # Максимальная ошибка

        # Адаптация шага
        if error < tol:
            t += h
            y = y_new
            solution.append((t, y.copy()))
            h = min(h_max, 0.9 * h * (tol / error)**0.2)  # Увеличиваем шаг
        else:
            h = max(h_min, 0.9 * h * (tol / error)**0.25)  # Уменьшаем шаг

    return solution

# Начальные условия 
y0 = [999, 1, 0, 0, 0]  # S(0), I(0), R(0), V(0), D(0)

# Интервал интегрирования
t_span = [0, 1000]

# Решение системы ОДУ
solution = dopri5(system, t_span, y0, tol=1e-6)

print(solution[0])
# Вывод результатов
# for t, y in solution:
#     print(f"t = {t:.2f}, Здоров. = {y[0]:.2f}, Заражен. = {y[1]:.2f}, Вылеч. = {y[2]:.2f}, Вакцин. = {y[3]:.2f}, Умер. = {y[4]:.2f}")

times = [point[0] for point in solution]
S_values = [point[1][0] for point in solution]
I_values = [point[1][1] for point in solution]
R_values = [point[1][2] for point in solution]
V_values = [point[1][3] for point in solution]
D_values = [point[1][4] for point in solution]

# Построение графиков
plt.figure(figsize=(12, 8))

plt.plot(times, S_values, label="S(t) - Восприимчивые")
plt.plot(times, I_values, label="I(t) - Инфицированные")
plt.plot(times, R_values, label="R(t) - Выздоровевшие")
plt.plot(times, V_values, label="V(t) - Вакцинированные")
plt.plot(times, D_values, label="D(t) - Умершие")

plt.xlabel("Время (t)")
plt.ylabel("Количество людей")
plt.title("Динамика эпидемии")
plt.legend()
plt.grid()
plt.show()