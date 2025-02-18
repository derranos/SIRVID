import matplotlib.pyplot as plt

# Коэффициенты метода DOPRI8
A = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1/18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1/48, 1/16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1/32, 0, 3/32, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5/16, 0, -75/64, 75/64, 0, 0, 0, 0, 0, 0, 0, 0],
    [3/80, 0, 0, 3/16, 3/20, 0, 0, 0, 0, 0, 0, 0],
    [29443841/614563906, 0, 0, 77736538/692538347, -28693883/1125000000, 23124283/1800000000, 0, 0, 0, 0, 0, 0],
    [16016141/946692911, 0, 0, 61564180/158732637, 22789713/633445777, 545815736/2771057229, -180193667/1043307555, 0, 0, 0, 0, 0],
    [39632708/573591083, 0, 0, -433636366/683701615, -421739975/2616292301, 100302831/723423059, 790204164/839813087, 800635310/3783071287, 0, 0, 0, 0],
    [246121993/1340847787, 0, 0, -37695042795/15268766246, -309121744/1061227803, -12992083/490766935, 6005943493/2108947869, 393006217/1396673457, 123872331/1001029789, 0, 0, 0],
    [-1028468189/846180014, 0, 0, 8478235783/508512852, 1311729495/1432422823, -10304129995/1701304382, -48777925059/3047939560, 15336726248/1032824649, -45442868181/3398467696, 3065993473/597172653, 0, 0],
    [185892177/718116043, 0, 0, -3185094517/667107341, -477755414/1098053517, -703635378/230739211, 5731566787/1027545527, 5232866602/850066563, -4093664535/808688257, 3962137247/1805957418, 65686358/487910083, 0]
]

B = [14005451/335480064, 0, 0, 0, 0, -59238493/1068277825, 181606767/758867731, 561292985/797845732, -1041891430/1371343529, 760417239/1151165299, 118820643/751138087, -528747749/2220607170, 1/4]  # Коэффициенты для основного решения
B_hat = [13451932/455176623, 0, 0, 0, 0, -808719846/976000145, 1757004468/5645159321, 656045339/265891186, -3867574721/1518517206, 465885868/322736535, 53011238/667516719, 2/45, 0]  # Коэффициенты для оценки ошибки
C = [0, 1/18, 1/12, 1/8, 5/16, 3/8, 59/400, 93/200, 5490023248/9719169821, 13/20, 1201146811/1299019798, 1, 1]  # Временные точки для стадий

# Параметры модели
beta = 0.5  # Коэффициент передачи
gamma = 0.1  # Коэффициент выздоровления
delta = 0.01  # Коэффициент смертности
alpha = 0.02  # Скорость вакцинации
sigma = 0.01  # Утрата иммунитета
N = 1000  # Общее население

# Система ОДУ
def system(t, y):
    S, I, R, V, D = y
    dSdt = -beta * I * S / N + sigma * R - alpha * S
    dIdt = beta * I * S / N - gamma * I - delta * I
    dRdt = gamma * I - sigma * R
    dVdt = alpha * S
    dDdt = delta * I
    return [dSdt, dIdt, dRdt, dVdt, dDdt]

# Реализация метода DOPRI8
def dopri8(f, t_span, y0, tol=1e-6, h_min=1e-4, h_max=1.0):
    t, y = t_span[0], y0
    h = h_max  # Начальный шаг
    solution = [(t, y.copy())]

    while t < t_span[1]:
        if t + h > t_span[1]:
            h = t_span[1] - t  # Корректируем шаг, чтобы не выйти за пределы интервала

        # Вычисляем стадии метода Рунге-Кутты
        k = [[0] * len(y0) for _ in range(13)]
        k[0] = f(t, y)
        for i in range(1, 13):
            y_stage = [y[j] + h * sum(A[i-1][m] * k[m][j] for m in range(i)) for j in range(len(y0))]
            k[i] = f(t + C[i] * h, y_stage)

        # Основное решение (8-го порядка)
        y_new = [y[j] + h * sum(B[m] * k[m][j] for m in range(13)) for j in range(len(y0))]

        # Оценка ошибки (7-го порядка)
        y_hat = [y[j] + h * sum(B_hat[m] * k[m][j] for m in range(13)) for j in range(len(y0))]
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
t_span = [0, 100]

# Решение системы ОДУ
# tol - допустимая ошибка
solution = dopri8(system, t_span, y0, tol=1e-6)

# Извлечение данных для построения графиков
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
plt.title("Динамика эпидемии (DOPRI8)")
plt.legend()
plt.grid()
plt.show()