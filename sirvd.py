from time import time
# Выполнили Кулаков Алексей 3381 Георгий Козлов 3381
import matplotlib.pyplot as plt

# Функция для ввода данных пользователем
def user_input(): 
    beta = float(input())   # Коэффициент передачи
    gamma = float(input())   # Коэффициент выздоровления
    delta = float(input())  # Коэффициент смертности
    alpha = float(input())  # Скорость вакцинации
    sigma = float(input())  # Утрата иммунитета
    n = int(input())      # Общее население
    inf = int(input()) # Начальное количество инфицированных
    h = float(input()) # Шаг времени
    t = int(input()) # Общее количество времени
    return [n, inf, h, t, beta, gamma, delta, alpha, sigma] # Возвращаем введённые данные массивом

# Функция для расчёта производных S, I, R, V, D
def derivatives(S, I, R, V, D, N, beta, gamma, delta, alpha, sigma):
    dS = -beta * S * I / N + sigma * R - alpha * S  # Производная для восприимчивых
    dI = beta * S * I / N - gamma * I - delta * I   # Производная для инфицированных
    dR = gamma * I - sigma * R                     # Производная для выздоровевших
    dV = alpha * S                                 # Производная для вакцинированных
    dD = delta * I                                 # Производная для умерших
    return [dS, dI, dR, dV, dD]                    # Возвращаем производные в виде списка

# Реализация одного шага метода Рунге-Кутты (4-го или 3/8 порядка)
def rk4_step(y, h, N, koef1, koef2, beta, gamma, delta, alpha, sigma):
    S, I, R, V, D = y  # Распаковываем текущие значения переменных
    k1 = derivatives(S, I, R, V, D, N, beta, gamma, delta, alpha, sigma)
    for i in range(5):
        k1[i] *= h  # Умножаем производные на шаг времени h
    k2 = derivatives(S + k1[0] / 2, I + k1[1] / 2, R + k1[2] / 2, V + k1[3] / 2, D + k1[4] / 2, N, beta, gamma, delta, alpha, sigma)
    for i in range(5):
        k2[i] *= h
    k3 = derivatives(S + k2[0] / 2, I + k2[1] / 2, R + k2[2] / 2, V + k2[3] / 2, D + k2[4] / 2, N, beta, gamma, delta, alpha, sigma)
    for i in range(5):
        k3[i] *= h
    k4 = derivatives(S + k3[0], I + k3[1], R + k3[2], V + k3[3], D + k3[4], N, beta, gamma, delta, alpha, sigma)
    for i in range(5):
        k4[i] *= h
    for i in range(5):  # Вычисляем новое значение по формуле Рунге-Кутты
        y[i] += (k1[i] + koef1 * k2[i] + koef1 * k3[i] + k4[i]) / koef2
    return y

# Основная функция программы
def main():
    print('Введите данные')  # Пользователь вводит параметры модели
    n, inf, h, t, beta, gamma, delta, alpha, sigma = user_input()
    print('Введите номер метода: 0/любое другое число (rk4, rk3/8)')  # Выбор метода Рунге-Кутты
    num = int(input())
    
    # Засекаем время выполнения
    start = time()
    
    # Устанавливаем коэффициенты в зависимости от выбранного метода
    koef1, koef2 = 2, 6
    if num:
        koef1, koef2 = 3, 8

    # Формируем временные точки
    t_values = []
    t0 = 0
    while t0 < t:
        t_values.append(t0)
        t0 += h

    # Инициализация списков для хранения результатов
    S_values, I_values, R_values, V_values, D_values = [], [], [], [], []
    y = [n - inf, inf, 0, 0, 0]  # Начальные значения переменных

    # Цикл интегрирования по времени
    for i in t_values:
        S_values.append(y[0])
        I_values.append(y[1])
        R_values.append(y[2])
        V_values.append(y[3])
        D_values.append(y[4])
        y = rk4_step(y, h, n, koef1, koef2, beta, gamma, delta, alpha, sigma)

    # Засекаем время окончания
    end = time()

    # Построение графиков
    plt.figure(figsize=(10, 6))
    plt.plot(t_values, S_values, label="Восприимчивые (S)")
    plt.plot(t_values, I_values, label="Инфицированные (I)")
    plt.plot(t_values, R_values, label="Выздоровевшие (R)")
    plt.plot(t_values, V_values, label="Вакцинированные (V)")
    plt.plot(t_values, D_values, label="Умершие (D)")
    plt.xlabel("Время")
    plt.ylabel("Популяция")
    plt.legend()
    plt.title("SIRVD-модель (RK)")
    plt.grid()
    plt.show()

    # Вывод времени выполнения
    print(end - start)

# Запуск программы
main()
