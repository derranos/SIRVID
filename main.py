import rk
import rk3_8
import dopri5
import dopri8

print("Выберите метод для вычисления параметров")

print("GOGOGOGOIIIDDDDAAA!!!")
line = "------------------"


beta = float(input("beta: "))   # Коэффициент передачи
gamma = float(input("gamma: "))   # Коэффициент выздоровления
delta = float(input("delta: "))  # Коэффициент смертности
alpha = float(input("alpha: "))  # Скорость вакцинации
sigma = float(input("sigma: "))  # Утрата иммунитета
n = int(input("N(people): "))      # Общее население
inf = int(input("INF: ")) # Начальное количество инфицированных
h = float(input("H: ")) # Шаг времени
t = int(input("t0 - t: ")) # Общее количество времени

print("GOIDA NUMBER 1) RK4")
rk.main(n, inf, h, t, beta, gamma, delta, alpha, sigma)
print(line)
print("GOIDA NUMBER 2) RK3/8")
rk3_8.main(n, inf, h, t, beta, gamma, delta, alpha, sigma)
print(line)
print("GOIDA NUMBER 3) DOPRI5")
dopri5.dorpi_5(beta, gamma, delta, alpha, sigma, n, t)
print(line)
print("GOIDA NUMBER 3) DOPRI8")
dopri8.dopri_8(beta, gamma, delta, alpha, sigma, n, t)

