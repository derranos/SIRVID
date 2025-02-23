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
s1,i1,r1,v1,d1 = rk.main(n, inf, h, t, beta, gamma, delta, alpha, sigma)
print(line)
print("GOIDA NUMBER 2) RK3/8")
s2,i2,r2,v2,d2 = rk3_8.main(n, inf, h, t, beta, gamma, delta, alpha, sigma)
print(line)
print("GOIDA NUMBER 3) DOPRI5")
s3,i3,r3,v3,d3 = dopri5.dorpi_5(beta, gamma, delta, alpha, sigma, n, t)
print(line)
print("GOIDA NUMBER 3) DOPRI8")
s4,i4,r4,v4,d4 = dopri8.dopri_8(beta, gamma, delta, alpha, sigma, n, t)
print(f'Разница RK4: {abs(s1 - s4), abs(i1 - i4), abs(r1 - r4), abs(v1 - v4), abs(d1-d4)}')
print(f'Разница RK3/8: {abs(s2 - s4), abs(i2 - i4), abs(r2 - r4), abs(v2 - v4), abs(d2-d4)}')
print(f'Разница DOPRI5: {abs(s3 - s4), abs(i3 - i4), abs(r3 - r4), abs(v3 - v4), abs(d3-d4)}')