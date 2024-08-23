# Импортируем библиотеки
from datetime import datetime, date, timedelta
import spacetrack.operators as op
from spacetrack import SpaceTrackClient
from pyorbital.orbital import Orbital
import numpy as np
import matplotlib.pyplot as plt

USERNAME = '*************************'
PASSWORD = '*************************'

# Средние энергии каналов в МэВ
E_mean = [0.047, 0.334, 0.506, 0.800, 1.265, 2.000, 3.162, 5.060, 8.000, 10.000]
# Эффективная площадь детектора см^2
square = 832
# Имена файлов
names = ['krf20090301_2_S2_bg.thr', 'krf20090302_2_S2_bg.thr', 'krf20090303_2_S2_bg.thr',
         'krf20090304_2_S2_bg.thr', 'krf20090305_2_S2_bg.thr', 'krf20090306_2_S2_bg.thr',
         'krf20090307_2_S2_bg.thr', 'krf20090308_2_S2_bg.thr', 'krf20090309_2_S2_bg.thr',
         'krf20090310_2_S2_bg.thr', 'krf20090311_2_S2_bg.thr', 'krf20090312_2_S2_bg.thr',
         'krf20090313_2_S2_bg.thr', 'krf20090314_2_S2_bg.thr', 'krf20090315_2_S2_bg.thr',
         'krf20090316_2_S2_bg.thr', 'krf20090317_2_S2_bg.thr', 'krf20090318_2_S2_bg.thr',
         'krf20090319_2_S2_bg.thr', 'krf20090320_2_S2_bg.thr', 'krf20090321_2_S2_bg.thr',
         'krf20090322_2_S2_bg.thr', 'krf20090323_2_S2_bg.thr', 'krf20090324_2_S2_bg.thr',
         'krf20090325_2_S2_bg.thr', 'krf20090326_2_S2_bg.thr', 'krf20090327_2_S2_bg.thr',
         'krf20090328_2_S2_bg.thr', 'krf20090329_2_S2_bg.thr',
         'krf20090330_2_S2_bg.thr', 'krf20090331_2_S2_bg.thr']
# Номера файлов с полным набором данных
days = [1, 2, 3, 4, 6, 19, 20, 21, 24, 28, 29, 31]
# Значения северной широты для выборки данных (два значения для сравнения результатов)
lat_num = [45, 30]
# Идентификатор КА Коронас-Фотон
identifier = 33504


# Описываем функцию получения данных формата TLE
def get_spacetrack_tle(sat_id, start_date, end_date, username, password):
    st = SpaceTrackClient(identity=username, password=password)
    # Определяем диапазон дат через оператор библиотеки:
    daterange = op.inclusive_range(start_date, end_date)
    # Выполнение запроса для диапазона дат:
    data = st.tle(norad_cat_id=sat_id, orderby='epoch desc', limit=1, format='tle',
                  epoch=daterange)

    if not data:
        return 0, 0

    tle_1 = data[0:69]
    tle_2 = data[70:139]
    return tle_1, tle_2


# Описываем функцию расчета интенсивности за одну порцию времени
def intensity_per_serving_time(file_name, line):

    # Достаем из файла значения количества отсчётов в список:
    with open(file_name) as f:
        counts = np.genfromtxt(f, dtype=np.int32, skip_header=line - 1,
                               skip_footer=21600 - line, usecols=np.arange(2, 12))

    r = sum(counts)  # Общее число отсчетов
    t_dead = r * 0.000006  # Мертвое время
    t_live = 4 - t_dead  # Живое время
    # суммарная средняя прошедшая энергия в эрг:
    erg = sum(np.array(E_mean) * np.array(counts)) * 1.60217733e-6

    intensity = erg / (t_live * square)  # средняя интенсивность

    return intensity, counts


# Описываем функцию вычисления средней интенсивности за 1 день
def intensity_average_calculation_for_day(sat_id, track_day, filename):
    # Для начала получаем TLE
    tle_1, tle_2 = get_spacetrack_tle(sat_id, track_day, track_day + timedelta(days=1),
                                      USERNAME, PASSWORD)

    # Если не получилось добыть
    if not tle_1 or not tle_2:
        print('Impossible to retrieve TLE')
        return

    # Создаём экземляр класса Orbital
    orb = Orbital("N", line1=tle_1, line2=tle_2)

    minutes = 0  # счетчик минут
    n: int = 0  # количество обработанных строк
    ins: float = 0  # суммарное значение интенсивности за день
    sum_counts = [0] * 10  # суммарное количество отсчетов за день

    while minutes < 1436:

        utc_hour = int(minutes // 60)
        utc_minutes = int((minutes - (utc_hour * 60)) // 1)
        utc_seconds = int(round((minutes - (utc_hour * 60) - utc_minutes) * 60))

        utc_time = datetime(track_day.year, track_day.month, track_day.day,
                            utc_hour, utc_minutes, utc_seconds)

        # Вычисляем положение КА на небесной сфере в географических координатах
        lon, lat, alt = orb.get_lonlatalt(utc_time)

        # Находим интенсивности во время положения КА в экваториальной области
        if 0 < lat < lat_num[1]:
            for p in range(15):
                # Определяем номер строки в файле для нужных секунд дня
                k = 60 + minutes * 15 + p + 1
                # Находим значение интенсивности
                ii, cc = intensity_per_serving_time(filename, k)

                ins = ins + ii
                sum_counts = np.array(sum_counts) + np.array(cc)
                n += 1

        minutes += 1

    # Средняя интенсивность за день
    i_average = ins / n

    # Среднее число отсчетов за день
    c_average = sum_counts / n
    c_average_int = [int(round(q)) for q in c_average]

    return i_average, c_average_int


sum_intensity_average = 0
sum_counts_average = [0] * 10

for day in days:

    i, c = intensity_average_calculation_for_day(identifier, date(2009, 3, day),
                                                 names[day - 1])

    sum_intensity_average = sum_intensity_average + i
    sum_counts_average = np.array(sum_counts_average) + np.array(c)
    print('Средняя интенсивность за ' + str(day) + ' день:\n', i)
    print('Среднее число отсчетов в каждом из энергетических окон детектора за '
          + str(day) + ' день:\n', c)

intensity_average = sum_intensity_average / len(days)
counts_average = sum_counts_average / len(days)
counts_average_int = [int(round(y)) for y in counts_average]
print('\nСредняя интенсивность в результате обработки данных '
      + str(len(days)) + '-ти дней месяца:\n', intensity_average)
print('Среднее число отсчетов в каждом из энергетических окон детктора:\n',
      counts_average_int)

# Построение столбчатой диаграммы
counts_average_int.pop(0)
counts_average_int.pop()
channels = np.array([2, 3, 4, 5, 6, 7, 8, 9])
fig, ax = plt.subplots()
ax.bar(channels, counts_average_int, color='black')
plt.xlabel('Номер окна')
plt.ylabel('Количество отсчетов')
plt.title('Среднее число отсчетов в каждом из энергетических окон')
plt.show()

