# Оценка интенсивности фонового излучения на орбите КА в экваториальной области

Программа, которая по данным временных историй за март 2009 года оценивает среднюю интенсивность излучения на орбите КА КОРОНАС-ФОТОН в экваториальной области.

## Введение

В эксперименте Конус-РФ используются два детектора, один из которых всегда направлен в солнечном, а второй – всегда в антисолнечном направлении, Д1 и Д2 соответственно. Используя данные, получаемые детектором Д2, возможно оценить фоновую интенсивность излучения на орбите.

Оба детектора измеряют излучение в двух энергетических диапазонах: 10 кэВ...1 МэВ и 280 кэВ...10 МэВ. В первом энергетическом диапазоне измерение происходит в двенадцати энергетических окнах с временем накопления информации 1 сек, во втором диапазоне – в десяти энергетических интервалах с временем накопления информации 4 сек. Было принято решение обрабатывать данные второго энергетического диапазона по двум причинам:

- Во-первых, о распределении фотонов при энергиях последнего энергетического окна ничего не известно, поэтому предполагаем для оценки интенсивности снизу, что все фотоны в этом окне имели энергию равную нижней границе интервала. Тогда становится понятным, что выгоднее взять более широкий диапазон для того, чтобы как можно ближе подвести оценку снизу к реальной интенсивности фонового излучения на орбите КА;
  
- Во-вторых, из соображений о скорости выполнения программы, поскольку файлы временных историй второго диапазона потенциально содержат в 4 раза меньше строчек информации и имеют меньшее количество энергетических окон.
Также в приполярных областях и при прохождении через бразильскую аномалию влияние радиации на детектор сильно повышается, происходит насыщение. Таким образом, нам интересны только значения в северном полушарии и в областях неподалёку от экватора. Был выбран диапазон широт от 0° до 30°-45° с. ш.

Чтобы оценить интенсивность, нам потребуются значения средних каждого из энергетических окон. Так как с ростом энергии фотона их количество сокращается, то удобно брать среднее геометрическое левой и правой границ ширины энергетического окна. Однако, как было сказано ранее, мы не сможем найти среднее последнего энергетического окна, и поскольку вклад в общую прошедшую энергию у этого окна больше, чем у других окон, мы возьмем его нижнюю границу для того, чтобы приблизиться к получению оценки интенсивности снизу. Таким образом, для оценки будем брать средние геометрические всех энергетических окон, кроме последнего, у которого указана лишь нижняя граница (см. Таблица 1).




