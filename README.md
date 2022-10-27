# gauss_block_method
## ТРЕБОВАНИЯ К ПРОГРАММАМ

1. Программа должна получать все параметры в качестве аргументов командной строки. Аргу-
    менты командной строки:

```
1) n– размерность матрицы,
2) m– размер блока,
3) r– количество выводимых значений в матрице,
4) s– задает номер формулы для инициализации матрицы, должен быть равен 0 при вводе
матрицы из файла
5) filename– имя файла, откуда надо прочитать матрицу. Этот аргумент отсутствует ,
еслиs!= 0.
```
```
Например, запуск
```
```
./a.out 4 3 4 0 a.txt
```
```
означает, что матрицу 4x4 надо прочитать из файла a.txt, использовать блочный алгоритм
с размером блока 3 и выводить не более 4-х строк и столбцов матрицы, а запуск
```
```
./a.out 2000 90 6 1
```
```
означает, что матрицу2000x2000надо инициализировать по формуле номер 1, использовать
блочный алгоритм с размером блока 90 и выводить не более 6-ти строк и столбцов матрицы.
```
2. В задачах, где требуется правая частьb, этот вектор строится после инициализации матрицы
    A= (ai,j)i,j= 1 ,...,n по формуле:
3. Ввод матрицы должен быть оформлен в виде подпрограммы, находящейся в отдельном файле.
4. Ввод матрицы из файла. В указанном файле находится матрица в формате:

```
a 1 , 1 ... a 1 ,n
a 2 , 1 ... a 2 ,n
.........
an, 1 ... an,n
гдеn- указанный размер матрицы,A= (ai,j)- матрица. Программа должна выводить сообще-
ние об ошибке, если указанный файл не может быть прочитан, содержит меньшее количество
данных или данные неверного формата.
```
5. Ввод матрицы и правой части по формуле. Элементai,j матрицыAполагается равным

```
ai,j=f(s,n,i,j), i,j= 1 ,...,n,
где f(s,n,i,j)- функция, которая возвращает значение (i,j)-го элемента n×nматрицы по
формуле номерs(аргумент командной строки). Функция f(s,n,i,j)должна быть оформлена
в виде отдельной подпрограммы.
```
```
```
n−max{i,j}+ 1 при s= 1
max{i,j} при s= 2
|i−j| при s= 3
1
i+j− 1
```
```
при s= 4
```

6. Решение системы (нахождение обратной матрицы) должно быть оформлено в виде функции,
    находящейся в отдельном файле и получающей в качестве аргументов

```
(a) размерность n матрицы A,
(b) размер блока m,
(c) матрицу A,
(d) правую часть b (если стоит задача решить линейную систему)
(e) векторx, в который будет помещено решение системы, если стоит задача решить линей-
ную систему, или матрицу X, в которую будет помещена обратная матрица, если стоит
задача обратить матрицу,
(f) дополнительные вектора, если алгоритму требуется дополнительная память.
```
```
Получать в этой подпрограмме дополнительную информацию извне через глобальные пере-
менные и т.п. запрещается.
```
7. Функция, реализующая задачу, возвращает ненулевое значение, если алгоритм решения
    неприменим к поданной на вход матрице A.
8. Функция, реализующая задачу, **не должна выделять или использовать дополнительную**
    **память**.
9. **Сложность работы функции** , реализующая задачу, не должна превышатьO(n^3 ).
10. Суммарный объем оперативной памяти, требуемой программе, не должен превышать:
- при вычислении решения системы:n^2 +O(n),
- при вычислении обратной матрицы: 2 n^2 +O(n).

```
Для выполнения этого требования после завершения алгоритма решения (нахождения обрат-
ной матрицы) вызывается подпрограмма инициализации матрицы (из файла или по формуле)
и вычисления вектораb(в задачах решения линейной системы).
```
11. Программа должна содержать подпрограмму вывода на экран прямоугольной матрицы l×n
    матрицы. Эта подпрограмма используется для вывода исходнойn×nматрицы после ее ини-
    циализации, а также для вывода на экран решения системы ( 1 ×n матрицы) или обратной
    n×nматрицы, если стоит задача обратить матрицу. Подпрограмма выводит на экран не бо-
    лее, чемrстрок и столбцовl×nматрицы, где r– параметр этой подпрограммы (аргумент
    командной строки). Каждая строка матрицы должна печататься на новой строке, каждый эле-
    мент матрицы выводится в строке по формату" %10.3e"(один пробел между элементами
    и экспоненциальный формат%10.3e).
12. Результатами работы программы являются 3 элемента:
    - Собственно вектор решенияx (в задачах нахождения решения линейной системы) или
       обратная матрицаA−^1 (в задачах нахождения обратной матрицы).
    - Два вещественных числаr 1 и r 2 , вычисляемых после вызова функции, реализующей
       задачу:
**- В задачах нахождения обратной матрицы**

Вычислениеr r1 и r2 должно быть оформлено в виде подпрограммы, вызываемой из функции
main. Эта подпрограмма не должна выделять или использовать дополнительную память.
```
13. Вывод результата работы функции в функцииmainдолжен производиться по формату:
    - Непосредственно вывод вектора решенияxили обратной матрицыA−^1 :
       **-** в задачах нахождения решения линейной системы вывод вектора решенияxпроиз-
          водится вызовом подпрограммы печати матрицы (см. пункт 11) размера 1 ×n(т.е. в
          строку и **по указанному там формату** )
       **-** в задачах нахождения обратной матрицы вывод обратной матрицы A−^1 произво-
          дится вызовом подпрограммы печати матрицы (см. пункт 11) размера n×n ( **по**
          **указанному там формату** )
Вывод не производится, если алгоритм решения неприменим к поданной на вход матри-
це A.
    - Отчет о результате и времени работы:
       printf (
       "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n",
       argv[0], task, r1, r2, t1, t2, s, n, m);

```
где
```
**-** argv[0]– первый аргумент командной строки (имя образа программы),
**-** task– номер задачи,
**-** r1=r 1 – вычисленное значениеr 1 (см. пункт 12), выводитсяr 1 =− 1 , если алго-
    ритм решения неприменим к поданной на вход матрице A,
**-** r2=r 2 – вычисленное значениеr 2 (см. пункт 12), выводитсяr 2 =− 1 , если алго-
    ритм решения неприменим к поданной на вход матрице A,
**-** t1– время работы функции, реализующей решение этой задачи, в секундах (с точ-
    ностью до сотых долей),
**-** t2– время работы функции, вычисляющей невязки решения (см. пункт 12), в се-
    кундах (с точностью до сотых долей), выводится t2= 0 , если алгоритм решения
    неприменим к поданной на вход матрице A,
**-** s,n,m– аргументы командной строки.

```
Вывод должен производиться в точности в таком формате , чтобы можно было автомати-
зировать обработку запуска многих тестов. Вывод отчета о результате и времени работы
должен производится всегда , даже если алгоритм решения неприменим к поданной на вход
матрице A.
```

## МЕТОД РЕШЕНИЯ СИСТЕМ ЛИНЕЙНЫХ УРАВНЕНИЙ
13. Метод Гаусса нахождения обратной матрицы с выбором главного элемента по строке.


