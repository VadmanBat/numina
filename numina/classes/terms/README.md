# Классы Term — временные базисные функции (оригинальная версия)

**Дата документа:** 28 апреля 2026  
**Автор:** на основе кода Vadim (05-07.01.2025)  
**Статус:** Исходная реализация (до модернизации)

---

## 1. Назначение и контекст

Группа классов `Term<Type>` (и `TermExpression<Type>`) предназначена для **представления слагаемых во временной области**, которые возникают при:

- Обратном преобразовании Лапласа (inverse Laplace transform)
- Решении линейных ОДУ с постоянными коэффициентами
- Анализе переходных процессов в системах управления (transfer-function)

Каждый `Term` — это **одно слагаемое** вида:

```
c · t^n · e^{α t} · cos(ω t + φ)     или     c · t^n · e^{r t}
```

где параметры зависят от корней характеристического уравнения (или полюсов передаточной функции).

**Ключевые особенности оригинальной реализации:**
- Полиморфизм через абстрактный базовый класс `Term<Type>`
- **Глобальное статическое время** `Term<Type>::time` (устанавливается перед вычислением)
- Ручное управление памятью (`new` / `delete` в `TermExpression`)
- Автоматический выбор конкретного класса в `TermExpression::emplace_back(c, r, n)`
- HTML-разметка в строковом представлении (`<sup>`)

---

## 2. Базовый класс `Term<Type>`

```cpp
template <typename Type>
class Term {
    virtual Type value() const = 0;                    // значение y(t) при текущем Term::time
    virtual std::vector<Term*> derivative() const = 0; // слагаемые производной
    virtual Term* clone() const = 0;
    virtual bool isPositive() const = 0;
    virtual std::string string() const = 0;            // с HTML <sup>
    virtual std::string unsignedString() const = 0;
    virtual Type derivativeConstant() const;           // добавка к константе при дифференцировании
    static Type time;                                  // ГЛОБАЛЬНОЕ время!
};
```

**Математическая роль:**  
Каждый потомок реализует **одно слагаемое** разложения.

---

## 3. Классификация всех терминов (по типу корня и степени)

| Класс              | Power n | Тип корня `r`          | Математическая формула                           | Когда используется |
|--------------------|---------|------------------------|--------------------------------------------------|--------------------|
| (константа)        | 0       | `r = 0`                | $ c $                                            | `init_value`       |
| `TimeTerm`         | 1       | `r = 0`                | $$c \cdot t$$                                    | ramp               |
| `UTimeTerm`        | ≥ 2     | `r = 0`                | $$ c \cdot t^{n} $$                              | полиномиальные     |
| `ExpTerm`          | 0       | `r` вещественное ≠ 0   | $$ c \cdot e^{r t} $$                            | экспонента         |
| `ExpTimeTerm`      | 1       | `r` вещественное ≠ 0   | $$ c \cdot t \cdot e^{r t} $$                    | эксп. ramp         |
| `ExpUTimeTerm`     | ≥ 2     | `r` вещественное ≠ 0   | $$ c \cdot t^{n} \cdot e^{r t} $$                | эксп. полином      |
| `CosTerm`          | 0       | `r = jω` (чисто мним.) | $$ A \cos(\omega t + \phi) $$                    | гармоника          |
| `CosTimeTerm`      | 1       | `r = jω`               | $$ A \cdot t \cdot \cos(\omega t + \phi) $$      | гармоника × t      |
| `CosUTimeTerm`     | ≥ 2     | `r = jω`               | $$ A \cdot t^{n} \cdot \cos(\omega t + \phi) $$  | гармоника × tⁿ     |
| `ExpCosTerm`       | 0       | `r = α ± jω`           | $$ A e^{\alpha t} \cos(\omega t + \phi) $$       | затухающие колеб.  |
| `ExpCosTimeTerm`   | 1       | `r = α ± jω`           | $$ A t e^{\alpha t} \cos(\omega t + \phi) $$     | затух. × t         |
| `ExpCosUTimeTerm`  | ≥ 2     | `r = α ± jω`           | $$ A t^{n} e^{\alpha t} \cos(\omega t + \phi) $$ | затух. × tⁿ        |

**Примечания по выбору класса (из `emplace_back`):**
- Если `r == 0` → полиномиальная ветка (`TimeTerm`/`UTimeTerm`)
- Если `r.real() == 0 && r.imag() != 0` → чисто гармоническая (`Cos*`)
- Если `r.imag() == 0 && r.real() != 0` → чисто экспоненциальная (`Exp*`)
- Иначе → комплексная (`ExpCos*`)
- Для `n == 0` и `r == 0` → константа в `init_value`

---

## 4. Детальное описание каждого класса

### 4.1. TimeTerm (линейный член, n=1, r=0)

**Формула:**
$$ y(t) = c \cdot t $$

**Конструктор:**
```cpp
explicit TimeTerm(Type c);
```

**value():**
```cpp
return coefficient * Term<Type>::time;
```

**derivative():**
```cpp
return {};   // производная — чистая константа
```
**derivativeConstant():** `return coefficient;` ← добавляется в `init_value` производной.

**string():** `"3.14 × t"`

---

### 4.2. UTimeTerm (полиномиальный член, n ≥ 2, r=0)

**Формула:**
$$ y(t) = c \cdot t^{n} \quad (n \ge 2) $$

**Конструктор:**
```cpp
explicit UTimeTerm(Type a, int n);
```

**value():**
```cpp
return coefficient * std::pow(Term<Type>::time, power);
```

**derivative():**
```cpp
if (power == 2)
    return {new TimeTerm(power * coefficient)};   // 2c · t
else
    return {new UTimeTerm(power * coefficient, power-1)};
```

**derivativeConstant():** `return coefficient;` (возможно, **ошибка** — для n≥2 не должно быть постоянного слагаемого в производной).

**string():** `"2.5 × t<sup>3</sup>"`

---

### 4.3. ExpTerm (экспоненциальный член, n=0)

**Формула:**
$$ y(t) = c \cdot e^{r t} $$

**Конструктор:**
```cpp
explicit ExpTerm(Type c, Type r);
```

**value():**
```cpp
return coefficient * std::exp(root * Term<Type>::time);
```

**derivative():**
```cpp
return {new ExpTerm(coefficient * root, root)};   // c·r · e^{r t}
```

**string():** `"1.2 × e<sup>-0.5 × t</sup>"`

---

### 4.4. ExpTimeTerm (экспонента × время, n=1)

**Формула:**
$$ y(t) = c \cdot t \cdot e^{r t} $$

**Конструктор:**
```cpp
explicit ExpTimeTerm(Type c, Type r);
```

**value():**
```cpp
return coefficient * std::exp(root * Term<Type>::time) * Term<Type>::time;
```

**derivative():**
```cpp
return {
    new ExpTerm(coefficient * root, root),   // от экспоненты
    new ExpTerm(coefficient, root)           // от t (даёт чистую экспоненту)
};
```

**string():** `"0.8 × e<sup>1.5 × t</sup> × t"`

---

### 4.5. ExpUTimeTerm (экспонента × полином, n ≥ 2)

**Формула:**
$$ y(t) = c \cdot t^{n} \cdot e^{r t} $$

**derivative() (для n==2):**
```cpp
return {
    new ExpTimeTerm(2 * coefficient, root),           // 2c t e^{rt}
    new ExpUTimeTerm(root * coefficient, root, 2)     // c r t² e^{rt}
};
```

**Общий случай (n > 2):**
- `power * coefficient` → понижение степени
- `root * coefficient` → та же степень (от дифференцирования экспоненты)

---

### 4.6. CosTerm (гармонический член, n=0)

**Формула:**
$$ y(t) = A \cos(\omega t + \phi) $$
где $$ A = 2 |c| $$, $$ \omega = \Im(r) $$, $$ \phi = \arg(c) $$

**Конструкторы:**
- Из комплексной пары: `CosTerm(const Comp& c, const Comp& r)`
- Прямой: `CosTerm(Type A, Type ω, Type φ)`

**value():**
```cpp
return amplitude * std::cos(omega * Term<Type>::time + phi);
```

**derivative():**
```cpp
return {new CosTerm(-amplitude * omega, omega, phi - π/2)};
```
(сдвиг фазы на -90° превращает cos в -sin)

**string():** `"5.0 × cos(2.0 × t + 0.3)"` или с `<sup>` нет (только для t)

---

### 4.7. CosTimeTerm (гармоника × t, n=1)

**Формула:**
$$ y(t) = A \cdot t \cdot \cos(\omega t + \phi) $$

**derivative():**
```cpp
return {
    new CosTerm(amplitude, omega, phi),                    // от t
    new CosTimeTerm(-amplitude * omega, omega, phi - π/2)  // от cos
};
```

---

### 4.8. CosUTimeTerm (гармоника × tⁿ, n ≥ 2)

Аналогично, но с `std::pow(t, power)` и рекурсивным понижением степени + добавление сдвинутой по фазе версии той же степени.

---

### 4.9. ExpCosTerm (затухающие/растущие колебания, n=0)

**Формула:**
$$ y(t) = A \cdot e^{\alpha t} \cdot \cos(\omega t + \phi) $$
где $$ A = 2 |c| $$, $$ \alpha = \Re(r) $$, $$ \omega = \Im(r) $$, $$ \phi = \arg(c) $$

**derivative():**
```cpp
return {
    new ExpCosTerm(amplitude * alpha, alpha, omega, phi),           // от экспоненты
    new ExpCosTerm(-omega * amplitude, alpha, omega, phi - π/2)     // от косинуса
};
```

---

### 4.10. ExpCosTimeTerm (затух. колебания × t, n=1)

**Формула:**
$$ y(t) = A \cdot t \cdot e^{\alpha t} \cdot \cos(\omega t + \phi) $$

**derivative() (оригинал):**
```cpp
return {
    new ExpCosTerm(amplitude, alpha, omega, phi),
    new ExpCosTimeTerm(amplitude, alpha, omega, phi),           // ← здесь пропущен множитель alpha!
    new ExpCosTimeTerm(-amplitude * omega, alpha, omega, phi - π/2)
};
```

> **Замечание:** Во второй строке должно быть `amplitude * alpha`, иначе математически неверно. Это известная ошибка оригинала.

---

### 4.11. ExpCosUTimeTerm (затух. колебания × tⁿ, n ≥ 2)

Самый сложный. Для `n == 2` возвращает **три** слагаемых:
1. `ExpCosTimeTerm(2 * amplitude, ...)` — от степени
2. `ExpCosUTimeTerm(alpha * amplitude, ...)` — от экспоненты
3. `ExpCosUTimeTerm(-amplitude * omega, ...)` — от косинуса (фазовый сдвиг)

Для `n > 2` — три слагаемых с понижением степени + alpha-компонента + omega-компонента.

---

## 5. TermExpression — «сборщик» выражения

```cpp
TermExpression(const std::vector<std::complex<Type>>& roots,
               const std::vector<std::complex<Type>>& coeffs,
               const std::vector<int>& powers);
```

**Логика выбора класса** (см. `emplace_back`):
- `r == 0` → `TimeTerm` / `UTimeTerm` / константа
- `r.real() == 0` → `Cos*`
- `r.imag() == 0` → `Exp*`
- иначе → `ExpCos*`

**Вычисление:**
```cpp
Type operator()(double t) const {
    Term<Type>::time = t;           // ГЛОБАЛЬНО!
    Type result = init_value;
    for (auto* term : terms) result += term->value();
    return result;
}
```

**Производная выражения:**
```cpp
TermExpression derivative() const {
    Type d_init = 0;
    std::vector<Term*> d_terms;
    for (auto* term : terms) {
        d_init += term->derivativeConstant();
        auto d = term->derivative();
        d_terms.insert(d_terms.end(), d.begin(), d.end());
    }
    return TermExpression(d_terms, d_init);
}
```

---

## 6. Известные особенности и потенциальные проблемы оригинальной версии

1. **Глобальное состояние** `Term::time` — не потокобезопасно, не функционально.
2. **Ручное управление памятью** — `TermExpression` владеет сырыми указателями, `delete` в деструкторе и assignment. Легко получить утечки/двойное удаление.
3. **Отсутствие copy-конструктора** у `TermExpression` (только `operator=`).
4. **Математическая неточность** в `ExpCosTimeTerm::derivative()` — пропущен множитель `alpha`.
5. **Подозрительный `derivativeConstant()` в `UTimeTerm`** — возвращает `coefficient` вместо `0`.
6. **HTML в string()** — удобно для веб, но не для логов/файлов.
7. **Специальные case'ы для power==2** — технический долг (чтобы не создавать `UTimeTerm` с power=1, которого нет).

---

## 7. Пример использования (псевдокод)

```cpp
std::vector<std::complex<double>> roots = {0.0, -1.5, std::complex<double>(-0.2, 3.0)};
std::vector<std::complex<double>> coeffs = {2.0, 1.0, std::complex<double>(0.5, 0.1)};
std::vector<int> powers = {0, 0, 1};

TermExpression<double> expr(roots, coeffs, powers);

double y = expr(0.5);           // y(0.5)
auto expr2 = expr.derivative(); // y'(t)
std::cout << expr.string();     // "2 + 1 × e<sup>-1.5 × t</sup> - 1.02 × e<sup>-0.2 × t</sup> × cos(3 × t - 0.197) × t"
```

---

## 8. Следующие шаги (план модернизации)

1. Убрать `static time` → `value(Type t)` и `operator()(Type t)`
2. Заменить сырые указатели на `std::unique_ptr<Term<Type>>`
3. Добавить полноценный `copy ctor`
4. Исправить математические ошибки в производных
5. Сделать строки без HTML (опционально `toHtmlString()`)
6. Привести стиль к `Polynomial` / `PolySolver` (move semantics, `[[nodiscard]]`, `noexcept`, `Type` алиасы и т.д.)
