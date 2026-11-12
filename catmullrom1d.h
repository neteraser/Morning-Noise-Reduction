
__inline double CatmullRom1D(double  t, double  p0, double  p1, double  p2, double  p3)
{
	double  t2 = t * t;
	double  t3 = t2 * t;
	return (0.5 * ((2.0 * p1) + (-p0 + p2) * t +
		(2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2 +
		(-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3));
}

__inline double CatmullRomDerivative1D(double  t, double  p0, double  p1, double  p2, double  p3) {
	double  t2 = t * t;

	return 0.5 * (
		(-p0 + p2) +
		2.0 * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t +
		3.0 * (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t2
		);
}

#define CATMULLROM_EPSILON 0.000005
//1e-6
#define CATMULLROM_EPSILON_SMALL 0.0000001
#define CATMULLROM_MAX_ITERATIONS 13
#define CATMULLROM_MAX_ITERATION_IN_ARRAY 6

__inline double InverseCatmullRom1D(double target_x, double p0, double p1, double p2, double p3)
{
	// Проверяем базовые границы (предполагаем, что функция монотонна на интервале)
	// Если значения на краях совпадают с целевым, возвращаем границы интервала
	if (fabs(p1 - target_x) < CATMULLROM_EPSILON) return 0.0;
	if (fabs(p2 - target_x) < CATMULLROM_EPSILON) return 1.0;

	// Инициализация границ для бинарного поиска (на случай сбоя метода Ньютона)
	double lower_t = 0.0;
	double upper_t = 1.0;

	// Начальное приближение — линейная интерполяция
	double t = (target_x - p1) / (p2 - p1);
	if (t < 0.0 || t > 1.0) {
		t = 1.0; // Защита от деления на ноль или выхода за границы
	}

	for (int i = 0; i < CATMULLROM_MAX_ITERATIONS; i++) {
		double current_x = CatmullRom1D(t, p0, p1, p2, p3);
		double diff = current_x - target_x;

		// Если точность достигнута, выходим
		if (fabs(diff) < CATMULLROM_EPSILON) {
			return t;
		}

		// Сужаем границы бинарного поиска в зависимости от монотонности
		if (p2 > p1) {
			if (diff > 0.0) upper_t = t;
			else lower_t = t;
		}
		else {
			if (diff > 0.0) lower_t = t;
			else upper_t = t;
		}

		// Считаем производную
		double derivative = CatmullRomDerivative1D(t, p0, p1, p2, p3);

		// Попытка сделать шаг методом Ньютона-Рафсона
		if (fabs(derivative) > 1e-9) {
			double next_t = t - diff / derivative;

			// Если шаг Ньютона не выходит за надежные границы бинарного поиска, принимаем его
			if (next_t > lower_t && next_t < upper_t) {
				t = next_t;
				continue;
			}
		}

		// Если метод Ньютона не сработал (вышел за границы или производная около 0),
		// делаем шаг бинарного поиска (Bisection) для стабильности
		t = 0.5 * (lower_t + upper_t);
	}

	return t; // Возвращаем наилучшее приближение после MAX_ITERATIONS
}

__inline int InverseCatmullRom1DInArray(int arraylength, double y_target, double p0, double p1, double p2, double p3) {

	if (fabs(p1 - y_target) < CATMULLROM_EPSILON) return 0;
	if (fabs(p2 - y_target) < CATMULLROM_EPSILON) return arraylength - 1;

	int i = 0;

	double t = fabs((y_target - p1) / (p2 - p1));
	if (t < 0.0 || t > 1.0)
		t = 0.5;

	double lower_t = 0.0;
	double upper_t = 1.0;

	double pt = -1.0;

	const double step = 0.25 / (double)(arraylength);

	for (int j = 0; j < CATMULLROM_MAX_ITERATION_IN_ARRAY; ++j)
	{
		t = fabs(t);

		if (std::abs(pt - t) < step)
			break;

		pt = t;

		double y_current = CatmullRom1D(t, p0, p1, p2, p3);
		double diff = y_current - y_target;

		if (std::abs(diff) > CATMULLROM_EPSILON)
		{
			if (p2 > p1) {
				if (diff > 0.0) upper_t = t;
				else lower_t = t;
			}
			else {
				if (diff > 0.0) lower_t = t;
				else upper_t = t;
			}

			double dy = CatmullRomDerivative1D(t, p0, p1, p2, p3);

			if (std::abs(dy) > CATMULLROM_EPSILON_SMALL)
			{
				double next_t = t - diff / dy;

				// Если шаг Ньютона не выходит за надежные границы бинарного поиска, принимаем его
				if (next_t > lower_t && next_t < upper_t) {
					t = next_t;
					continue;
				}
			}
		}
		else
		{
			break;
		}

		t = 0.5 * std::abs(lower_t + upper_t);

		if (t < 0.0)
			t = 0.0;
		if (t > 1.0)
			t = 1.0;

	}

	i = (int)truncf(fabs(t) * (double)(arraylength - 1));

	if (i > arraylength - 1)
		i = arraylength - 1;
	if (i < 0)
		i = 0;
	return i;
}
