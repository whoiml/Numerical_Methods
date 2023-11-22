#include "LSM.h"

int main() {
	vector <double> x = { 1.1, 1.4, 1.7, 2.1, 2.6, 4.7, 6.1, 7.0, 10.0, 12.8, 16.5, 20.8, 40.6 };
	vector <double> y = { 25.0, 22.7, 22.1, 19.8, 17.0, 12.3, 10.7, 10.0, 8.2, 6.7, 5.6, 5.0, 3.5 };

	// y = ax^(-1/e) - моя функция 
	/*	Вариаций того, что может быть в этих скобках - много,
		поэтому примем тут в main, что конкретно для нас (-1/e) = b
		далее просто вычислим 'e' из 'b'  */

	// Решение для y=ax^b
	vector<double> parameters = solvePowerFunctionApproximation(x, y);

	cout << "a = " << parameters[0] << "\te = " << -1 / parameters[1] << "\n";
	cout << "\ny = (" << parameters[0] << ")x^(-1/" << -1 / parameters[1] << ")\n";
	return 0;
}
