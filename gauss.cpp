#include <math.h>
#include "gauss.h"

namespace {

	const double EPS = 1e-12;

	bool equal(double a, double b)
	{
		return fabs(b-a) <= EPS;
	}

	void swap_double(double &a, double &b)
	{
		double tmp = a;
		a = b;
		b = tmp;
	}

	std::vector<double> add(std::vector<double> a, std::vector<double> b)
	{
		std::vector<double> c;
		c.resize(a.size());
		for (int i = 0; i < c.size(); ++i)
			c[i] = a[i] + b[i];
		return c;
	}

	std::vector<double> mult(std::vector<double> a, double k)
	{
		std::vector<double> c;
		c.resize(a.size());
		for (int i = 0; i < c.size(); ++i)
			c[i] = a[i]*k;
		return c;
	}

}

std::vector<std::vector<double> > Gauss(
	std::vector<std::vector<double> > a, std::vector<double> b)
{
	std::vector<std::vector<double> > ans;
	int m, n;

	if (a.empty())
		throw;
	if (a[0].empty())
		throw;
	if (a.size() != b.size())
		throw;
	m = a.size();
	n = a[0].size();
	for (int i = 1; i < m; ++i) {
		if (a[i].size() != n)
			throw;
	}

	/* Прямой ход метода Гаусса */

	std::vector<int> jj;
	int j = 0;
	int r = 0;
	for (int i = 0; i < m; ++i) {
		double max_abs;
		int k_max;
		while (j <= n-1) {
			//находим максимальный по модулю элемент
			max_abs = 0;
			k_max = i;
			for (int k = i; k < m; ++k) {
				if (fabs(a[k][j]) > max_abs) {
					max_abs = fabs(a[k][j]);
					k_max = k;
				}
			}
			if (!equal(max_abs, 0))
				break;
			++j;
		}
		if (j > n-1)
			break;
		++r;
		jj.push_back(j);
		if (k_max != i) {
			//поменять местами i-ю строчку с k_max-й
			for (int l = j; l < n; ++l)
				swap_double(a[i][l], a[k_max][l]);
			swap_double(b[i], b[k_max]);
		}
		//делим i-е уравнение на a[i][j]
		for (int l = j+1; l < n; ++l) {
			a[i][l] /= a[i][j];
		}
		b[i] /= a[i][j];
		a[i][j] = 1;
		//путём элементарных преобразований
		//обнулить элементы a[i+1][j], a[i+2][j], ..., a[m-1][j]
		for (int k = i+1; k < m; ++k) {
			//умножаем i-е уравнение на a[k][j] и вычитаем из k-го уравнения
			for (int l = j+1; l < n; ++l)
				a[k][l] -= a[i][l]*a[k][j];
			b[k] -= b[i]*a[k][j];
			a[k][j] = 0;
		}
		++j;
	}

	/* Проверка системы на совместность */

	bool flag = true;
	for (int i = r; i < m; ++i) {
		if (!equal(b[i], 0)) {
			flag = false;
			break;
		}
	}
	if (!flag) {
		//система несовместна
		return ans;
	}

	/* Определение свободных переменных */

	int free_vars_count = n-r;
	ans.resize(n);
	for (int j = 0; j < n; ++j)
		ans[j].resize(free_vars_count + 1);
	if (r == 0) {
		for (int j = 0; j < n; ++j)
			ans[j][j+1] = 1;
	}
	else {
		int c = 0;
		for (int j = 0; j < jj[0]; ++j) {
			++c;
			ans[j][c] = 1;
		}
		for (int i = 0; i < r-1; ++i) {
			for (int j = jj[i]+1; j < jj[i+1]; ++j) {
				++c;
				ans[j][c] = 1;
			}
		}
		for (int j = jj[r-1]+1; j < n; ++j) {
			++c;
			ans[j][c] = 1;
		}
	}

	/* Обратный ход метода Гаусса */

	for (int i = r-1; i >= 0; --i) {
		ans[jj[i]][0] = b[i];
		for (j = jj[i]+1; j < n; ++j)
			ans[jj[i]] = add(ans[jj[i]], mult(ans[j], -a[i][j]));
	}

	return ans;
}
