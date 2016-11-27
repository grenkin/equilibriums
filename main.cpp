#include <fstream>
#include <vector>
#include <math.h>
#include "gauss.h"

const double EPS = 1e-12;

std::ifstream fin("input.txt");
std::ofstream fout("output.txt");

std::vector<std::vector<double> > mu0, k, v0, qmin, qmax;
std::vector<double> z0;
int m, n;
double D;

std::vector<std::vector<double> > z_bar;
std::vector<bool> I_star, J_star;
std::vector<double> z_bar_max;

double equal(double a, double b)
{
	return fabs(b-a) <= EPS;
}

double less_or_equal(double a, double b)
{
	return a < b + EPS;
}

double q_tilde(int i, int j, double z)
{
	return qmax[i][j]*z/(z + (k[i][j] + z)/v0[i][j]*(qmax[i][j]-qmin[i][j])*D);
}

double q_tilde_inv(int i, int j, double q)
{
	return k[i][j]/v0[i][j]*(qmax[i][j]-qmin[i][j])*D*q/
		(qmax[i][j] - q*(1 + 1/v0[i][j]*(qmax[i][j]-qmin[i][j])*D));
}

double mu(int i, int j, double q)
{
	return mu0[i][j]*(1 - qmin[i][j]/q);
}

double mu_inv(int i, int j, double mu)
{
	return qmin[i][j]/(1 - mu/mu0[i][j]);
}

double mu_tilde_inv(int i, int j, double mu)
{
	return q_tilde_inv(i, j, mu_inv(i, j, mu));
}

/*bool check(std::vector<double> y, std::vector<double> z)
{
	for (int j = 0; j < n; ++j) {
		double sum = 0;
		for (int i = 0; i < m; ++i)
			sum += q_tilde(i, j, z[j])*y[i];
		if (!equal(sum, z0[j]-z[j])) {
			//fout << sum << " " << z0[j]-z[j]  << "\n";
			return false;
		}
	}
	return true;
}*/

void J_star_work()
{
	std::vector<int> vars;
	std::vector<int> vars_numbers;
	vars_numbers.resize(m);
	for (int i = 0; i < m; ++i) {
		if (I_star[i]) {
			vars.push_back(i);
			vars_numbers[i] = vars.size()-1;
		}
		else
			vars_numbers[i] = -1;
	}
	std::vector<int> equations;
	for (int j = 0; j < n; ++j) {
		if (J_star[j])
			equations.push_back(j);
	}
	std::vector<std::vector<double> > a;
	std::vector<double> b;
	a.resize(equations.size());
	b.resize(equations.size());
	for (int j = 0; j < equations.size(); ++j) {
		a[j].resize(vars.size());
		for (int i = 0; i < vars.size(); ++i)
			a[j][i] = q_tilde(vars[i], equations[j], z_bar_max[equations[j]]);
		b[j] = z0[equations[j]] - z_bar_max[equations[j]];
	}
	std::vector<std::vector<double> > yy = Gauss(a, b);
	if (!yy.empty()) {
		if (yy[0].size() > 1) {
			//решений бесконечно много
			//fout << "Решений бесконечно много\n";
			int free_vars_count = yy[0].size() - 1;
			fout << "y = (";
			for (int i = 0; i < m; ++i) {
				//fout << y[i];
				if (!I_star[i])
					fout << "0";
				else {
					//yy[vars_numbers[i]]
					fout << yy[vars_numbers[i]][0];
					for (int j = 1; j <= free_vars_count; ++j) {
						if (!equal(yy[vars_numbers[i]][j], 0)) {
							if (yy[vars_numbers[i]][j] > 0)
								fout << " + " << yy[vars_numbers[i]][j] << "a[" << j << "]";
							else
								fout << " - " << -yy[vars_numbers[i]][j] << "a[" << j << "]";
						}
					}
				}
				if (i != m-1)
					fout << ", ";
			}
			fout << ")\n";
			fout << "y[i] >= 0\n";
			for (int j = 0; j < n; ++j) {
				if (!J_star[j]) {
					bool flag = false;
					for (int i = 0; i < m; ++i) {
						if (I_star[i]) {
							if (flag)
								fout << " + ";
							flag = true;
							fout << q_tilde(i, j, z_bar_max[j]) << "y[" << i+1 << "]";
						}
					}
					fout << " <= ";
					fout << z0[j] - z_bar_max[j];
					fout << "\n";
				}
			}
			fout << "z = (";
			for (int j = 0; j < n; ++j) {
				if (J_star[j])
					fout << z_bar_max[j];
				else
					fout << "?";
				if (j != n-1)
					fout << ", ";
			}
			fout << ")\n";
			fout << "\n";
		}
		else {
			//решение единственно
			std::vector<double> y;
			y.resize(m);
			for (int i = 0; i < m; ++i) {
				if (!I_star[i])
					y[i] = 0;
				else
					y[i] = yy[vars_numbers[i]][0];
			}

			//проверяем неотрицательность y[i]
			bool flag = true;
			for (int i = 0; i < m; ++i) {
				if ((y[i] < 0) && !equal(y[i], 0)) {
					flag = false;
					break;
				}
			}
			if (!flag)
				return;

			std::vector<double> z;
			z.resize(n);
			for (int j = 0; j < n; ++j) {
				if (J_star[j])
					z[j] = z_bar_max[j];
				else {
					double l = 0, r = z0[j];
					while (r-l > 1e-10) {
						double q = (l+r)/2;
						double sum = 0;
						for (int i = 0; i < vars.size(); ++i)
							sum += q_tilde(vars[i], j, q)*y[vars[i]];
						if (sum > z0[j]-q)
							r = q;
						else
							l = q;
					}
					z[j] = (l+r)/2;
				}

			}

			//провяряем неравенства для z[j]
			flag = true;
			for (int j = 0; j < n; ++j) {
				if ((z[j] < z_bar_max[j]) && !equal(z[j], z_bar_max[j])) {
					flag = false;
					break;
				}
			}
			if (!flag)
				return;

			//fout << check(y, z) << "\n";

			fout << "y = (";
			for (int i = 0; i < m; ++i) {
				fout << y[i];
				if (i != m-1)
					fout << ", ";
			}
			fout << ")\n";
			fout << "z = (";
			for (int j = 0; j < n; ++j) {
				fout << z[j];
				if (j != n-1)
					fout << ", ";
			}
			fout << ")\n";

			bool stable = true;
			for (int i = 0; i < m; ++i) {
				bool flag = false;
				for (int j = 0; j < n; ++j) {
					double q = q_tilde(i, j, z[j]);
					if (less_or_equal(mu(i, j, q), D)) {
						flag = true;
						break;
					}
				}
				if (!flag) {
					stable = false;
					break;
				}
			}
			if (stable)
				fout << "Равновесие устойчиво\n";
			else
				fout << "Равновесие неустойчиво\n";
			fout << "\n";
		}
	}

/*	fout << "I_star: ";
	for (int i = 0; i < m; ++i)
		fout << I_star[i] << " ";
	fout << "\n";
	fout << "J_star: ";
	for (int j = 0; j < n; ++j)
		fout << J_star[j] << " ";
	fout << "\n";
	fout << "\n";*/
}

void rec_J_star(int j)
{
	if (j == n) {
		bool correct = true;
		for (int i = 0; i < m; ++i) {
			if (I_star[i]) {
				bool flag = false;
				for (int j = 0; j < n; ++j) {
					if (J_star[j]) {
						if (equal(z_bar[i][j], z_bar_max[j])) {
							flag = true;
							break;
						}
					}
				}
				if (!flag) {
					correct = false;
					break;
				}
			}
		}
		if (correct)
			J_star_work();
	}
	else {
		J_star[j] = false;
		rec_J_star(j+1);
		J_star[j] = true;
		rec_J_star(j+1);
	}
}

void I_star_work()
{
	bool flag = true;
	for (int i = 0; i < m; ++i) {
		if (I_star[i]) {
			flag = false;
			break;
		}
	}
	if (flag) {
		//все y[i] равны 0
		fout << "y = (";
		for (int i = 0; i < m; ++i) {
			fout << "0";
			if (i != m-1)
				fout << ", ";
		}
		fout << ")\n";
		fout << "z = (";
		for (int j = 0; j < n; ++j) {
			fout << z0[j];
			if (j != n-1)
				fout << ", ";
		}
		fout << ")\n";
		bool stable = true;
		for (int i = 0; i < m; ++i) {
			bool flag1 = false;
			for (int j = 0; j < n; ++j) {
				double q = q_tilde(i, j, z0[j]);
				if (less_or_equal(mu(i, j, q), D)) {
					flag1 = true;
					break;
				}
			}
			if (!flag1) {
				stable = false;
				break;
			}
		}
		if (stable)
			fout << "Равновесие устойчиво\n";
		else
			fout << "Равновесие неустойчиво\n";
		fout << "\n";
	}
	else {
		for (int j = 0; j < n; ++j) {
			z_bar_max[j] = 0;
			for (int i = 0; i < m; ++i) {
				if (I_star[i]) {
					if (z_bar[i][j] > z_bar_max[j])
						z_bar_max[j] = z_bar[i][j];
				}
			}
		}
		rec_J_star(0);
	}

/*	for (int i = 0; i < m; ++i)
		std::cout << I_star[i] << " ";
	std::cout << std::endl;*/
}

void rec_I_star(int i)
{
	if (i == m)
		I_star_work();
	else {
		I_star[i] = false;
		rec_I_star(i+1);
		I_star[i] = true;
		rec_I_star(i+1);
	}
}

int main()
{
	fin >> m >> n;
	if (m < 1) {
		fout << "m < 1";
		return 0;
	}
	if (n < 1) {
		fout << "n < 1";
		return 0;
	}
	z0.resize(n);
	mu0.resize(m);
	k.resize(m);
	v0.resize(m);
	qmin.resize(m);
	qmax.resize(m);
	for (int i = 0; i < m; ++i) {
		mu0[i].resize(n);
		k[i].resize(n);
		v0[i].resize(n);
		qmin[i].resize(n);
		qmax[i].resize(n);
	}
	fin >> D;
	for (int i = 0; i < n; ++i) {
		fin >> z0[i];
		if (z0[i] <= 0) {
			fout << "z0[" << i+1 << "] <= 0";
			return 0;
		}
	}
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			fin >> mu0[i][j];
			if (mu0[i][j] <= 0) {
				fout << "mu0[" << i+1 << "][" << j+1 << "] <= 0";
				return 0;
			}
			/*if (mu0[i][j] <= D) {
				fout << "mu0[" << i+1 << "][" << j+1 << "] <= D";
				return 0;
			}*/
		}
	}
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			fin >> k[i][j];
			if (k[i][j] <= 0) {
				fout << "k[" << i+1 << "][" << j+1 << "] <= 0";
				return 0;
			}
		}
	}
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			fin >> v0[i][j];
			if (v0[i][j] <= 0) {
				fout << "v0[" << i+1 << "][" << j+1 << "] <= 0";
				return 0;
			}
		}
	}
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			fin >> qmin[i][j];
			if (qmin[i][j] <= 0) {
				fout << "qmin[" << i+1 << "][" << j+1 << "] <= 0";
				return 0;
			}
		}
	}
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			fin >> qmax[i][j];
			if (qmax[i][j] <= 0) {
				fout << "qmax[" << i+1 << "][" << j+1 << "] <= 0";
				return 0;
			}
		}
	}

	z_bar.resize(m);
	for (int i = 0; i < m; ++i) {
		z_bar[i].resize(n);
		for (int j = 0; j < n; ++j) {
			z_bar[i][j] = mu_tilde_inv(i, j, D);
			if (z_bar[i][j] <= 0) {
				fout << "z_bar[" << i+1 << "][" << j+1 << "] <= 0";
				return 0;
			}
		}
	}

	I_star.resize(m);
	J_star.resize(n);
	z_bar_max.resize(n);
	rec_I_star(0);

	return 0;
}
