#include<iostream>
#include<vector>
#include<algorithm>
#include<chrono>
#include<fstream>
#include<cmath>
# define M_PI 3.14159265358979323846
//Исследование схемы Кранка Николсона для уравнения теплопроводности на точность
using namespace std;
class equation {
public:
	double u(double t, double x)
	{
		this->x = x;
		this->t = t;
		return (-2 * pow(x, 4) - 3 * pow(t, 3) + 3 * pow(t, 2) * x + exp(x));
	}
	double f(double t, double x) {
		this->x = x;
		this->t = t;
		return (-9 * pow(t, 2) + 6 * x * t + a * (24 * pow(x, 2) - exp(x)));
	}
	void showSolution(ofstream& fout, int k) {
		if (k == 1) {
			for (int i = 0; i <= n + 1; ++i) {
				fout<<(k - 1)* tau << " " << i * h << " " << u((k - 1) * tau, i * (h)) << endl;
			}
			fout << endl;
		}
		fout <<(k)* tau << " " << 0 * h << " " << u(k * tau, 0 * (h)) << endl;
		for (int i = 1; i <= n; i++) {
			fout << k* tau << " " << i * h << " " << sol[i - 1] << endl;
		}
		fout <<(k)* tau << " " <<(n + 1) * h << " " << u(k * tau, (n + 1) * (h)) << endl;
		fout << endl;
	}
	void fill_sol(int j) {
		sol.resize(n);
		for (int i = 1; i <= n; ++i) {
			sol[i - 1] = u(j * tau, i * h);
		}
	}
	void fill_bs(int j) // For the first layout
	{
		if (j == 0) {
			fill_sol(j);
			bs.resize(n);
		}

		for (int i = 1; i <= n; ++i) {
			if (i == 1)
				bs[i - 1] = (sol[i - 1] / tau) + a * (((u(j * tau, (i - 1) * h) - 2 * sol[i - 1] + sol[i]) / (2 * pow(h, 2)))) + f((j * tau) + (tau / 2), i * h) - u((j + 1) * tau, (i - 1) * h) * (koef2);
			else if (i == n)
				bs[i - 1] = (sol[i - 1] / tau) + a * (((sol[i - 2] - 2 * sol[i - 1] + u(j * tau, (i + 1) * h)) / (2 * pow(h, 2)))) + f((j * tau) + (tau / 2), i * h) - u((j + 1) * tau, (i + 1) * h) * (koef2);
			else {
				bs[i - 1] = (sol[i - 1] / tau) + a * (((sol[i - 2] - 2 * sol[i - 1] + sol[i]) / (2 * pow(h, 2)))) + f((j * tau) + (tau / 2), i * h);
			}
		}

	}


	void Tridiagonal_method(ofstream& fout){
	
		m = (int)(1 / tau);
		n = (int)(1 / h) - 1;
		for (int k = 1; k <= m; ++k) {
			fill_bs(k - 1);
			vector<double> y1(n);
			vector <double> a(n);
			vector <double> b(n);
			y1[0] = koef1;
			a[0] = -koef2 / y1[0];
			b[0] = bs[0] / y1[0];
			for (int i = 1; i < n - 1; i++) {
				y1[i] = koef1 + koef2 * a[i - 1];
				a[i] = -koef2 / y1[i];
				b[i] = (bs[i] - koef2 * b[i - 1]) / y1[i];

			}
			y1[n - 1] = koef1 + koef2 * a[n - 2];
			b[n - 1] = (bs[n - 1] - koef2 * b[n - 2]) / y1[n - 1];
			sol[n - 1] = (b[n - 1]);
			
			for (int i = (n - 2); i >= 0; i--) {
				sol[i] = (a[i] * sol[i + 1] + b[i]);
			}
			
			if (k == 1) {
				maxi = abs(sol[0] - u(k * tau, h));
			}
			
			for (int i = 1; i <= n; ++i) {
				if (abs(u(k * tau, i * h) - sol[i - 1]) > maxi) {
					maxi = abs(sol[i - 1] - u(k * tau, i * h));
				}
			}
			showSolution(fout, k);
		}
	}
private:
	double x, t, h, tau;
	double a = 0.019;
	vector<double> bs;// Right part for the Tridiagonal method
	int n, m;
	vector<double> sol;// Current solution for the Tridiagonal method
	long double maxi;
	double koef1, koef2;
	double sum;

public:
	double cmax() {
		return maxi;
	}

	void get_h_and_tau() {
		cout << " Enter h: ";
		cin >> h;
		cout << " Enter tau: ";
		cin >> tau;
		koef1 = (1 / tau) + (a / pow(h, 2));
		koef2 = -(a / (2 * pow(h, 2)));

	}
	void get_exact_solution(ofstream& fout) {

		for (int j = 0; j <= m; j++) {

			for (int i = 0; i <= n + 1; i++) {
				fout<<j* tau<<" "<<i * h <<" " << u(j * tau, i * (h))<< endl;
			}
			fout << endl;
		}
		cout << " ======== " << endl;
		cout <<cmax() << endl;

	}
	
};

int main() {
	
	ofstream fout;
	fout.open("C:/Users/Ekaterina/Desktop/output2.txt");
	ofstream fout1;
	fout1.open("C:/Users/Ekaterina/Desktop/output1.txt");
	equation e;
	e.get_h_and_tau();
	auto start = chrono::high_resolution_clock::now();
	e.Tridiagonal_method(fout1);

	e.get_exact_solution(fout);
	auto end = chrono::high_resolution_clock::now();
	chrono::duration<float> dur = end - start;
	cout << "Our time: " << dur.count() << endl;
	fout.close();
	system("pause");
	return 0;
}