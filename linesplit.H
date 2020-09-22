#ifndef LINESPLIT_H
#define LINESPLIT_H

namespace linesplit
{
    template <class T>
    void linear(double f_int, T start, T end)
    {
        double len = *(end - 1) - *start;
        int n = end - start;
        double a = 2 * len / (n - 2) / (n - 1) - 2 * f_int / (n - 2);
        for (int i = 0;i < n - 2;i++)
        {
            *(start + i + 1) = *(start + i) + a * i + f_int;
        }
    }
    template <class T>
    void linear(T start, T end, double f_int)
    {
        double len = *(end - 1) - *start;
        int n = end - start;
        double a = 2 * len / (n - 2) / (n - 1) - 2 * f_int / (n - 2);
        for (int i = 0;i < n - 2;i++)
        {
            *(end - i - 2) = *(end - i - 1) - a * i - f_int;
        }
    }
    template <class T>
    void linear(T start, T end)
    {
        linear(start, end, (*(end - 1) - *start) / (end - start - 1));
    }

    double solveK(int n, double df, double s)
    {
        double k = 2;
        double f;
        double ff;
        do
        {
            f = df * pow(k, n) - s * k - df + s;
            ff = df * n * pow(k, n - 1) - s;
            k = k - f / ff;
        } while (fabs(f) > 1e-10);
        return k;
    }
    template <class T>
    void geometric(T start, T end, double f_int)
    {
        double len = *(end - 1) - *start;
        int n = end - start;
        double k = solveK(n - 1, f_int, len);
        double kp = 1;
        for (int i = 0;i < n - 2;i++)
        {
            *(end - i - 2) = *(end - i - 1) - kp * f_int;
            kp *= k;
        }
    }
    template <class T>
    void geometric(double f_int, T start, T end)
    {
        double len = *(end - 1) - *start;
        int n = end - start;
        double k = solveK(n - 1, f_int, len);
        double kp = 1;
        for (int i = 0;i < n - 2;i++)
        {
            *(start + i + 1) = *(start + i) + kp * f_int;
            kp *= k;
        }
    }



} // namespace linesplit


#endif