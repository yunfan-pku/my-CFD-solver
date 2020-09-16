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

} // namespace linesplit


#endif