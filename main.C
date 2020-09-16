#include <iostream>
#include <vector>
#include "linesplit.H"
using namespace std;
using namespace linesplit;
int main()
{
    vector<double> x(51), y(51);
    vector<vector<double> > phi(51);
    for (int i = 0; i < 51; i++)
        phi[i].resize(51);

    x[0] = -50;
    x[15] = 0;
    x[35] = 1;
    x[50] = 51;

    linear(x.begin(), x.begin() + 16, 1.0 / 20);
    linear(x.begin() + 15, x.begin() + 36);
    linear(1.0 / 20, x.begin() + 35, x.begin() + 51);
    
    y[0] = 0;
    y[50] = 50;
    linear(0.006, y.begin(), y.end());

    for (int i = 0;i < 51;i++)
    {
        cout << i << " " << y[i] << endl;
    }



    for (int i = 0; i < 51; i++)
        for (int j = 0; j < 51; j++) {
            phi[i][j] = i * 51 + j;
        }
    /*
  for (int i = 0; i < 51; i++)
      for (int j = 0; j < 51; j++)
      {
          cout<<phi[i][j]<<endl;
      }
  */
    int t;
    cin >> t;
    return 0;
}
