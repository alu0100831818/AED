#include <iostream>

using namespace std;

#include "matrix_t.hpp"
#include "equation_sys_t.hpp"


int main(void)
{
    AED::equation_sys_t equation_system;
    AED::matrix_t m;
    cin >> equation_system;  //equation_system.from_stream(cin); esto queda mal...
    equation_system.to_stream(cout);
    equation_system.solut_inverse(m);
    return 0;
}