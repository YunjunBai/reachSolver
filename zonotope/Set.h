/* 
 Description:
    Generates a linear system object according to the following first
    order differential equations:
       x' = A x + B u + c
       y = C x + D u + k
 Inputs:
    name - name of system
    A - state matrix
    B - input matrix
    c - constant input
    C - output matrix
    D - throughput matrix
    k - output offset
 Outputs:
     obj - Generated Object

 Reference:
   CORA ../contSet/contSet/contSet.m 
 */
#include <cstddef>

namespace reachSolver {

class Set{
  protected:
    int id;
    size_t dimension;
  public:
    Set();
    ~Set(){};
    void display();
};

}