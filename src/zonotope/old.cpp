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

 */

