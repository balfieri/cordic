// come back to this later
// these are likely just gonna be routines at the std:: level
//
#if 0
    // machine learning
    T    tanh_backprop( const T& x, const T& x_backprop ) const;        // (1-x^2) * x_backprop
    T    sigmoid( const T& x ) const;                                   // 1/(1 + exp(-x)) = exp(x)/(exp(x) + 1)
    T    sigmoid_backprop( const T& x, const T& x_backprop ) const;     // x * (1-x) * x_backprop
    T    relu( const T& x ) const;                                      // (x > 0) x : 0
    T    relu_backprop( const T& x, const T& x_backprop ) const;        // (x > 0) x_backprop : 0
#endif
