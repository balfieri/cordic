// later, come back and define two classes:
// 
// Rand_Uniform - most of the routines below could be used for this
// Rand_Dists - would like to find a template somewhere else that's already written
//
#if 0
    // random numbers
    T    rand_seed( uint64_t seed0, uint64_t seed1=0xdeadbeefbabecafe );  // set random seed(s)
    T    rand_raw32_mwc( void );                                        // 32 bits from Frank Marsaglia's MWC RNG     (good fast&small)
    T    rand_raw32_jsf( void );                                        // 32 bits from Robert Jenkins' JSF RNG       (better f&s)
    T    rand_raw64_jsf( void );                                        // 64 bits from Robert Jenkins' JSF RNG       (better f&s)
    T    rand_raw16_sfc( void );                                        // 16 bits from Chris Doty-Humphrey's SFC RNG (best   f&s)
    T    rand_raw32_sfc( void );                                        // 32 bits from Chris Doty-Humphrey's SFC RNG (best   f&s)
    T    rand_raw64_sfc( void );                                        // 64 bits from Chris Doty-Humphrey's SFC RNG (best   f&s)
    T    rand_raw128_hc128( void );                                     //128 bits from eSTREAM's HC128 RNG (larger, but crypto-secure)
    T    rand_raw256_hc256( void );                                     //256 bits from eSTREAM's HC256 RNG (larger, but crypto-secure)
    static void rand_raw32_fn_set( T (*raw32_fn)(void) );               // set default raw32() function to use for following routines:
    static void rand_raw64_fn_set( T (*raw64_fn)(void) );               // set default raw64() function to use for following routines:
    T    rand_uniform( void );                                          // return uniform random in range [0.0, 1.0)  (1.0 excluded)
    T    rand_gaussian( const T& mu, const T& std );                    // return gaussian random with mu and std using uniform_fn
#endif    
