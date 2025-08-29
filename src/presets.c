#include <flint/arb.h>
#include "../src/presets.h"

void initializeLambda(enum Preset preset, arb_t lambda, arb_t phi, arb_t E, long prec) {
    // c = 0.1
    //small X (pi(X)~10000, X~104729)
    if (preset==smallX1) {
        arb_set_str(lambda, "1.4", prec);
        arb_set_str(phi, "0.22946", prec);
        arb_set_str(E, "1.1756", prec);
    }
    //large X (pi(X)~ 540,000, X~8003537)
    else if (preset==mediumX1) {
        arb_set_str(lambda, "1.1", prec);
        arb_set_str(phi, "0.23362", prec);
        arb_set_str(E, "1.1952", prec);
    }
    //extremely large X (pi(X)~ 5million, X~86028121)
    else if (preset==largeX1) {
        arb_set_str(lambda, "1", prec);
        arb_set_str(phi, "0.23504", prec);
        arb_set_str(E, "1.2019", prec);
    }

    // c = 0.2
    //small X (pi(X)~10000, X~104729)
    else if (preset==smallX2) {
        arb_set_str(lambda, "1.8", prec);
        arb_set_str(phi, "0.2241", prec);
        arb_set_str(E, "1.1505", prec);
    }
    //large X (pi(X)~ 540,000, X~8003537)
    else if (preset==mediumX2) {
        arb_set_str(lambda, "1.4", prec);
        arb_set_str(phi, "0.22946", prec);
        arb_set_str(E, "1.1756", prec);
    }
    //extremely large X (pi(X)~ 5million, X~86028121)
    else if (preset==largeX2) {
        arb_set_str(lambda, "1.3", prec);
        arb_set_str(phi, "0.23083", prec);
        arb_set_str(E, "1.1821", prec);
    }
}