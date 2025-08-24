#include <flint/arb.h>
#include "../src/presets.h"

void initializeLambda(enum Preset preset, arb_t lambda, arb_t phi, arb_t E, long prec) {
    // c = 0.1
    //small X (pi(X)~10000, X~104729)
    if (preset==smallX1) {
        arb_set_str(lambda, "1.4526", prec);
        arb_set_str(phi, "0.22910585", prec);
        arb_set_str(E, "1.0627878", prec);
    }
    //large X (pi(X)~ 540,000, X~8003537)
    else if (preset==mediumX1) {
        arb_set_str(lambda, "1.1122", prec);
        arb_set_str(phi, "0.23344604", prec);
        arb_set_str(E, "1.09531399", prec);
    }
    //extremely large X (pi(X)~ 5million, X~86028121)
    else if (preset==largeX1) {
        arb_set_str(lambda, "1", prec);
        arb_set_str(phi, "0.235032993", prec);
        arb_set_str(E, "1.10741261", prec);
    }

    // c = 0.2
    //small X (pi(X)~10000, X~104729)
    else if (preset==smallX2) {
        arb_set_str(lambda, "1.87", prec);
        arb_set_str(phi, "0.224", prec);
        arb_set_str(E, "1.02", prec);
    }
    //large X (pi(X)~ 540,000, X~8003537)
    else if (preset==mediumX2) {
        arb_set_str(lambda, "1.45", prec);
        arb_set_str(phi, "0.229", prec);
        arb_set_str(E, "1.061", prec);
    }
    //extremely large X (pi(X)~ 5million, X~86028121)
    else if (preset==largeX2) {
        arb_set_str(lambda, "1.3", prec);
        arb_set_str(phi, "0.231", prec);
        arb_set_str(E, "1.08", prec);
    }
}