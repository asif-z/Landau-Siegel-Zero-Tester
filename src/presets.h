#ifndef PRESETS_H
#define PRESETS_H

enum Preset {
    smallX1,
    mediumX1,
    largeX1,
    smallX2,
    mediumX2,
    largeX2
  };

void initializeLambda(enum Preset preset, arb_t lambda, arb_t phi, arb_t E, long prec);

#endif