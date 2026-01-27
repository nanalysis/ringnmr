
package org.comdnmr.modelfree;

import org.nmrfx.chemistry.relax.OrderPar;

public record SimonsModelFitResult(
    OrderPar orderPar,
    double[][] replicateData,
    Double validationValue,
    SimonsScore[] scores
) {}
