package org.comdnmr.data;

import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.optim.PointValuePair;

public record CMAESFitResult(PointValuePair result, CMAESOptimizer optimizer) {}
