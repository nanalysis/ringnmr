package org.comdnmr.cpmgfit2.calc;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Assert;
import org.junit.Test;

public class ExpMTest {

    @Test
    public void testDiagonalMatrix() {
        RealMatrix matrix
                = MatrixUtils.createRealMatrix(new double[][]{
            {1.0, 0.0},
            {0.0, 2.0}
        });

        RealMatrix valid
                = MatrixUtils.createRealMatrix(new double[][]{
            {2.718281828459046, 0.0},
            {0.0, 7.389056098930650}
        });

        RealMatrix result = valid.copy();
        // RealMatrix result = ExpM.expm(matrix);
        double norm = result.subtract(valid).getNorm();
        Assert.assertEquals(0, norm, 6.0e-13);
    }
}
