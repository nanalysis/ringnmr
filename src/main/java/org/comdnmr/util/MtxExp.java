/*
 * CoMD/NMR Software : A Program for Analyzing NMR Dynamics Data
 * Copyright (C) 2018-2019 Bruce A Johnson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.util;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.LUDecomposition;

/**
 *
 * @author Martha Beckwith
 *
 * computes the exponential of a matrix, M: exp[M]
 *
 * based on the Matlab function c8mat_expm1() from
 * https://people.sc.fsu.edu/~jburkardt/m_src/matrix_exponential/r8mat_expm1.m
 *
 */
public class MtxExp {

    public static double[][] matrixExp(double[][] A1) {
        RealMatrix A = new Array2DRowRealMatrix(A1);
        RealMatrix E = matrixExp(A);
        double[][] E1 = new double[E.getRowDimension()][E.getColumnDimension()];
        for (int i = 0; i < E.getRowDimension(); i++) {
            for (int j = 0; j < E.getColumnDimension(); j++) {
                E1[i][j] = E.getEntry(i, j);
            }
        }
        return E1;
    }

    public static RealMatrix matrixExp(RealMatrix A) {
        //RealMatrix A = A1.copy();
        int e = (int) (Math.log(A.getNorm()) / Math.log(2));  // log base 2 of the infinity norm of A.
        int s = Math.max(0, e + 1);

        A = A.scalarMultiply(1 / Math.pow(2, s));

        RealMatrix X = A.copy();
        double c = 0.5;
        RealMatrix E = MatrixUtils.createRealIdentityMatrix(A.getRowDimension()).add(A.scalarMultiply(c));
        RealMatrix D = MatrixUtils.createRealIdentityMatrix(A.getRowDimension()).subtract(A.scalarMultiply(c));
        double q = 6;
        boolean p = true;

        for (int k = 2; k <= q; k++) {
            c = c * (q - k + 1) / (k * (2 * q - k + 1));
            X = A.multiply(X);
            RealMatrix cX = X.scalarMultiply(c);
            E = E.add(cX);

            if (p == true) {
                D = D.add(cX);
            } else {
                D = D.subtract(cX);
            }
            p = !p;
        }

        LUDecomposition LUD = new LUDecomposition(D);
        E = LUD.getSolver().solve(E);

        for (int k = 1; k <= s; k++) {
            E = E.multiply(E);
        }
        return E;
    }
}
