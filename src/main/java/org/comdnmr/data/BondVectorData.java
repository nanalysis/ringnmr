/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author brucejohnson
 */
public class BondVectorData {

    double[][] vectors;
    List<String> names = new ArrayList<>();
    Map<String, Integer> indexMap = new HashMap<>();

    public BondVectorData(double[][] vectors, List<String> names, Map<String, Integer> indexMap) {
        this.vectors = vectors;
        this.names.addAll(names);
        this.indexMap.putAll(indexMap);
    }
    
    public Integer getIndex(String name) {
        return indexMap.get(name);
    }
    
    public double[][] getVectors() {
        return vectors;
    }

}
