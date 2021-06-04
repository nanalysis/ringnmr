/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.data;

import java.util.Set;

/**
 *
 * @author brucejohnson
 */
public interface ValueSet {

    public String getName();
    public Set<String> getResidueNumberStrs();

}
