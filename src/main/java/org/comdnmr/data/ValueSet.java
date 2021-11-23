/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.comdnmr.data;

import java.util.Set;
import org.nmrfx.chemistry.relax.ResonanceSource;

/**
 *
 * @author brucejohnson
 */
public interface ValueSet {

    public String getName();

    public Set<ResonanceSource> getDynamicsSources();

}
