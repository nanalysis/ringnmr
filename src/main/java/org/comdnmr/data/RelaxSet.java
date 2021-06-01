
package org.comdnmr.data;

import java.util.ArrayList;
import java.util.List;
import org.nmrfx.chemistry.relax.RelaxationValues;

/**
 *
 * @author brucejohnson
 */
public class RelaxSet implements ValueSet {
    List<RelaxationValues> relaxValueList = new ArrayList<>();
    
    public void addValue(RelaxationValues relaxValues) {
        relaxValueList.add(relaxValues);
    }
    
    public List<RelaxationValues> getValues() {
        return relaxValueList;
    }
}
