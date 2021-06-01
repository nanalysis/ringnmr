package org.comdnmr.data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.nmrfx.chemistry.relax.RelaxationValues;

/**
 *
 * @author brucejohnson
 */
public class RelaxSet implements ValueSet {

    Map<String, RelaxationValues> map = new HashMap<>();
    List<RelaxationValues> relaxValueList = new ArrayList<>();

    public void addValue(RelaxationValues relaxValues) {
        relaxValueList.add(relaxValues);
        map.put(String.valueOf(relaxValues.getAtom().getResidueNumber()), relaxValues);
    }

    public List<RelaxationValues> getValues() {
        return relaxValueList;
    }

    public Set<String> getResidueNumberStrs() {
        return map.keySet();
    }

}
