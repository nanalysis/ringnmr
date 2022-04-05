package org.comdnmr.data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.nmrfx.chemistry.relax.RelaxationValues;
import org.nmrfx.chemistry.relax.ResonanceSource;

/**
 *
 * @author brucejohnson
 */
public class RelaxSet implements ValueSet {

    final String name;

    public RelaxSet(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    Map<ResonanceSource, RelaxationValues> map = new HashMap<>();
    List<RelaxationValues> relaxValueList = new ArrayList<>();

    public void addValue(RelaxationValues relaxValues) {
        relaxValueList.add(relaxValues);
        map.put(relaxValues.getResonanceSource(), relaxValues);
    }

    public List<RelaxationValues> getValues() {
        return relaxValueList;
    }

    public RelaxationValues getValue(ResonanceSource resonanceSource) {
        return map.get(resonanceSource);
    }

    @Override
    public Set<ResonanceSource> getDynamicsSources() {
        return map.keySet();
    }

}
