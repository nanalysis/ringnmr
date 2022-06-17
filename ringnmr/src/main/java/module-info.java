module org.comdnmr {
    exports org.comdnmr.data;
    exports org.comdnmr.fit;
    exports org.comdnmr.modelfree;
    exports org.comdnmr.util;
    exports org.comdnmr.eqnfit;
    exports org.comdnmr.modelfree.models;
    requires commons.math3;
    requires smile.data;
    requires smile.core;
    requires smile.interpolation;
    requires smile.math;
    requires java.logging;
    requires java.desktop;
    requires com.google.common;
    requires jython.slim;
    requires java.prefs;
    requires ojalgo;
    requires org.yaml.snakeyaml;
    requires org.nmrfx.core;
    requires ejml.fat;
    requires javafx.base;
    requires javafx.graphics;
}
