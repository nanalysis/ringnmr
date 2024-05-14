package org.comdnmr.ringnmr;

import javafx.scene.control.Menu;
import javafx.scene.control.MenuItem;
import javafx.stage.Stage;
import javafx.stage.StageStyle;
import org.comdnmr.gui.PyController;
import org.nmrfx.plugin.api.EntryPoint;
import org.nmrfx.plugin.api.NMRFxPlugin;
import org.nmrfx.plugin.api.PluginFunction;

import java.util.Set;
import java.util.function.Function;

public class RINGNmrPlugin implements NMRFxPlugin {
    public static PyController ringNMRController;
    Function<String, String> nmrfxFunction;

    @Override
    public Set<EntryPoint> getSupportedEntryPoints() {
        return Set.of(EntryPoint.MENU_PLUGINS);
    }

    @Override
    public void registerOnEntryPoint(EntryPoint entryPoint, Object object) {
        if (entryPoint != EntryPoint.MENU_PLUGINS) {
            throw new IllegalArgumentException("Only " + EntryPoint.MENU_PLUGINS + " is supported by RingNMR");
        }
        Menu menu;
        if (object instanceof PluginFunction pluginListener) {
            menu = (Menu) pluginListener.guiObject();
            nmrfxFunction = pluginListener.pluginFunction();
        } else if ((object instanceof Menu)) {
            menu = (Menu) object;
        } else {
            throw new IllegalArgumentException("Expected a menu, but received " + (object == null ? "null" : object.getClass().getName()) + " instead");
        }

        menu.getItems().add(createDynamicsMenu());
    }

    private Menu createDynamicsMenu() {
        Menu dynamicsMenu = new Menu("Dynamics");
        MenuItem ringNMRMenuItem = new MenuItem("Show RINGNMRGui");
        ringNMRMenuItem.setOnAction(e -> showRING());
        dynamicsMenu.getItems().addAll(ringNMRMenuItem);
        return dynamicsMenu;
    }

    private void showRING() {
        if (ringNMRController == null) {
            Stage stage = new Stage(StageStyle.DECORATED);
            ringNMRController = PyController.create(stage);
            ringNMRController.setNMRFxFunction(nmrfxFunction);
        }
        Stage stage = ringNMRController.getStage();
        stage.toFront();
        stage.show();
    }
}
