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
package org.comdnmr.gui;

import javafx.application.Application;
import javafx.application.HostServices;
import javafx.application.Platform;
import javafx.scene.control.Alert;
import javafx.scene.control.ButtonType;
import javafx.stage.Stage;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;
import org.python.util.InteractiveInterpreter;

public class MainApp extends Application {

    public static Stage primaryStage;
    public static ScriptEngine engine;
    public static InteractiveInterpreter interpreter = new InteractiveInterpreter();
    public static PreferencesController preferencesController;
    public static ConsoleRedirect console;
    public static HostServices hostServices;
    public static final String APP_NAME = "RING NMR Dynamics";

    public static void main(String[] args) throws ScriptException {
        engine = new ScriptEngineManager().getEngineByName("jython");
        engine.put(ScriptEngine.ARGV, args);
        engine.eval("from exporter import Writer");
        launch(args);
    }

    @Override
    public void start(Stage stage) throws Exception {
        PyController controller = PyController.create(stage);
        Platform.setImplicitExit(true);
        hostServices = getHostServices();
        stage.setTitle(APP_NAME + " " + ""); // fixme getVersion

        Parameters parameters = getParameters();
        System.out.println(parameters.getRaw());

//        interpreter.exec("from pyproc import *\ninitLocal()\nfrom gscript import *\nnw=NMRFxWindowScripting()\nfrom dscript import *\nfrom pscript import *\nimport os");
//        interpreter.set("argv", parameters.getRaw());
//        interpreter.exec("parseArgs(argv)");
        // Dataset.addObserver(this);
//        if (defaultFont == null) {
//            loadFont();
//        }
    }

    public static InteractiveInterpreter getInterpreter() {
        return interpreter;
    }

    public static void setConsoleController(ConsoleRedirect controller) {
        console = controller;
    }

    public static void checkExit() {
        Alert alert = new Alert(Alert.AlertType.CONFIRMATION, "Exit program", ButtonType.CANCEL, ButtonType.YES);
        alert.showAndWait().ifPresent(response -> {
            if (response == ButtonType.YES) {
                Platform.exit();
                System.exit(0);
            }
        });
    }
}
